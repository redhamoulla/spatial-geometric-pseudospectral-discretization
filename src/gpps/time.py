"""Continuous polynomial state / projected effort time integration.

For a constant canonical J without ports/damping this is HBVM(k,s).
For s=1 and exact integration it is AVF. For k=s it is Gauss collocation.
The implementation reports quadrature and nonlinear-solve defects separately.
"""
from dataclasses import dataclass
import numpy as np
from .polynomials import gauss_rule, temporal_basis


@dataclass
class Step:
    state: np.ndarray
    work: float
    dissipation: float
    balance_residual: float
    chain_defect: float
    solve_power_defect: float
    nonlinear_residual: float
    iterations: int
    stage_states: np.ndarray
    stage_efforts: np.ndarray
    stage_inputs: np.ndarray
    stage_outputs: np.ndarray
    tangent: object = None


class ProjectedGauss:
    def __init__(self, degree=1, quadrature=None, tolerance=3e-13, max_iterations=18):
        if degree < 1:
            raise ValueError("Polynomial degree must be positive")
        quadrature = degree if quadrature is None else quadrature
        if quadrature < degree:
            raise ValueError("Need quadrature >= degree for an exact temporal Gram matrix")
        if tolerance <= 0 or max_iterations < 1:
            raise ValueError("Invalid nonlinear solver settings")
        self.s,self.k,self.tolerance,self.max_iterations=degree,quadrature,tolerance,max_iterations
        self.c,self.w=gauss_rule(quadrature)
        self.P,self.I=temporal_basis(degree,self.c)

    def step(self, system, x, t, dt, u=None, tangent=False):
        if dt <= 0 or not np.isfinite(dt):
            raise ValueError("dt must be finite and positive")
        x=np.asarray(x,dtype=float)
        if x.shape != (system.dimension,) or not np.all(np.isfinite(x)):
            raise ValueError("Invalid initial state")
        s,d=self.s,len(x)
        L=system.J-system.R
        energy=system.hamiltonian
        inputs=np.array([system.input(u,t+dt*c) for c in self.c])

        def evaluate(a, need_jac=False):
            states=x+dt*self.I@a
            gradients=np.array([energy.gradient(z) for z in states])
            coeff=self.P.T@(self.w[:,None]*gradients)
            efforts=self.P@coeff
            velocities=efforts@L.T+inputs@system.B.T
            residual=a-self.P.T@(self.w[:,None]*velocities)
            if not need_jac:
                return residual,states,coeff,efforts
            hess=np.array([energy.hessian(z) for z in states])
            # dg_i/da_j = dt * sum_l b_l P_i(c_l) I_j(c_l) H''(X_l).
            dg=dt*np.einsum('k,ki,kj,kab->iajb',self.w,self.P,self.I,hess)
            jac=np.eye(s*d)-np.einsum('ab,ibjc->iajc',L,dg).reshape(s*d,s*d)
            return residual,states,coeff,efforts,jac,hess

        a=np.zeros((s,d))
        a[0]=system.rhs(x,t,u)
        iterations=0
        for iterations in range(self.max_iterations+1):
            residual,states,coeff,efforts=evaluate(a)
            norm=float(np.linalg.norm(residual.ravel(),ord=np.inf))
            target=self.tolerance*max(1.0,float(np.linalg.norm(a.ravel(),ord=np.inf)))
            if norm <= target:
                break
            if iterations == self.max_iterations:
                raise RuntimeError(f"Newton did not converge at t={t}: residual={norm:g}, target={target:g}")
            jac=evaluate(a,True)[4]
            delta=np.linalg.solve(jac,-residual.ravel()).reshape(s,d)
            scale=1.0
            for _ in range(16):
                trial=a+scale*delta
                new_norm=np.linalg.norm(evaluate(trial)[0].ravel(),ord=np.inf)
                if np.isfinite(new_norm) and (new_norm < norm or new_norm <= target):
                    a=trial
                    break
                scale*=.5
            else:
                raise RuntimeError(f"Newton line search failed at t={t}")
        x1=x+dt*a[0]
        outputs=efforts@system.B+inputs@system.S.T
        work=float(dt*np.einsum('k,ki,ki',self.w,inputs,outputs))
        diss=float(dt*np.einsum('k,ki,ij,kj',self.w,efforts,system.R,efforts))
        delta_h=energy.value(x1)-energy.value(x)
        chain=delta_h-float(dt*np.sum(coeff*a))
        solve_power=float(dt*np.sum(coeff*residual))
        derivative=None
        if tangent:
            _,_,_,_,jac,hess=evaluate(a,True)
            dgdx=np.einsum('k,ki,kab->iab',self.w,self.P,hess)
            rhs=np.einsum('ab,ibc->iac',L,dgdx).reshape(s*d,d)
            dadx=np.linalg.solve(jac,rhs).reshape(s,d,d)
            derivative=np.eye(d)+dt*dadx[0]
        return Step(x1,work,diss,float(delta_h-work+diss),float(chain),solve_power,norm,
                    iterations,states,efforts,inputs,outputs,derivative)


def integrate(system,x0,dt,steps,method,u=None,observer=None):
    """Return the energy ledger and nodal trajectory; no energy projection is applied."""
    if steps < 1:
        raise ValueError("steps must be positive")
    x=np.asarray(x0,dtype=float).copy()
    trajectory=[x.copy()]
    h0=system.hamiltonian.value(x)
    rows=[{"time":0.0,"energy":h0,"work":0.0,"dissipation":0.0,
           "balance":0.0,"step_balance":0.0,"chain_defect":0.0,
           "nonlinear_residual":0.0,"newton_iterations":0}]
    work=diss=0.0
    for n in range(steps):
        result=method.step(system,x,n*dt,dt,u)
        x=result.state
        work+=result.work
        diss+=result.dissipation
        energy=system.hamiltonian.value(x)
        rows.append({"time":(n+1)*dt,"energy":energy,"work":work,"dissipation":diss,
                     "balance":energy-h0-work+diss,"step_balance":result.balance_residual,
                     "chain_defect":result.chain_defect,"nonlinear_residual":result.nonlinear_residual,
                     "newton_iterations":result.iterations})
        trajectory.append(x.copy())
        if observer is not None:
            observer(n,result)
    return np.asarray(trajectory),rows


def exact_linear_flow(system,x,t):
    """Independent matrix-exponential reference via an energy-coordinate eigendecomposition.

Only closed undamped positive-definite quadratic systems are accepted.
"""
    if np.linalg.norm(system.R) != 0 or not hasattr(system.hamiltonian,"Q"):
        raise ValueError("Requires a conservative quadratic model")
    Q=system.hamiltonian.Q
    vals,vec=np.linalg.eigh(Q)
    if vals.min() <= 0:
        raise ValueError("Q must be positive definite")
    sqrtQ=(vec*np.sqrt(vals))@vec.T
    invsqrtQ=(vec*(1/np.sqrt(vals)))@vec.T
    K=sqrtQ@system.J@sqrtQ
    # i K is Hermitian up to measured assembly roundoff.
    freqs,U=np.linalg.eigh(1j*K)
    z=U@(np.exp(-1j*freqs*t)*(U.conj().T@(sqrtQ@x)))
    return np.real(invsqrtQ@z)
