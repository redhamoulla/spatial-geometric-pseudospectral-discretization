"""Energy closures and constant-structure port-Hamiltonian models."""
from dataclasses import dataclass
import numpy as np
from .polynomials import gauss_rule


@dataclass
class PHSystem:
    J: np.ndarray
    R: np.ndarray
    B: np.ndarray
    S: np.ndarray
    hamiltonian: object
    name: str = "pH system"

    def __post_init__(self):
        self.J, self.R, self.B, self.S = [np.asarray(a, dtype=float) for a in (self.J,self.R,self.B,self.S)]
        d = self.J.shape[0]
        if self.J.shape != (d,d) or self.R.shape != (d,d) or self.B.shape[0] != d:
            raise ValueError("Incompatible system matrix dimensions")
        m = self.B.shape[1]
        if self.S.shape != (m,m):
            raise ValueError("Invalid feedthrough shape")
        tol = 2e-9 * max(1.0, np.linalg.norm(self.J), np.linalg.norm(self.R))
        if np.linalg.norm(self.J+self.J.T) > tol or np.linalg.norm(self.S+self.S.T) > tol:
            raise ValueError("J and S must be skew-symmetric")
        if np.linalg.norm(self.R-self.R.T) > tol or np.min(np.linalg.eigvalsh(self.R)) < -tol:
            raise ValueError("R must be symmetric positive semidefinite")

    @property
    def dimension(self):
        return self.J.shape[0]

    def input(self, u, t):
        value = np.zeros(self.B.shape[1]) if u is None else np.asarray(u(t),dtype=float)
        if value.shape != (self.B.shape[1],) or not np.all(np.isfinite(value)):
            raise ValueError("Input callback returned an invalid vector")
        return value

    def rhs(self, x, t=0.0, u=None):
        return (self.J-self.R) @ self.hamiltonian.gradient(x) + self.B @ self.input(u,t)


class QuadraticEnergy:
    degree = 2

    def __init__(self, Q):
        self.Q = np.asarray(Q)

    def value(self, x):
        return float(.5 * x @ self.Q @ x)

    def gradient(self, x):
        return self.Q @ x

    def hessian(self, x):
        return self.Q


class FieldEnergy:
    """Exactly integrated polynomial closure (quartic line or cubic shallow water)."""
    def __init__(self, mesh, kind="quartic", alpha=1.0, gravity=1.0):
        self.n, self.kind, self.alpha, self.gravity = mesh.n, kind, alpha, gravity
        self.degree = 4 if kind == "quartic" else 3
        count = max(mesh.n+1, (self.degree * (mesh.n-1) + 2)//2)
        z, self.w = gauss_rule(count, mesh.interval)
        self.V = mesh.flow.values(z)

    def _fields(self, x):
        return self.V @ x[:self.n], self.V @ x[self.n:]

    def value(self, x):
        q, p = self._fields(x)
        density = (.5*(q*q+p*p)+.25*self.alpha*q**4 if self.kind == "quartic"
                   else .5*self.gravity*q*q+.5*q*p*p)
        return float(self.w @ density)

    def gradient(self, x):
        q, p = self._fields(x)
        eq, ep = ((q+self.alpha*q**3,p) if self.kind == "quartic"
                  else (self.gravity*q+.5*p*p,q*p))
        return np.r_[self.V.T@(self.w*eq), self.V.T@(self.w*ep)]

    def hessian(self, x):
        q, p = self._fields(x)
        a,b,c = ((1+3*self.alpha*q*q,np.zeros_like(q),np.ones_like(q))
                 if self.kind == "quartic" else (np.full_like(q,self.gravity),p,q))
        def mass(v):
            return self.V.T @ ((self.w*v)[:,None]*self.V)
        return np.block([[mass(a),mass(b)],[mass(b),mass(c)]])


def transmission_line(mesh, capacitance=3.0, inductance=2.0, resistance=0.0, quartic=0.0):
    if capacitance <= 0 or inductance <= 0 or resistance < 0 or quartic < 0:
        raise ValueError("Positive constitutive parameters and nonnegative damping are required")
    n, p = mesh.n, mesh.ports
    z = np.zeros((n,n))
    Q = np.block([[mesh.mass/capacitance,z],[z,mesh.mass/inductance]])
    energy = QuadraticEnergy(Q) if quartic == 0 else FieldEnergy(mesh,"quartic",alpha=quartic)
    if quartic and (capacitance != 1 or inductance != 1):
        raise ValueError("The normalized quartic line requires C=L_ind=1")
    R = np.block([[z,z],[z,resistance*np.linalg.solve(mesh.mass,np.eye(n))]])
    return PHSystem(p.J,R,p.B,p.S,energy,"transmission line")


def shallow_water(mesh, gravity=1.0):
    if gravity <= 0:
        raise ValueError("gravity must be positive")
    p = mesh.ports
    return PHSystem(p.J,np.zeros_like(p.J),p.B,p.S,
                    FieldEnergy(mesh,"shallow_water",gravity=gravity),"shallow water")


def block_diag(*arrays):
    out = np.zeros((sum(a.shape[0] for a in arrays),sum(a.shape[1] for a in arrays)))
    i=j=0
    for a in arrays:
        out[i:i+a.shape[0],j:j+a.shape[1]]=a
        i,j=i+a.shape[0],j+a.shape[1]
    return out


class SumEnergy:
    def __init__(self, systems):
        self.energies = [s.hamiltonian for s in systems]
        self.cuts = np.cumsum([0]+[s.dimension for s in systems])
        self.degree = max(e.degree for e in self.energies)

    def value(self,x):
        return sum(e.value(x[a:b]) for e,a,b in zip(self.energies,self.cuts[:-1],self.cuts[1:]))

    def gradient(self,x):
        return np.concatenate([e.gradient(x[a:b]) for e,a,b in zip(self.energies,self.cuts[:-1],self.cuts[1:])])

    def hessian(self,x):
        return block_diag(*[e.hessian(x[a:b]) for e,a,b in zip(self.energies,self.cuts[:-1],self.cuts[1:])])


def interconnect(systems, K):
    """Eliminate u=K y+v, retaining the boundary feedthrough S."""
    J,R,B,S = [block_diag(*[getattr(s,key) for s in systems]) for key in ("J","R","B","S")]
    K = np.asarray(K)
    if K.shape != S.shape or np.linalg.norm(K+K.T) > 1e-12:
        raise ValueError("The interconnection matrix must be skew-symmetric")
    A = np.eye(len(K))-K@S
    if np.linalg.cond(A) > 1e12:
        raise ValueError("Singular interconnection: an implicit DAE treatment is required")
    F = np.linalg.solve(A,np.eye(len(K)))
    system = PHSystem(J+B@F@K@B.T,R,B@F,S@F,SumEnergy(systems),"interconnected lines")
    def internal_ports(e,v):
        u = F@(K@B.T@e+v)
        return u,B.T@e+S@u
    return system,internal_ports
