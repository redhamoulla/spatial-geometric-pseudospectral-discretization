import unittest
import numpy as np
from gpps import MixedDiscretization,ProjectedGauss,transmission_line,shallow_water
from gpps.systems import PHSystem,interconnect
from gpps.time import integrate,exact_linear_flow


class QuarticOscillator:
    degree=4
    def value(self,x): return .5*np.sum(x*x)+.25*x[0]**4
    def gradient(self,x): return np.array([x[0]+x[0]**3,x[1]])
    def hessian(self,x): return np.diag([1+3*x[0]**2,1])


def oscillator():
    return PHSystem(np.array([[0.,1.],[-1.,0.]]),np.zeros((2,2)),np.zeros((2,0)),np.zeros((0,0)),QuarticOscillator())


class TimeTests(unittest.TestCase):
    def test_quadratic_energy_and_time_orders(self):
        m=MixedDiscretization(4)
        sys=transmission_line(m)
        x=np.r_[np.sin(np.pi*m.z/2),np.zeros(4)]
        target=exact_linear_flow(sys,x,1)
        for s in (1,2):
            errors=[]
            for steps in (10,20,40):
                path,rows=integrate(sys,x,1/steps,steps,ProjectedGauss(s,s))
                errors.append(np.linalg.norm(path[-1]-target))
                self.assertLess(max(abs(r["balance"]) for r in rows),1e-11)
            self.assertGreater(np.log2(errors[-2]/errors[-1]),2*s-.25)

    def test_polynomial_exactness_and_underintegration_negative_control(self):
        sys=oscillator()
        x=np.array([.8,.7])
        low=ProjectedGauss(1,1).step(sys,x,0,.2)
        exact=ProjectedGauss(1,2).step(sys,x,0,.2)
        self.assertGreater(abs(low.balance_residual),1e-6)
        self.assertLess(abs(exact.balance_residual),1e-12)
        self.assertAlmostEqual(low.balance_residual,low.chain_defect+low.solve_power_defect,places=13)

    def test_symplecticity_is_not_energy_preservation(self):
        sys=oscillator()
        x=np.array([.8,.7])
        midpoint=ProjectedGauss(1,1).step(sys,x,0,.2,tangent=True)
        avf=ProjectedGauss(1,2).step(sys,x,0,.2,tangent=True)
        def defect(step): return np.linalg.norm(step.tangent.T@sys.J@step.tangent-sys.J)
        self.assertLess(defect(midpoint),1e-12)
        self.assertGreater(defect(avf),1e-6)

    def test_analytic_tangent_against_finite_difference(self):
        sys=oscillator();x=np.array([.8,.7]);method=ProjectedGauss(2,4)
        tangent=method.step(sys,x,0,.1,tangent=True).tangent
        eps=1e-6
        finite=np.column_stack([(method.step(sys,x+eps*e,0,.1).state-method.step(sys,x-eps*e,0,.1).state)/(2*eps) for e in np.eye(2)])
        np.testing.assert_allclose(tangent,finite,rtol=1e-7,atol=1e-9)

    def test_nonlinear_forced_dissipative_balance(self):
        m=MixedDiscretization(5)
        sys=transmission_line(m,1,1,resistance=.2,quartic=1)
        x=np.r_[.3*np.sin(np.pi*m.z/2),np.zeros(5)]
        u=lambda t:np.array([.2*np.sin(t),.1*np.cos(2*t)])
        _,rows=integrate(sys,x,.02,30,ProjectedGauss(2,4),u)
        self.assertLess(max(abs(r["balance"]) for r in rows),2e-11)
        self.assertGreater(rows[-1]["dissipation"],0)

    def test_shallow_water_rest_and_cubic_closure(self):
        m=MixedDiscretization(6);sys=shallow_water(m)
        rest=np.r_[np.ones(6),np.zeros(6)]
        u=lambda t:np.array([1.,0.])
        self.assertLess(np.linalg.norm(sys.rhs(rest,0,u)),1e-11)
        x=rest.copy();x[:6]+=.02*np.sin(np.pi*m.z/2)
        r=ProjectedGauss(2,3).step(sys,x,0,.02,u)
        self.assertLess(abs(r.balance_residual),1e-12)

    def test_interconnection_cancels_internal_power(self):
        meshes=[MixedDiscretization(4),MixedDiscretization(3)]
        parts=[transmission_line(m,1,1) for m in meshes]
        K=np.zeros((4,4));K[1,2]=1;K[2,1]=-1
        sys,ports=interconnect(parts,K)
        g=np.random.default_rng(3).normal(size=sys.dimension)
        u,y=ports(g,np.zeros(4))
        self.assertAlmostEqual(u@y,0,places=11)
        self.assertAlmostEqual(u[1]*y[1]+u[2]*y[2],0,places=11)
        self.assertLess(np.linalg.norm(sys.J+sys.J.T),1e-10)

    def test_nonlinear_solve_failure_is_not_silently_accepted(self):
        with self.assertRaises(RuntimeError):
            ProjectedGauss(2,4,tolerance=1e-15,max_iterations=1).step(oscillator(),np.array([2.,1.]),0,1.)

    def test_invalid_quadrature_rejected(self):
        with self.assertRaises(ValueError):
            ProjectedGauss(3,2)


if __name__ == "__main__":
    unittest.main()
