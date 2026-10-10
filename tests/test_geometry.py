import unittest
import numpy as np
from gpps.spatial import MixedDiscretization,tensor_diagnostics
from gpps.systems import transmission_line


class GeometryTests(unittest.TestCase):
    def test_stokes_and_dirac_for_arbitrary_nodes(self):
        for family in ("gauss","lobatto","uniform","jittered"):
            for n in (1,2,5,8):
                with self.subTest(family=family,n=n):
                    m=MixedDiscretization(n,family=family)
                    d=m.diagnostics()
                    self.assertLess(d["stokes_relative"],2e-10)
                    self.assertLess(d["dirac_relative"],2e-10)
                    self.assertEqual(d["pairing_rank"],n)
                    self.assertEqual(d["dirac_image_rank"],2*n+2)

    def test_rank_one_effort_kernel_and_trace_completion(self):
        m=MixedDiscretization(7,family="jittered")
        _,_,v=np.linalg.svd(m.M.T,full_matrices=True)
        kernel=v[-1]
        self.assertLess(np.linalg.norm(m.M.T@kernel),1e-13)
        self.assertGreater(abs(m.trace[0]@kernel),1e-3)
        self.assertGreater(abs(m.trace[1]@kernel),1e-3)

    def test_direct_boundary_power_not_only_skewness(self):
        m=MixedDiscretization(6)
        rng=np.random.default_rng(51)
        g,u=rng.normal(size=12),rng.normal(size=2)
        eq=np.linalg.solve(m.left_completion,np.r_[g[:6],u[0]])
        ep=np.linalg.solve(m.right_completion,np.r_[g[6:],u[1]])
        flow=np.r_[-m.D@ep,-m.D@eq]
        y=np.array([m.trace[0]@ep,-m.trace[1]@eq])
        self.assertAlmostEqual(g@flow,u@y,places=10)
        np.testing.assert_allclose(y,m.ports.B.T@g+m.ports.S@u,atol=2e-11)

    def test_paper_spectrum_n8(self):
        m=MixedDiscretization(8,interval=(0,2),family="uniform")
        s=transmission_line(m)
        ev=np.linalg.eigvals(s.J@s.hamiltonian.Q)
        omega=np.sort(ev.imag[ev.imag>0])
        published=np.array([.32063745754047,.96191238151097,1.60321563740732,
                            2.24771134219675,2.94830702957369,3.98374558826630,
                            6.38606332778657,18.76179259703309])
        np.testing.assert_allclose(omega,published,rtol=2e-10,atol=1e-11)

    def test_space_time_boundary_and_complex(self):
        d=tensor_diagnostics(MixedDiscretization(5,family="uniform"),MixedDiscretization(3,family="jittered"))
        self.assertLess(d["stokes_relative"],1e-12)
        self.assertLess(d["complex_relative"],1e-14)

    def test_bad_nodes_rejected(self):
        with self.assertRaises(ValueError):
            MixedDiscretization(2,flow_nodes=[.1,.1])


if __name__ == "__main__":
    unittest.main()
