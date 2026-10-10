"""Mixed P_N / P_(N-1) construction, quotient, and boundary causality.

Sign convention: q_t=-d_z e_p, p_t=-d_z e_q.
u=(e_q(a),e_p(b)); y=(e_p(a),-e_q(b)); Hdot=u.T@y.
"""
from dataclasses import dataclass
import numpy as np
from .polynomials import LagrangeBasis, nodes, gauss_rule


@dataclass
class PortMatrices:
    J: np.ndarray
    B: np.ndarray
    C: np.ndarray
    S: np.ndarray


class MixedDiscretization:
    def __init__(self, n, interval=(0.0, 1.0), family="gauss",
                 flow_nodes=None, effort_nodes=None):
        if n < 1:
            raise ValueError("n must be positive")
        self.n, self.interval = n, interval
        a, b = interval
        self.z = nodes(n, family, interval) if flow_nodes is None else np.asarray(flow_nodes)
        self.zeta = nodes(n + 1, family, interval) if effort_nodes is None else np.asarray(effort_nodes)
        if len(self.z) != n or len(self.zeta) != n + 1:
            raise ValueError("Need n flow nodes and n+1 effort nodes")
        self.flow = LagrangeBasis(self.z, interval)
        self.effort = LagrangeBasis(self.zeta, interval)
        zq, wq = gauss_rule(n + 1, interval)
        phi, psi = self.effort.values(zq), self.flow.values(zq)
        self.M = phi.T @ (wq[:, None] * psi)
        self.D = self.effort.derivatives(self.z)
        self.mass = psi.T @ (wq[:, None] * psi)
        self.effort_mass = phi.T @ (wq[:, None] * phi)
        self.trace = self.effort.values([a, b])
        self.boundary = np.outer(self.trace[1], self.trace[1]) - np.outer(self.trace[0], self.trace[0])
        self.left_completion = np.vstack((self.M.T, self.trace[0]))
        self.right_completion = np.vstack((self.M.T, self.trace[1]))
        self.ports = self._ports()

    def _ports(self):
        n = self.n
        # Full efforts reconstructed from reduced effort and one boundary trace.
        rq = np.zeros((n + 1, 2 * n + 2))
        rp = np.zeros_like(rq)
        rq[:n, :n] = np.eye(n)
        rp[:n, n:2*n] = np.eye(n)
        rq[n, 2*n], rp[n, 2*n + 1] = 1, 1
        eq = np.linalg.solve(self.left_completion, rq)
        ep = np.linalg.solve(self.right_completion, rp)
        f = np.vstack((-self.D @ ep, -self.D @ eq))
        y = np.vstack((self.trace[0] @ ep, -self.trace[1] @ eq))
        # Keep raw matrices: do NOT repair skewness after assembly.
        return PortMatrices(f[:, :2*n], f[:, 2*n:], y[:, :2*n], y[:, 2*n:])

    def diagnostics(self):
        p = self.ports
        stokes = self.M @ self.D + self.D.T @ self.M.T - self.boundary
        extended = np.block([[p.J, p.B], [-p.C, -p.S]])
        scale = max(1.0, np.linalg.norm(extended))
        n = self.n
        zero = np.zeros_like(self.D)
        flows = np.block([[zero, -self.D], [-self.D, zero]])
        reduced = np.block([[self.M.T, zero], [zero, self.M.T]])
        traces = np.block([[self.trace, np.zeros_like(self.trace)],
                           [np.zeros_like(self.trace), self.trace]])
        image_matrix = np.vstack((flows, reduced, traces))
        return {
            "stokes_relative": float(np.linalg.norm(stokes) / max(1.0, np.linalg.norm(self.boundary))),
            "dirac_relative": float(np.linalg.norm(extended + extended.T) / scale),
            "pairing_rank": int(np.linalg.matrix_rank(self.M)),
            "dirac_image_rank": int(np.linalg.matrix_rank(image_matrix)),
            "condition_flow": self.flow.condition,
            "condition_effort": self.effort.condition,
            "condition_completion": float(max(np.linalg.cond(self.left_completion), np.linalg.cond(self.right_completion))),
        }


def tensor_diagnostics(space, time):
    """Two independent checks: boundary integration and d_1 d_0=0."""
    sx, st = space.M @ space.D, time.M @ time.D
    gx, gt = space.effort_mass, time.effort_mass
    k = np.kron(st, gx) + np.kron(gt, sx)
    bd = np.kron(time.boundary, gx) + np.kron(gt, space.boundary)
    d0 = np.vstack((np.kron(time.D, np.eye(space.n + 1)),
                    np.kron(np.eye(time.n + 1), space.D)))
    d1 = np.hstack((-np.kron(np.eye(time.n), space.D),
                    np.kron(time.D, np.eye(space.n))))
    return {"stokes_relative": float(np.linalg.norm(k+k.T-bd) / max(1.0, np.linalg.norm(bd))),
            "complex_relative": float(np.linalg.norm(d1@d0) / max(1.0, np.linalg.norm(d1)*np.linalg.norm(d0)))}
