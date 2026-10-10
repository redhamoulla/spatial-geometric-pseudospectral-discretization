"""Polynomial coordinates in a Legendre basis; no monomial Vandermonde."""
import numpy as np
from numpy.polynomial.legendre import Legendre, legvander


def nodes(n, family="gauss", interval=(0.0, 1.0), seed=20261010):
    if n < 1:
        raise ValueError("At least one point is required")
    a, b = interval
    if not b > a:
        raise ValueError("The interval must have positive length")
    if family == "gauss":
        x = np.polynomial.legendre.leggauss(n)[0]
    elif family == "lobatto":
        x = (np.array([0.0]) if n == 1 else
             np.r_[-1.0, Legendre.basis(n - 1).deriv().roots(), 1.0])
    elif family == "uniform":
        x = np.linspace(-1, 1, n + 2)[1:-1]
    elif family == "jittered":
        rng = np.random.default_rng(seed + n)
        x = np.linspace(-1, 1, n + 2)[1:-1] + rng.uniform(-.35, .35, n) * 2 / (n + 1)
    else:
        raise ValueError(f"Unknown node family: {family}")
    return a + (b - a) * (x + 1) / 2


def gauss_rule(n, interval=(0.0, 1.0)):
    a, b = interval
    x, w = np.polynomial.legendre.leggauss(n)
    return a + (b - a) * (x + 1) / 2, w * (b - a) / 2


class LagrangeBasis:
    def __init__(self, points, interval=(0.0, 1.0)):
        self.points = np.asarray(points, dtype=float)
        self.a, self.b = interval
        if self.points.ndim != 1 or not len(self.points) or not self.b > self.a:
            raise ValueError("Invalid nodes or interval")
        if not np.all(np.isfinite(self.points)) or len(np.unique(self.points)) != len(self.points):
            raise ValueError("Nodes must be finite and distinct")
        if np.min(self.points) < self.a or np.max(self.points) > self.b:
            raise ValueError("Nodes must lie in the interval")
        self.n = len(self.points)
        v = legvander(self._map(self.points), self.n - 1)
        self.condition = float(np.linalg.cond(v))
        self.coefficients = np.linalg.solve(v, np.eye(self.n))

    def _map(self, z):
        return 2 * (np.asarray(z) - self.a) / (self.b - self.a) - 1

    def values(self, z):
        return legvander(self._map(np.atleast_1d(z)), self.n - 1) @ self.coefficients

    def derivatives(self, z):
        x = self._map(np.atleast_1d(z))
        d = np.column_stack([Legendre.basis(j).deriv()(x) for j in range(self.n)])
        return (2 / (self.b - self.a)) * d @ self.coefficients


def temporal_basis(degree, c):
    """Orthonormal shifted Legendre P_j and I_j(c)=integral_0^c P_j."""
    c = np.asarray(c)
    p, integ = [], []
    for j in range(degree):
        poly = Legendre.basis(j) * np.sqrt(2 * j + 1)
        primitive = poly.integ()
        p.append(poly(2 * c - 1))
        integ.append((primitive(2 * c - 1) - primitive(-1)) / 2)
    return np.column_stack(p), np.column_stack(integ)
