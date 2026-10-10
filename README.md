# Spatial geometric pseudospectral discretization

Theory and reproducible simulations of the spatial method introduced by
[Moulla–Lefèvre–Maschke (2012)](https://doi.org/10.1016/j.jcp.2011.10.008),
coupled to projected-effort time integration.

**[Read the complete manuscript](docs/manuscript.md)**: analysis of the paper,
spatial construction, proofs, time discretization, interconnection, simulations
and references, in one document with a table of contents.
Equations render directly on GitHub; figures are provided as SVG.

## Reproduce the simulations

Python 3.10 or later. NumPy handles the computations; ReportLab generates
the vector figures.

~~~bash
git clone https://github.com/redhamoulla/spatial-geometric-pseudospectral-discretization.git
cd spatial-geometric-pseudospectral-discretization
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[figures]"
python -m unittest discover -s tests -v
python experiments/reproduce.py --check
~~~

To run without generating figures:

~~~bash
python -m pip install -e .
python experiments/reproduce.py --check --no-figures --output /tmp/gpps-results
~~~

The reference run passes **15 unit tests and 22 simulation criteria**.
The eight frequencies published in 2012 are reproduced with a maximum
discrepancy of $`1.3\,10^{-12}`$; measured temporal orders are 2 and 4.
The [results chapter](docs/manuscript.md#5-numerical-protocol-and-results)
reports the parameters, figures, errors and limitations of each experiment.
The thresholds allow for platform-dependent roundoff.

## Repository structure

| Location | Contents |
|---|---|
| [docs/manuscript.md](docs/manuscript.md) | Complete scientific manuscript |
| src/gpps/ | Polynomials, spatial reduction, Hamiltonians, ports and time integrator |
| experiments/ | Simulation campaign and figure generation |
| tests/ | Geometry, energy and convergence checks |
| results/tables/ | Numerical data in CSV format |
| results/figures/ | Six vector figures in SVG format |
| results/summary.json | Measurements and validation criteria |
| results/manifest.json | Numerical environment and source hashes |
| .github/workflows/ | Automated tests and simulation validation |

## Scope and citation

Research version 0.1 covers smooth one-dimensional systems, constant structure
matrices, quadratic, cubic or quartic energies, and synchronous interconnections.
Energy preservation does not guarantee symplecticity, positive Saint-Venant
water depth or shock handling.

The spatial method is attributed to the 2012 authors. The AVF, HBVM and
Petrov–Galerkin time discretizations have their own precedents, discussed
in the [references](docs/manuscript.md#6-references-and-relation-to-prior-work).
See [CITATION.cff](CITATION.cff) and [CONTRIBUTING.md](CONTRIBUTING.md).
