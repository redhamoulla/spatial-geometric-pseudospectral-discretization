# Reference results

These files were produced by running experiments/reproduce.py.
Numerical results are computed, not prefilled. Only the frequencies published
in 2012 are stored as an independent reference.

- summary.json: measurements and 22 validation criteria.
- manifest.json: software versions and source hashes.
- tables/: balances, convergence, spectra, profiles and condition numbers.
- figures/: SVG plots generated from the CSV files.

Units are normalized except for the explicit constants used to reproduce
the 2012 results. Interval length and inductance are distinct parameters.
Reproduction does not require bit-for-bit agreement across platforms:
it checks the numerical criteria and records source hashes.

[Results and discussion](../docs/manuscript.md#5-numerical-protocol-and-results).
