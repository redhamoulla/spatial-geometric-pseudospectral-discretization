# Contributing

Changes to the geometry must retain explicit sign conventions, dimensions
and power pairings. Do not silently repair a matrix by skew-symmetrizing it
after assembly.

For a new Hamiltonian, provide its value, gradient, Hessian and constitutive
quadrature. Distinguish structural, temporal and spatial errors.
Energy preservation alone is not a trajectory accuracy test.

~~~bash
python -m unittest discover -s tests -v
python experiments/reproduce.py --check --no-figures --output /tmp/gpps-review
~~~

Document positive and negative results. Cite existing methods; do not describe
an observed property as proved.

## Documentation and equations

Write documentation in English and keep the scientific discussion in
[the single manuscript](docs/manuscript.md).

Use fenced math blocks for display equations, so that Markdown preserves
matrix row separators and LaTeX escapes. Put equation identifiers such as
**(T1)** in the surrounding prose. Avoid LaTeX equation tags: they can render
as vertical stacks in GitHub's browser output.

Use GitHub's dollar-and-backtick delimiters for inline mathematics containing
underscores or backslashes. Check the rendered GitHub preview, including
matrices and numbered equations, after changing mathematical notation.
See [GitHub's math formatting documentation](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).
