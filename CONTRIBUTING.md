# Contributing

Open an issue before changing statistical assumptions or output schemas. Keep
biological labels configurable, never use target hybrids to select diagnostic
markers, and add a synthetic regression test for every estimator change.

Pull requests should pass `python -m py_compile workflow/scripts/*.py`, `pytest`,
and `snakemake --lint`.
