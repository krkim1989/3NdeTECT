.PHONY: test check dry-run smoke-data smoke smoke-installed

test:
	python -m pytest -q

check:
	python -m py_compile workflow/scripts/*.py workflow/lib.py
	python -m pytest -q
	snakemake --configfile tests/data/dryrun.yaml --lint
	snakemake --configfile tests/data/dryrun.yaml --cores 1 --dry-run --quiet

dry-run:
	snakemake --configfile tests/data/dryrun.yaml --cores 1 --dry-run

smoke-data:
	bash tests/smoke/prepare_indexes.sh

smoke: smoke-data
	snakemake --configfile tests/smoke/config.yaml --use-conda --cores 4
	python tests/smoke/verify_smoke.py

# For clean-room/CI tests where core and ploidy executables are already on PATH.
smoke-installed: smoke-data
	snakemake --configfile tests/smoke/config.yaml --cores 4
	python tests/smoke/verify_smoke.py
