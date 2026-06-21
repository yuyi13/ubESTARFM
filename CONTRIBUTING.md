# Contributing

Create changes on a feature branch and keep the published
`legacy/ubESTARFM.R` file byte-for-byte unchanged.

Run R verification:

```bash
Rscript -e 'Rcpp::compileAttributes(".")'
Rscript -e 'roxygen2::roxygenise(".")'
R CMD build .
R CMD check --no-manual ubESTARFM_3.0.0.9000.tar.gz
```

Run Python and cross-language verification:

```bash
uv pip install --python "$(command -v python3)" -e 'python[test]'
python3 -m ruff check python tests
python3 -m pytest python/tests tests -v
```

All new or revised R, Python, and shell scripts must follow the workspace script
header standard. Maintained plotting code uses Arial where the plotting system
supports explicit font selection.
