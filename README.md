# Unbiased ESTARFM

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![R 4.3+](https://img.shields.io/badge/R-4.3%2B-276DC3?style=flat&logo=R&logoColor=white)](https://www.r-project.org/)
[![Python 3.12+](https://img.shields.io/badge/Python-3.12%2B-3776AB?style=flat&logo=Python&logoColor=white)](https://www.python.org/)
[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.8017282-1682D4?style=flat)](https://doi.org/10.5281/zenodo.8017282)
[![Article](https://img.shields.io/badge/Article-Remote%20Sensing%20of%20Environment-0A70A3?style=flat)](https://doi.org/10.1016/j.rse.2023.113784)
[![Dataset](https://img.shields.io/badge/Dataset-CSIRO%20Data%20Access%20Portal-00A9CE?style=flat)](https://doi.org/10.25919/rrpg-m948)

ubESTARFM generates fine-resolution land surface temperature estimates from two
fine/coarse reference pairs and a coarse target raster. It applies the local
bias correction described by Yu et al. (2023).

The maintained implementation trains reference-pair candidate relationships
once and reuses them for any number of target dates. The published R
implementation is preserved unchanged in
[`legacy/ubESTARFM.R`](legacy/ubESTARFM.R).

## R quick start

```bash
Rscript -e 'install.packages(c("terra", "Rcpp"))'
R CMD INSTALL .
```

```r
library(ubESTARFM)

model <- ubestarfm_train(
  "inst/extdata/Landsat_LST_cloudrm_20160205.tif",
  "inst/extdata/Landsat_LST_cloudrm_20160308.tif",
  "inst/extdata/MOD11A1_LST_cloudrm_20160205.tif",
  "inst/extdata/MOD11A1_LST_cloudrm_20160308.tif",
  workers = 4L
)

prediction <- ubestarfm_predict(
  model,
  "inst/extdata/MOD11A1_LST_cloudrm_20160218.tif"
)
```

## Python quick start

```bash
uv pip install ./python
```

```python
from ubestarfm import predict, train

model = train(
    "inst/extdata/Landsat_LST_cloudrm_20160205.tif",
    "inst/extdata/Landsat_LST_cloudrm_20160308.tif",
    "inst/extdata/MOD11A1_LST_cloudrm_20160205.tif",
    "inst/extdata/MOD11A1_LST_cloudrm_20160308.tif",
    workers=4,
)
prediction = predict(
    model,
    "inst/extdata/MOD11A1_LST_cloudrm_20160218.tif",
)
```

## Batch prediction

Use `ubestarfm_predict_batch()` in R or `predict_batch()` in Python. Candidate
search is performed once during training, while temporal differences and
target-date missing values remain prediction-specific.

The packed candidate cache is approximately 56 MiB for the bundled 400 by 400
example with a 51 by 51 window. Exact per-pixel weight caching would require
multiple gigabytes and is intentionally avoided.

## Documentation

- [R tutorial](docs/tutorials/r-single-and-batch.md)
- [Python tutorial](docs/tutorials/python-single-and-batch.md)
- [Model caching](docs/tutorials/model-caching.md)
- [Published R API migration](docs/migration/from-published-r-api.md)
- [Repository layout migration](docs/migration/repository-layout.md)
- [Benchmarks](benchmarks/)
- [Archived paper-processing scripts](legacy/paper-processing/)

## Citation

Yu, Y., Renzullo, L. J., McVicar, T. R., Malone, B. P., and Tian, S. (2023).
Generating daily 100 m resolution land surface temperature estimates
continentally using an unbiased spatiotemporal fusion approach.
*Remote Sensing of Environment*, 297, 113784.
https://doi.org/10.1016/j.rse.2023.113784

## Data

The bundled rasters are demonstration subsets. The full study dataset is
available from the
[CSIRO Data Access Portal](https://doi.org/10.25919/rrpg-m948).
