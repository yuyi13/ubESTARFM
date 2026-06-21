# R: single and batch prediction

Install directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("yuyi13/ubESTARFM")
library(ubestarfm)
```

For development from a local clone, run `R CMD INSTALL .`.

Train the two reference dates once:

```r
library(ubestarfm)

model <- ubestarfm_train(
  fine_1 = "inst/extdata/Landsat_LST_cloudrm_20160205.tif",
  fine_2 = "inst/extdata/Landsat_LST_cloudrm_20160308.tif",
  coarse_1 = "inst/extdata/MOD11A1_LST_cloudrm_20160205.tif",
  coarse_2 = "inst/extdata/MOD11A1_LST_cloudrm_20160308.tif",
  window_radius = 25L,
  patch_size = 200L,
  workers = 4L
)
```

Predict one date:

```r
prediction <- ubestarfm_predict(
  model,
  "inst/extdata/MOD11A1_LST_cloudrm_20160218.tif",
  output_path = "examples/outputs/r_single_20160218.tif"
)
```

For multiple dates, pass every target to one batch call. Candidate search is
not repeated:

```r
predictions <- ubestarfm_predict_batch(
  model,
  list(target_1, target_2, target_3),
  output_paths = c(output_1, output_2, output_3),
  workers = 4L
)
```

Save the model with `ubestarfm_save_model()` and restore it later with
`ubestarfm_load_model()`.
