# Migrating from the published R API

The published `ubESTARFM()` function remains available. It now delegates to the
maintained train/predict API:

```r
model <- ubestarfm_train(
  fine_1,
  fine_2,
  coarse_1,
  coarse_2,
  window_radius = 25L,
  patch_size = 200L,
  method = "zero_bias"
)

prediction <- ubestarfm_predict(
  model,
  coarse_target,
  value_range = c(250, 350)
)
```

Use `ubestarfm_predict_batch()` to reuse the same model for multiple target
dates. The published `tmp_path` argument is accepted by `ubESTARFM()` for source
compatibility but is intentionally unused; modern code no longer discovers and
averages every GeoTIFF in a shared temporary directory.

Method names in the modern API are `"zero_bias"` and `"baseline"`. The
compatibility wrapper continues to accept `"zero bias"` and `"baseline"`.

The maintained implementation handles target-date missing values during
prediction. Candidate weights are renormalized after target filtering, and the
published low-candidate fallback is used when five or fewer candidates remain.
