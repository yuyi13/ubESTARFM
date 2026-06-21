# Model caching and memory

The failed experimental fast implementation stored a nested object for every
pixel and duplicated candidate coordinates, values, and weights. On the bundled
400 by 400 example, an exact weight cache would require multiple gigabytes.

The maintained engines store a fixed packed candidate bitset per pixel. With a
51 by 51 window, the full reusable model is approximately 56 MiB in memory.
Weights are reconstructed during prediction so target-date missing values can
be filtered correctly.

Use a saved model when target rasters arrive in later sessions:

```r
ubestarfm_save_model(model, "reference_pair.rds")
model <- ubestarfm_load_model("reference_pair.rds")
```

```python
save_model(model, "reference_pair.npz")
model = load_model("reference_pair.npz")
```

For an immediately available batch, `ubestarfm_fit_predict_batch()` and
`fit_predict_batch()` provide a concise train-once workflow.
