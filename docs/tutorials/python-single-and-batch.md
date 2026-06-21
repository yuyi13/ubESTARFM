# Python: single and batch prediction

Install the package from the repository root:

```bash
uv pip install ./python
```

Train once and predict one target:

```python
from ubestarfm import predict, train

model = train(
    "inst/extdata/Landsat_LST_cloudrm_20160205.tif",
    "inst/extdata/Landsat_LST_cloudrm_20160308.tif",
    "inst/extdata/MOD11A1_LST_cloudrm_20160205.tif",
    "inst/extdata/MOD11A1_LST_cloudrm_20160308.tif",
    window_radius=25,
    patch_size=200,
    workers=4,
)

prediction = predict(
    model,
    "inst/extdata/MOD11A1_LST_cloudrm_20160218.tif",
    output_path="examples/outputs/python_single_20160218.tif",
)
```

Use `predict_batch(model, targets, output_paths=...)` to reuse the same
candidate model for many target dates. `save_model()` writes a compressed NPZ
archive and `load_model()` validates and restores it.
