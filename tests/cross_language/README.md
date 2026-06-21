# Cross-language acceptance

`test_cross_language.py` trains the bundled reference pair independently in R
and Python, predicts the February 18 target, and compares both outputs with the
published tracked result.

Run:

```bash
python3 -m pytest tests/cross_language/test_cross_language.py -v
```

`compare_outputs.py` is a standalone diagnostic command for comparing any two
GeoTIFF outputs.
