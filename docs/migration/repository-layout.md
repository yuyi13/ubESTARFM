# Repository layout migration

The maintained repository uses descriptive paths instead of numbered workflow
directories.

| Published path | Modern path |
|---|---|
| `0_algorithm/ubESTARFM.R` | `legacy/ubESTARFM.R` |
| `0_algorithm/example.R` | `legacy/examples/example.R` |
| `0_algorithm/visualise.R` | `legacy/examples/visualise.R` |
| `1_test_data/` | `inst/extdata/` |
| `2_tmp_path/` | Removed; temporary directories are managed internally |
| `3_output/fused_result.tif` | `tests/fixtures/legacy/fused_result.tif` |
| `3_output/visualisation.png` | `examples/outputs/legacy_visualisation.png` |
| `4_lst_processing_scripts/` | `legacy/paper-processing/` |
| `5_ozflux_lst/` | `data/ozflux/` |
| `README.pdf` | `legacy/README.pdf` |

The published algorithm is preserved byte-for-byte. The former generated fused
result is retained as an immutable numerical acceptance fixture rather than an
example output. New examples write to `examples/outputs/`, whose generated
contents are ignored by Git.
