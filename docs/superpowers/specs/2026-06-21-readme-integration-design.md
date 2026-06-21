# README Integration Design

**Date:** 2026-06-21  
**Status:** Approved

## Purpose

Restore the scientific narrative, demonstration material, figures, provenance,
conference information, and references from the published README while retaining
the modern R/Python package, train-once workflow, batch prediction, caching, and
documentation guidance.

## Structure

The README will use this order:

1. Title and modern badges.
2. Contents.
3. Overview.
4. Background and local-bias-correction figure.
5. Repository structure.
6. Installation for R and Python.
7. Demo and usage:
   - Train once.
   - Predict one target.
   - Predict a batch.
   - Run the maintained R and Python examples.
   - Display the archived demonstration visualization.
8. Model caching and performance.
9. Published R compatibility.
10. LST processing scripts and experimental-design figure.
11. Updated OzFlux processing note.
12. Data and permalinks.
13. Documentation links.
14. Citation in prose and BibTeX.
15. Conference talk.
16. References.

## Content Rules

- Preserve the substance of the original Overview, Background, LST processing,
  OzFlux note, permalink, citation, conference talk, and references sections.
- Update claims that the algorithm is R-only; the maintained implementation is
  available in R and Python.
- Use current package commands and dependencies rather than the deprecated
  `raster`, `foreach`, and `doParallel` setup.
- Explain that reference-pair candidate relationships are trained once and
  reused across target dates.
- Keep the README self-contained while linking detailed tutorials, benchmarks,
  migration notes, and archived materials.
- State clearly that archived paper-processing scripts require external study
  data and are retained for reference rather than direct execution.

## Path Translation

| Published README path | Current repository path |
|---|---|
| `0_algorithm/ubESTARFM.R` | `legacy/ubESTARFM.R` |
| `0_algorithm/example.R` | `legacy/examples/example.R` |
| `0_algorithm/visualise.R` | `legacy/examples/visualise.R` |
| `1_test_data/` | `inst/extdata/` |
| `3_output/visualisation.png` | `examples/outputs/legacy_visualisation.png` |
| `4_lst_processing_scripts/` | `legacy/paper-processing/` |
| `5_ozflux_lst/` | `data/ozflux/` |

## Validation

- Every Contents link must resolve to a README heading.
- Every referenced repository path must exist.
- R and Python quick-start snippets must use current public APIs.
- Badges, DOI links, article links, dataset links, figures, citation details,
  conference information, and references must remain valid.
- The change is limited to README content and this design/plan documentation.
