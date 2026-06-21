# Unbiased ESTARFM (ubESTARFM) in R and Python

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![R 4.3+](https://img.shields.io/badge/R-4.3%2B-276DC3?style=flat&logo=R&logoColor=white)](https://www.r-project.org/)
[![Python 3.12+](https://img.shields.io/badge/Python-3.12%2B-3776AB?style=flat&logo=Python&logoColor=white)](https://www.python.org/)
[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.8017282-1682D4?style=flat)](https://doi.org/10.5281/zenodo.8017282)
[![Article](https://img.shields.io/badge/Article-Remote%20Sensing%20of%20Environment-0A70A3?style=flat)](https://doi.org/10.1016/j.rse.2023.113784)
[![Dataset](https://img.shields.io/badge/Dataset-CSIRO%20Data%20Access%20Portal-00A9CE?style=flat)](https://doi.org/10.25919/rrpg-m948)

## Contents

- [Overview](#overview)
- [Background](#background)
- [Repository structure](#repository-structure)
- [Installation](#installation)
- [Demo and usage](#demo-and-usage)
  - [Train once](#train-once)
  - [Predict one target](#predict-one-target)
  - [Predict a batch](#predict-a-batch)
  - [Run the maintained examples](#run-the-maintained-examples)
  - [Demonstration result](#demonstration-result)
- [Model caching and performance](#model-caching-and-performance)
- [Published R compatibility](#published-r-compatibility)
- [LST processing scripts](#lst-processing-scripts)
  - [Important note for in-situ LST from OzFlux](#important-note-for-in-situ-lst-from-ozflux)
- [Data and permalinks](#data-and-permalinks)
- [Documentation](#documentation)
- [Citation](#citation)
- [Conference talk](#conference-talk)
- [References](#references)

## Overview

This repository contains the unbiased Enhanced Spatial and Temporal Adaptive
Reflectance Fusion Model (ubESTARFM), described in
[Yu et al. (2023)](https://doi.org/10.1016/j.rse.2023.113784). It also retains
the scripts used to process and evaluate the land surface temperature (LST)
data in that study.

The maintained implementation is available in R and Python. It separates the
algorithm into two stages: train the relationships from two fine/coarse
reference pairs once, then reuse that model to predict one or many intervening
target dates. The published R implementation remains preserved byte-for-byte
in [`legacy/ubESTARFM.R`](legacy/ubESTARFM.R).

## Background

Fine-spatial-resolution LST data are crucial for studying heterogeneous
landscapes such as agricultural and urban areas. Well-known spatiotemporal
fusion methods including the Spatial and Temporal Adaptive Reflectance Fusion
Model (STARFM; Gao et al., 2006) and Enhanced STARFM (ESTARFM; Zhu et al.,
2010) were originally developed to fuse surface reflectance. Direct
application to LST can be unsuitable because LST has strong sub-daily
dynamics.

ubESTARFM is an ESTARFM variant designed for those temporal dynamics. It
applies local bias correction to the central pixel and similar
fine-resolution pixels within each moving window, using the mean of the
corresponding coarse-resolution pixels as the reference. This linear scaling
places systematic fine-resolution bias on the same level as the corresponding
coarse-resolution data while retaining fine-resolution variation and spatial
detail.

![Schematic of the ubESTARFM local bias correction](figures/local-bias-correction.png)

## Repository structure

The repository now uses descriptive paths rather than the numbered directories
in the published release:

| Purpose | Current path |
|---|---|
| Maintained R package | [`R/`](R/) and [`src/`](src/) |
| Maintained Python package | [`python/`](python/) |
| Bundled demonstration rasters | [`inst/extdata/`](inst/extdata/) |
| Maintained examples | [`examples/R/`](examples/R/) and [`examples/python/`](examples/python/) |
| Published R algorithm | [`legacy/ubESTARFM.R`](legacy/ubESTARFM.R) |
| Published example and visualization scripts | [`legacy/examples/`](legacy/examples/) |
| Archived paper-processing scripts | [`legacy/paper-processing/`](legacy/paper-processing/) |
| Updated OzFlux LST data | [`data/ozflux/`](data/ozflux/) |
| Tutorials and migration notes | [`docs/`](docs/) |
| Cross-language and regression tests | [`tests/`](tests/) |

See the [repository-layout migration guide](docs/migration/repository-layout.md)
for the complete old-to-new path mapping.
