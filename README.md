# Unbiased ESTARFM (ubESTARFM) in R

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Programming](https://img.shields.io/badge/-R%20Programming-3776AB?style=flat&logo=R&logoColor=white)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8017282.svg)](https://doi.org/10.5281/zenodo.8017282)
[![Twitter Follow](https://img.shields.io/twitter/follow/yy_nash13?style=social)](https://twitter.com/yy_nash13)

## Contents

- [Overview](#overview)
- [Background](#background)
- [Usage](#usage)
- [Permalink](#permalink)
- [LST processing scripts](#lst-processing-scripts)
- [IMPORTANT update for *in-situ* LST from OzFlux](#important-update-for-in-situ-lst-from-ozflux)
- [To cite ubESTARFM](#to-cite-ubestarfm)
- [References](#references)

## Overview

This is the repository for the ubESTARFM algorithm. The algorithm is described in detail in our paper published in *Remote Sensing of Environment* ([Yu et al., 2023](https://doi.org/10.1016/j.rse.2023.113784)). Apart from the [algorithm](#usage), we also provided the [scripts for processing and fusing satellite LST data](#lst-processing-scripts), as well as **AN IMPORTANT UPDATE** for processing [*in-situ* LST data from OzFlux sites](#important-update-for-in-situ-lst-from-ozflux).

## Background

Fine spatial resolution land surface temperature (LST) data are crucial to study heterogeneous landscapes (e.g., agricultural and urban). Some well-known spatiotemporal fusion methods like the Spatial and Temporal Adaptive Reflectance Fusion Model (STARFM; Gao et al., 2006) and the Enhanced STARFM (ESTARFM; Zhu et al., 2010), which were originally developed to fuse surface reflectance data, may not be suitable for direct application in LST studies due to the high sub-diurnal dynamics of LST. To address this, we proposed a variant of ESTARFM, referred to as the unbiased ESTARFM (ubESTARFM), specifically designed to accommodate the high temporal dynamics of LST to generate fine-resolution LST estimates. 

In ubESTARFM, we implement a local bias correction on the central pixel and similar fine-resolution pixels within the moving window using the mean value of corresponding coarse-resolution pixels as reference. By applying this linear scaling approach, we can scale the systematic biases of the fine-resolution data to a same level of the corresponding coarse resolution data in each moving window, while maintaining the variation and spatial details of fine resolution data.

![](figures/local-bias-correction.png)

## Usage

The ubESTARFM algorithm is written in R. We recommend users to use a multi-core processor that can allow ubESTARFM to run in parallel and to be more efficient.

Please install essential R packages before running ubESTARFM. 

```
install.packages('raster')
install.packages('foreach')
install.packages('doParallel')
```

To see an example of ubESTARFM, **please make sure you are under the directory** `ubESTARFM/`, then simply run the following via the command line:

```
Rscript 0_code/example.R
```

This will run ubESTARFM on a small subset of data (Yanco site) using 4 cores and generate a `fused_result.tif` in the directory `3_output/`.

Have a look at the result:

```
Rscript 0_code/visualise.R
```

This will generate a `visualisation.png` in output that looks like:

![](3_output/visualisation.png)

Please note the data included in this repository are for demonstration purposes only.

## Permalink

If you are interested in having a comprehensive assessment of ubESTARFM, please refer to the dataset published in the [CSIRO Data Access Portal](https://doi.org/10.25919/b77m-8n31), which contains the full set of data (12 OzFlux sites across Australia for the period of 2013-2021) used in our RSE paper.

The published link of this code is at [Zenodo](https://doi.org/10.5281/zenodo.8017282). You can also find an archived version at [ResearchGate](https://www.researchgate.net/publication/371376456_Unbiased_ESTARFM_ubESTARFM).

## LST processing scripts

The scripts for processing and fusing satellite LST are archived in `4_lst_processing` for **reference purposes only**. The scripts are ordered in sequence 0-9, which follows the experimental design as below. However, it is unlikely you can run the scripts directly as the input data are massive and not available here.

![](figures/experimental-design.png)

## IMPORTANT update for *in-situ* LST from OzFlux

TBC

## To cite ubESTARFM

If you found this repository helpful, please kindly consider citing:

- Yu, Y., Renzullo, L. J., McVicar, T. R., Malone, B. P. and Tian, S., 2023. Generating daily 100 m resolution land surface temperature estimates continentally using an unbiased spatiotemporal fusion approach. *Remote Sensing of Environment, 297*, 113784. https://doi.org/10.1016/j.rse.2023.113784

## References

- Gao, F., Masek, J., Schwaller, M. and Hall, F., 2006. On the blending of the Landsat and MODIS surface reflectance: Predicting daily Landsat surface reflectance. *IEEE Transactions on Geoscience and Remote Sensing, 44*, 2207-2218. https://doi.org/10.1109/TGRS.2006.872081

- Zhu, X., Chen, J., Gao, F., Chen, X. and Masek, J. G., 2010. An enhanced spatial and temporal adaptive reflectance fusion model for complex heterogeneous regions. *Remote Sensing of Environment, 114*, 2610-2623. https://doi.org/10.1016/j.rse.2010.05.032
