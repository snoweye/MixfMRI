# MixfMRI

* **License:** [![License](http://img.shields.io/badge/license-MPL%202-orange.svg?style=flat)](https://www.mozilla.org/MPL/2.0/)
* **Download:** [![Download](http://cranlogs.r-pkg.org/badges/MixfMRI)](https://cran.r-project.org/package=MixfMRI)
* **Status:** [![Build Status](https://travis-ci.org/snoweye/MixfMRI.png)](https://travis-ci.org/snoweye/MixfMRI) [![Appveyor Build status](https://ci.appveyor.com/api/projects/status/32r7s2skrgm9ubva?svg=true)](https://ci.appveyor.com/project/snoweye/MixfMRI)
* **Authors:** Wei-Chen Chen and Ranjan Maitra

MixfMRI is an R package clustering fMRI data using parallel model-based
clustering. The package mplements methods in Chen and Maitra (2021),
arXiv:2102.03639 ([pdf](https://arxiv.org/pdf/2102.03639)).
The developed methods include 2D and 3D clustering analyses
(for p-values with voxel locations) and
segmentation analyses (for p-values alone) for fMRI data where p-values
indicate significant level of activation responding to stimulate
of interesting. The analyses are mainly identifying active
voxel/signal associated with normal brain behaviors.


## Installation

For installing MixfMRI:
* See "./INSTALL" for details.

## Usage

For simulation studies (with MPI and pbdMPI):
* See "./inst/workflow/simulation/create_simu.txt" for simulations
  that support the methods developed and evaluated
  in Chen and Maitra (2021), arXiv:2102.03639.

For testing run time (with MPI and pbdMPI):
* See "./inst/workflow/timing/" for testing a main function.

## Workflow

For large scale workflows and examples (with MPI and pbdMPI):
* "./inst/workflow/spmd/*_common/" for workflow main scripts.
* "./inst/workflow/spmd/2d_V/" for running a 2D_V exmaple.
* "./inst/workflow/spmd/2d_I/" for running a 2D_I exmaple.
* "./inst/workflow/spmd/2d_noX/" for running a 2D_noX exmaple.
* "./inst/workflow/spmd/3d_V/" for running a 3D_V exmaple.
* "./inst/workflow/spmd/3d_I/" for running a 3D_I exmaple.
* "./inst/workflow/spmd/3d_noX/" for running a 3D_noX exmaple.

