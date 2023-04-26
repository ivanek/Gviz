<img src="vignettes/Gviz-logo.png" align="right" alt="" width="120" />

# _Gviz_ - Plotting data and annotation information along genomic coordinates

## Software status

<!-- badges: start -->
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![How long since the package was first in a released Bioconductor version](https://bioconductor.org/shields/years-in-bioc/Gviz.svg)](https://bioconductor.org/packages/Gviz) 
[![Bioconductor-devel Downloads](https://bioconductor.org/shields/downloads/devel/Gviz.svg)](https://bioconductor.org/packages/stats/bioc/Gviz/)
[![Support site activity in last 6 months: agged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts](https://bioconductor.org/shields/posts/Gviz.svg)](https://support.bioconductor.org/t/gviz/)
[![License:
Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0)
[![GitHub issues](https://img.shields.io/github/issues/ivanek/Gviz)](https://github.com/ivanek/Gviz/issues)
<!-- badges: end -->

&nbsp;

|                     | Bioc ([release](https://bioconductor.org/packages/release/bioc/html/Gviz.html)) | Bioc ([devel](https://bioconductor.org/packages/devel/bioc/html/Gviz.html)) |
|:--------------------|----------------------------------------------------------------------------:|--------------------------------------------------------------------------------:|
| OS                  | [![Platforms](https://bioconductor.org/shields/availability/release/Gviz.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/Gviz/) | [![Platforms](https://bioconductor.org/shields/availability/devel/Gviz.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/Gviz/) |
| Bioc Last Update    | [![Bioconductor-release Last Commit](https://bioconductor.org/shields/lastcommit/release/bioc/Gviz.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/Gviz/) | [![Bioconductor-devel Last Commit](https://bioconductor.org/shields/lastcommit/devel/bioc/Gviz.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/Gviz/) |
| Bioc Status         | [![Bioconductor-release Build Status](https://bioconductor.org/shields/build/release/bioc/Gviz.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/Gviz) | [![Bioconductor-devel Build Status](https://bioconductor.org/shields/build/devel/bioc/Gviz.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/Gviz) |
| GitHub Last Commit  | [![GitHub last commit (Bioconductor-release)](https://img.shields.io/github/last-commit/ivanek/Gviz/RELEASE_3_17)](https://github.com/ivanek/Gviz/tree/RELEASE_3_17) | [![GitHub last commit (Bioconductor-devel)](https://img.shields.io/github/last-commit/ivanek/Gviz/devel)](https://github.com/ivanek/Gviz/tree/devel/) |
| GitHub Actions      | [![R build status](https://github.com/ivanek/Gviz/workflows/R-CMD-check-bioc/badge.svg?branch=RELEASE_3_17)](https://github.com/ivanek/Gviz/actions) | [![R build status](https://github.com/ivanek/Gviz/workflows/R-CMD-check-bioc/badge.svg?branch=devel)](https://github.com/ivanek/Gviz/actions) |
| Coverage            | [![Codecov.io (Bioconductor-release)](https://codecov.io/github/ivanek/Gviz/coverage.svg?branch=RELEASE_3_17)](https://codecov.io/gh/ivanek/Gviz/branch/RELEASE_3_17) | [![Codecov.io (Bioconductor-devel)](https://codecov.io/github/ivanek/Gviz/coverage.svg?branch=devel)](https://codecov.io/github/ivanek/Gviz) |

## Authors

- Florian Hahne
- Steffen Durinck
- Robert Ivanek
- Arne Mueller
- Steve Lianoglou
- Ge Tan 
- Lance Parsons
- Shraddha Pai

## Overview

![Gviz UCSC like screenshot](vignettes/Gviz-example.png)

Genomic data analyses requires integrated visualization of known genomic information and new experimental data. Gviz uses the [biomaRt](https://bioconductor.org/packages/biomaRt/) and the [rtracklayer](https://bioconductor.org/packages/rtracklayer/) packages to perform live annotation queries to [Ensembl](https://www.ensembl.org/) and [UCSC](https://genome.ucsc.edu) and translates this to e.g. gene/transcript structures in viewports of the grid graphics package. This results in genomic information plotted together with your data.

## Installation

#### Release version

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Gviz", version = "release")
```

#### Developmental version

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Gviz", version = "devel")
```

#### Github

```
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("ivanek/Gviz")
```
## Usage

For detailed instructions check the package vignette 
([release](https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html) 
or 
[developmental](https://bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.html) 
version). Or check the [GitHub.io page](https://ivanek.github.io/Gviz/).

## Citation 

- Hahne F, Ivanek R (2016). "Statistical Genomics: Methods and Protocols." In Mathé E, Davis S (eds.), chapter Visualizing Genomic Data Using Gviz and Bioconductor, 335–351. Springer New York, New York, NY. ISBN 978-1-4939-3578-9, doi: [10.1007/978-1-4939-3578-9_16](https://dx.doi.org/10.1007/978-1-4939-3578-9_16).
