<img src="vignettes/Gviz-logo.png" align="right" alt="" width="120" />

# _Gviz_ - Plotting data and annotation information along genomic coordinates

## Software status


[![How long since the package was first in a released Bioconductor version](https://bioconductor.org/shields/years-in-bioc/Gviz.svg)](https://bioconductor.org/packages/Gviz) 
[![Support site activity in last 6 months: agged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts](https://bioconductor.org/shields/posts/Gviz.svg)](https://support.bioconductor.org/t/gviz/)


| Platform  |      OS   |    Last Update |      Status      |      Coverage    | Downloads |
|:----------|:---------:|:--------------:|:----------------:|:----------------:|:---------:|
| Travis CI |   Linux   | ![GitHub last commit](https://img.shields.io/github/last-commit/ivanek/Gviz) | [![Travis CI](https://travis-ci.org/ivanek/Gviz.svg?branch=master)](https://travis-ci.org/ivanek/Gviz) | [![Codecov.io coverage status](https://codecov.io/github/ivanek/Gviz/coverage.svg?branch=master)](https://codecov.io/github/ivanek/Gviz) | NA |
| Bioc ([devel](http://bioconductor.org/packages/devel/bioc/html/Gviz.html))     | [![Platforms](https://bioconductor.org/shields/availability/devel/Gviz.svg)](https://bioconductor.org/packages/Gviz#archives) | [![Bioconductor-devel Last Commit](https://bioconductor.org/shields/lastcommit/devel/bioc/Gviz.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/Gviz/) | [![Bioconductor-devel Build Status](https://bioconductor.org/shields/build/devel/bioc/Gviz.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/Gviz) | NA | [![Bioconductor-devel Downloads](https://bioconductor.org/shields/downloads/devel/Gviz.svg)](http://bioconductor.org/packages/stats/bioc/Gviz/) |
| Bioc ([release](http://bioconductor.org/packages/release/bioc/html/Gviz.html)) | [![Platforms](https://bioconductor.org/shields/availability/release/Gviz.svg)](https://bioconductor.org/packages/Gviz#archives) | [![Bioconductor-release Last Commit](https://bioconductor.org/shields/lastcommit/release/bioc/Gviz.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/Gviz/) | [![Bioconductor-release Build Status](https://bioconductor.org/shields/build/release/bioc/Gviz.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/Gviz) | NA | [![Bioconductor-release Downloads](https://bioconductor.org/shields/downloads/release/Gviz.svg)](http://bioconductor.org/packages/stats/bioc/Gviz/) |



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

Genomic data analyses requires integrated visualization of known genomic information and new experimental data. Gviz uses the [biomaRt](http://bioconductor.org/packages/biomaRt/) and the [rtracklayer](http://bioconductor.org/packages/rtracklayer/) packages to perform live annotation queries to [Ensembl](https://www.ensembl.org/) and [UCSC](https://genome.ucsc.edu) and translates this to e.g. gene/transcript structures in viewports of the grid graphics package. This results in genomic information plotted together with your data.

## Installation

### Release version

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Gviz", version = "release")
```

### Developmental version

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Gviz", version = "devel")
```

### Github

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("ivanek/Gviz")
```


## Citation 

- Hahne F, Ivanek R (2016). "Statistical Genomics: Methods and Protocols." In Mathé E, Davis S (eds.), chapter Visualizing Genomic Data Using Gviz and Bioconductor, 335–351. Springer New York, New York, NY. ISBN 978-1-4939-3578-9, doi: [10.1007/978-1-4939-3578-9_16](http://dx.doi.org/10.1007/978-1-4939-3578-9_16).
