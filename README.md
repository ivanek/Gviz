<img src="vignettes/Gviz-logo.png" align="right" alt="" width="120" />

# _Gviz_ - Plotting data and annotation information along genomic coordinates

## Software status

| In BioC | Travis CI | Bioc ([devel](http://bioconductor.org/packages/devel/bioc/html/Gviz.html)) | Bioc ([release](http://bioconductor.org/packages/release/bioc/html/Gviz.html)) | Coverage |
|:--------|:---------:|:--------------:|:----------------:|:--------:|
| [![In Bioc](https://bioconductor.org/shields/years-in-bioc/Gviz.svg)](http://bioconductor.org/packages/Gviz) | [![Travis CI status](https://travis-ci.org/ivanek/Gviz.svg?branch=master)](https://travis-ci.org/ivanek/Gviz) | [![Bioconductor-devel Build Status](https://bioconductor.org/shields/build/devel/bioc/Gviz.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/Gviz) | [![Bioconductor-release Build Status](https://bioconductor.org/shields/build/release/bioc/Gviz.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/Gviz) | [![Codecov.io coverage status](https://codecov.io/github/ivanek/Gviz/coverage.svg?branch=master)](https://codecov.io/github/ivanek/Gviz) |

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
