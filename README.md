# MetaScope [![](https://img.shields.io/badge/bioconductor-3.18-blue)](https://www.bioconductor.org/packages/release/data/experiment/html/MetaScope.html)

## What is MetaScope?

MetaScope is a complete R-based 16S, metagenomic, and metatranscriptomic profiling package that can accurately identify the composition of microbes at a strain-level resolution within a sample. MetaScope can be considered as an updated and expanded R translation of [PathoScope 2.0](https://microbiomejournal.biomedcentral.com/articles/10.1186/2049-2618-2-33), a Python-based metagenomic profiling package also created by the Johnson lab. 

A few improvements made in MetaScope over PathoScope include using the BAM file format instead of the SAM file format for significantly less disk space usage, removing all dependencies to NCBI’s now defunct GI sequence annotations, and properly filtering reads that align to filter reference genomes. Functions to analyze host microbiome data are also planned to be added in future updates to the package.

## Documentation
Documentation and tutorials for MetaScope are available at our [website](https://wejlab.github.io/metascope-docs/).

## Installation

MetaScope requires R Version 4.2.

Install the development version of the package from Github:

```
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
devtools::install_github("wejlab/MetaScope")
```

Install the release version of the package from Bioconductor:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MetaScope")
```
