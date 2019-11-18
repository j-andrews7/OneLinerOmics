# OneLinerOmics

[![Docs](https://img.shields.io/badge/docs-site-brightgreen?style=plastic)](https://j-andrews7.github.io/OneLinerOmics)

This R package contains several functions that compress major RNA-seq and ChIP-seq analysis steps into one-liners (okay, maybe two, I'm not perfect).

## Objective
This package serves as a complement to an in-preparation manuscript and as a record of how the data was analyzed in addition to serving as a basis for the Payton Lab bulk RNA-seq and ChIP-seq analyses moving forwards.

## What It Can Do
Ultimately, this package tries to make generic RNA-seq and ChIP-seq analyses as high-level and straightforward as possible. 
It does so by requiring only the path to a samplesheet containing sample metadata and file locations and the metadata variable of interest. 
It will then perform RNA or ChIP-seq analyses with `DESeq2` or `DiffBind`, respectively, for all possible comparisons for the variable of interest and save the results for each analysis to disk.

You can feed it multiple p-value/FDR/fold change thresholds, and it will generate fairly nice plots for differentially expressed genes or differentially bound regions including PCA/MA plots, heatmaps, volcano plots, GO/pathway enrichment analyses, boxplots, etc., saving all of them as PDFs.

## Installation

```r
require("devtools")
devtools::install_github("j-andrews7/OneLinerOmics")
```

## Usage

Please see the [docs](https://j-andrews7.github.io/OneLinerOmics) for a full reference and vignettes with examples (including steps like gene count quantification, peak calling, etc). 
90% of the functionality can be viewed with `?RunDESeq2` and `?RunDiffBind` after installation.
The [salmon](https://combine-lab.github.io/salmon/getting_started/), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), and [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) manuals and vignettes may also be helpful for understanding what's going on under the hood and interpreting the output.

## Issues

Please direct all complaints to [management](https://github.com/j-andrews7/OneLinerOmics/issues). 