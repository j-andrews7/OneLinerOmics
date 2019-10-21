# AndrewsBCellLymphoma

[![Docs](https://img.shields.io/badge/docs-site-brightgreen?style=plastic)](https://j-andrews7.github.io/AndrewsBCellLymphoma)

This R package contains several functions that compress major RNA-seq and ChIP-seq analysis steps into one-liners (or two, I'm not perfect).
It is meant to serve as a complement to the (future) publication and as a record of how the data was analyzed in addition to serving as a basis for the Payton Lab bulk RNA-seq and ChIP-seq analyses moving forwards.

## Installation

```r
require("devtools")
devtools::install_github("j-andrews7/AndrewsBCellLymphoma")
```

## Usage

Please see the [docs](https://j-andrews7.github.io/AndrewsBCellLymphoma) for a full reference and vignettes with examples (including steps like gene count quantification). 
The [salmon](https://combine-lab.github.io/salmon/getting_started/), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), and [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) manuals and vignettes may also be helpful for understanding what's going on under the hood and interpreting the output.

## Issues

Please direct all complaints to [management](https://github.com/j-andrews7/AndrewsBCellLymphoma/issues). 