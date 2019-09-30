# AndrewsBCellLymphoma

This R package contains several functions that compress major RNA-seq and ChIP-seq analysis steps into one-liners.
It is meant to serve as a complement to the (future) publication and as a record of how the data was analyzed.

## Installation

```r
require(devtools)
devtools::install_github("j-andrews7/AndrewsBCellLymphoma")
```

## Usage

Please see the [vignette](). 
The [salmon](https://combine-lab.github.io/salmon/getting_started/), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), and [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) manuals and vignettes may also be helpful for understanding what's going on under the hood and interpreting the output.