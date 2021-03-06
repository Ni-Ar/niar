# My personal `R` package

Here, I collect R functions that I’ve created over time. This is not really stable or robust, the primary purpose of this package is to neatly organise my R functions in one place and have them ready to use whenever. Feedback is welcome! If you find an error please [open a GitHub issue](https://github.com/Ni-Ar/niar/issues/new).

## Installation

### From this GitHub Repository

If interested in giving it a try:

```R
devtools::install_github("Ni-Ar/niar")
```

### Load the local repository (only for CRG users)

If you have a [CRG](https://www.crg.eu/) user account, instead of installing, you can load the R package at the beginning of your session. Log in with your credentials on the [CRG RStudio Server IDE](http://rstudio4.linux.crg.es/) (R version 4.0.3) and run:
```R
devtools::load_all(path = '/users/mirimia/narecco/software/R/niar')
```
This command will  just load the package (from my local repository on the CRG cluster).
For this step to work properly you must make sure all the  all the required dependencies have been already installed as just loading the package won’t install the dependencies. You can install most of the required packages (this step takes a while) with :

```R
install.packages(c('devtools', 'matrixStats', 'ggplot2', 'ggrepel', 'scales', 'patchwork',
                   'dplyr', 'tidyr', 'forcats', 'stringr', 'BiocManager', 'XICOR')) 
```
and the following Bioconductor packages:

```R
BiocManager::install("Biostrings")
```

If you get an error with this method it’s probably cause by the fact that the dependencies are already in the `.libPaths()`. Remember that once you log out of your `R` session the `niar` functions won’t be availbale anymore and you’ll need to load them again next time. 

#### Troubleshooting plots in RStudio server

In order to visualise the plots you might need to select the right graphics device, especially if you get an error that says something like:

```R
Error in diff.default(from) : 
  Shadow graphics device error: r error 4 (R code execution error)
In grDevices:::png("/tmp/Rtmp....",  :
  unable to open connection to X11 display ''
```

To solve this go to: Tools Menu (on top of the windos) > Global Options > General section > Graphics tab > and in the Graphic Device Backend drop down menu select `Cairo` and then click “Apply”. Now the plots should be correctly displayed.

## Quick Start

Currently this package contains:

- one function to perform Principal Component Analysis (PCA)
- some functions to fetch and parse data analysed with [vast-tools](https://github.com/vastgroup/vast-tools) for alternatively spliced events and gene expression with some plot functions to quicly glimpse into the data (e.g. `plot_corr_gene_expr_psi()`).
- [Biomart](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) handy functions
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) handy functions

#### PCA

The easiest way to make a PCA assuming `mat` is your numerical matrix is:

```R
showme_PCA2D(mat)
```

To know more you can type:

```R
?showme_PCA2D()
```

The underlying function is `prcomp` and you can pass extra arguments with `...` for example:

```R
showme_PCA2D(mat, scale. = T, center = F)
```

Extra info can be added from a dataframe  ( `mt`) which serves as a metadata. To specify which column of the dataframe contains the `colnames` of the matrix `mat` use `mcol`. In the following example the `mt` contains a column called `sample_name`:

```R
showme_PCA2D(mat = mat, mt = mt, mcol = "sample_name", show_variance = T, show_stats = T)
```

To show the PCA loadings:

```R
showme_PCA2D(mat = mat, n_loadings = 12)
```

More details can be found in the vignette linked below.

#### vast-tools

Since I use `vast-tools` quite often I made a couple of functions to easily import the text output tables into `R`.   Namely, `grep_psi()`  or `grep_gene_expression()` import the PSI of an AS events or gene expression levels that can then be turned into a long-format dataframe with the accompanying tidy functions `tidy_vst_psi()` or `tidy_vst_expr()`. These functions work great with the `magrittr` pipe (`%>%`) or the base `R` pipe operator (`|>`) as in:

```R
grep_psi(inclusion_tbl = file.path(dir_location, "INCLUSION_LEVELS_FULL-hg38-n-v251.tab"), 
         vst_id = c("HsaEX0000001", "HsaEX0000002"), tmp_dir = tempdir()) %>%
    tidy_vst_psi() 
```

which will return a dataframe. These functions are basically “hacks” that call the system `grep` command write to a temporary file that is then read into R and removed from the system. A better way would  probably be to implement the functions in `Rcpp`.

## Vignettes

[Link for PCA](https://htmlpreview.github.io/?https://github.com/Ni-Ar/niar/blob/main/doc/Introduction_Dim_Reduction.html).
