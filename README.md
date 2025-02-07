# Nice Interactive Analysis in `R` package

Here I collect many `R` functions that I’ve created over time. The primary purpose of this package is to neatly organise my R code in one place and have it ready to use whenever possible. Plus this forces me to document my functions and it's easier to share the analysis code. The code could be subject to frequent changes and is "always" in development. Feedback is welcome!
If you find an error please [open a GitHub issue](https://github.com/Ni-Ar/niar/issues/new).

## Requirements
 - `R` >=4.1.0
 - Mac or Linux operating system (not tested on Windows)

## Installation

```R
devtools::install_github("Ni-Ar/niar")
```

####  Dependencies

If you encounter installation issues try installing the R package dependencies before with: 
```R
install.packages(c('devtools', 'matrixStats', 'BiocManager', 'XICOR', 
                   'ggplot2', 'ggrepel', 'scales', 'patchwork',
                   'MetBrewer', 'ggalluvial', 'ggfittext', 'ggseqlogo', 'seqinr',
                   'dplyr', 'tidyr', 'tibble', 'forcats', 'stringr'))
```
That step takes a while.

```R
install.packages('Cairo') 
```

If you get an error installing `Cairo` you might need to first install the cairographics C library on your operating system from [here](https://www.cairographics.org/download/).

The following Bioconductor packages:

```R
BiocManager::install("Biostrings")
BiocManager::install("biomaRt")
BiocManager::install("DESeq2")
BiocManager::install("csaw")
BiocManager::install("msa")
```

#### Troubleshooting plots in RStudio server

To visualise the plots you might need to select the right graphics device, especially if you get an error that says something like:
```R
Error in diff.default(from) : 
  Shadow graphics device error: r error 4 (R code execution error)
In grDevices:::png("/tmp/Rtmp....",  :
  unable to open connection to X11 display ''
```

To solve this, go to the Tools Menu (at the top of the window) > Global Options > General section > Graphics tab > and select Cairo from the Graphic Device Backend drop-down menu. Then click “Apply.” The plots should now be correctly displayed.

## Overview

Currently, this package contains:
- one function to perform Principal Component Analysis (PCA) in 2D with lots of options to enrich visualisation and exploration. See the [vignette](#Vignettes) below for more details.
- several functions to fetch and parse data analysed with [vast-tools](https://github.com/vastgroup/vast-tools) for alternatively spliced events and gene expression. There are also plotting functions to quickly glimpse into the data (e.g. `plot_corr_gene_expr_psi()`).
- Some publically available datasets have been packaged in *ad-hoc* functions to quickly plot and explore the data:
	- [Mouse Development ENCODE](https://www.encodeproject.org/mouse-development-matrix/?type=Experiment&status=released&related_series.@type=OrganismDevelopmentSeries&replicates.library.biosample.organism.scientific_name=Mus+musculus) data `plot_mouse_tissue_devel()` which uses data I preprocessed fetched with `get_mouse_tissue_devel_tbl()`. See the [vignette](#Vignettes) below for more details.
- Some [Biomart](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) handy functions for quick gene IDs conversions (e.g. `ensembl_id_2_gene_name()`).
- Some [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) wrappers.
- Some [rMATS](https://github.com/Xinglab/rmats-turbo) wrappers.
- [EpiProfile](https://github.com/zfyuan/EpiProfile2.0_Family) (Histone Mass Spectrometry) post-processing functions.
- Multiple sequence analysis from fasta format to generate PWMs, visualise logos or Jensen-Shannon divergence.
- Generic genomics files handling (e.g. bed, gtf).

More examples grouped by topic are listed below:

### PCA

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

Extra info can be passed from a metadata dataframe with  `mt =`. To specify which column of the dataframe contains the `colnames` of the matrix `mat` use `mcol`. In the following example the `mt` contains a column called `sample_name`:
```R
showme_PCA2D(mat = mat, mt = mt, mcol = "sample_name", show_variance = T, show_stats = T)
```

To show the PCA loadings:
```R
showme_PCA2D(mat = mat, n_loadings = 12)
```

More details can be found in the [vignette](#Vignettes) below.

### vast-tools

Since I use `vast-tools` quite often I made functions to easily import the output tables into `R`. Namely, `grep_psi()` or `grep_gene_expression()` import the PSI of an AS events or gene expression levels respectively and parse the data into a long-format dataframe with the accompanying tidy functions `tidy_vst_psi()` or `tidy_vst_expr()`. These functions work great with the `magrittr` pipe (`%>%`) or the base `R` pipe operator (`|>`) as in:
```R
grep_psi(inclusion_tbl = file.path(dir_location, "INCLUSION_LEVELS_FULL-hg38-n-v251.tab"), 
         vst_id = c("HsaEX0000001", "HsaEX0000002")) |>
    tidy_vst_psi() 
```

These functions are basically “hacks” that call the system `grep` command, and write to a temporary file that is then read into R and removed from the system. Maybe a better way would probably be to implement the functions in `Rcpp`.

## Vignettes

[Link for PCA](https://htmlpreview.github.io/?https://github.com/Ni-Ar/niar/blob/main/doc/Introduction_Dim_Reduction.html).

[Link for mouse ENCODE AS data AS exploration](https://htmlpreview.github.io/?https://github.com/Ni-Ar/niar/blob/main/doc/Mouse_Development_ENCODE_AS_exploration.html).


## To do

- [ ] Make a vignette for `biomaRt` functions 
- [ ] Make a vignette for vast-tools utility and plotting functions, especially correlations. 
- [ ] Maybe add the mouse ENCODE data (fetched with `get_mouse_tissue_devel_tbl`) to the package?
- [ ] Make a vignette for logo analysis
