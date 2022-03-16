# My personal `R` package

Here, I collect R functions that I’ve created over time. This is not really stable or robust, the primary purpose of this package is to neatly organise my R functions in one place and have them ready to use in an easy way. Feedback is welcome as this is my first `R` package.

## Installation 

### From this GitHub Repository

If interested in giving it a try:

```R
devtools::install_github("Ni-Ar/niar")
```

### Load the local repository (only for CRG users)

If you have a CRG user account, instead of installing, you can load the R package at the beginning of your session. Log in with your credentials on the [CRG RStudio Server IDE](http://rstudio4.linux.crg.es/) (R version 4.0.3). Install all the required dependencies as just loading the package won’t do that for you. This step takes a while.

```R
install.packages( c('devtools', 'matrixStats','ggplot2','ggrepel','scales',
                   'patchwork','dplyr','tidyr','forcats') ) 
```

Then just load the package (from my local repository on the CRG cluster) with:

```R
devtools::load_all(path = '/users/mirimia/narecco/software/R/niar')
```

Which should not throw you an error if all the dependencies are already in the `.libPaths()`. Once you log out of your R session the functions of the R package won’t be availbale anymore.

## Examples

Currently there’s only one function to plot PCA. The easiest way to make a PCR assuming `mat` is your numerical matrix is:

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

Extra info can be added from a dataframe  ( `mt`) which serves as a metadata. To specify which column of the dataframe contains the `colnames` of the matrix `mat` use `mcol`. In the following example the `mt` contains a column called `sample_name`

```R
showme_PCA2D(mat = mat,  mt = mt, mcol = "sample_name",
             show_variance = T, show_stats = T)
```

To show the PCA loadings:

```R
showme_PCA2D(mat = mat, n_loadings = 12)
```

## Vignettes

[Link for PCA](https://htmlpreview.github.io/?https://github.com/Ni-Ar/niar/blob/main/doc/Introduction_Dim_Reduction.html).
