# My personal `R` package

Here, I collect the R function that I’ve created over time. This is not really stable or robust, but it’s easy to use, correct, and useful.

Feedback is welcome. This is my first `R` package.

If interested in giving it a try:

```R
devtools::install_github("Ni-Ar/niar")
```

Currently there’s only one function to plot PCA. The easiest way to make a PCR assuming `mat` is your numerical matrix:

```R
showme_PCA2D(mat)
```

The underlying function is `prcomp` and you can pass extra arguments with `...` for example:

```R
showme_PCA2D(mat, x = 2, scale. = T, center = F)
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

[Link for PCA](https://htmlpreview.github.io/?https://github.com/Ni-Ar/niar/blob/main/doc/Introduction.html).
