# bascule

### Installing the package

The `bascule` R package can be installed directly from GitHub using
`devtools`:

``` r
devtools::install_github("caravagnalab/bascule")
```

Once installed, load the package in your R session:

``` r
library(bascule) 
```

### Installing Python dependencies

The `bascule` package relies on the Python module `pybascule`. The
recommended approach is to install these Python dependencies inside a
dedicated conda environment.

This can be done from R using the `reticulate` package through the
following steps.

- if Miniconda is not already available on your system, you can install
  it using through
  [`reticulate`](https://rstudio.github.io/reticulate/reference/install_miniconda.html):

``` r
install_miniconda(path = miniconda_path(), update = TRUE, force = FALSE)
```

- create a new conda environment with Python version 3.10:

``` r
reticulate::conda_create(envname="bascule", python_version="3.10")
```

- activate the new conda environment within the current R session:

``` r
reticulate::use_condaenv("bascule")
```

- install the `pybascule` Python module using pip inside the conda
  environment:

``` r
reticulate::conda_install(envname="bascule", packages="pybascule", pip=TRUE)
```

## Fitting a `bascule` object

Before fitting a `bascule` object with the
[`fit()`](https:%3A/caravagnalab.github.io/bascule/reference/fit.md) or
the
[`fit_clustering()`](https:%3A/caravagnalab.github.io/bascule/reference/fit_clustering.md)
functions, make sure that the correct conda environment is active so
that the `pybascule` module can be accessed:

``` r
reticulate::use_condaenv("bascule")
```

Functions that rely on `pybascule` include an additional argument, `py`,
which allows you to explicitly pass the imported Python module. This
ensures that the correct version of `pybascule` is used.

``` r
reticulate::use_condaenv("bascule")
py = reticulate::import("pybascule")
x = fit(..., py=py)
# or
x = fit_clustering(..., py=py)
```
