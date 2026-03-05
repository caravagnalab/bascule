# Load the input `conda` environment.

Function to load the input `conda` environment. The function will raise
an error if a `Python` version has already been attached to the
`reticulate` package. In that case, it will be necessary to restart the
`R` session and to load the desired environment \*\*before\*\* calling
`lineaGT` function interfacing with `Python` -
`filter_dataset() and fit()`.

## Usage

``` r
load_conda_env(envname = "bascule-env")
```
