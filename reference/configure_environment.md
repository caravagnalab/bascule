# Configure the reticulate environment

Function to configure the Python dependencies on R. If a Python
environment is not available, the function will check if there is a
version of `conda` or `miniconda`, otherwise it will install
`miniconda`, on which install the Python package `pybascule`.

## Usage

``` r
configure_environment(envname = "bascule-env", use_default = FALSE)
```

## Arguments

- env_name:

  name of the `conda` environment to use, if available.
