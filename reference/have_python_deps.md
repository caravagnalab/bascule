# Check if `Python` packages are installed in the environment.

Function to check if one or more `Python` packages are present in a
`conda` environment.

## Usage

``` r
have_python_deps(envname = "", py_pkgs = c("pybascule"))
```

## Arguments

- envname:

  the name of the environment to check. If empty, the function will
  check the currently loaded environment.

- py_pkgs:

  a list or vector of `Python` packages.

## Value

a list of Boolean. For each input package, `TRUE` if the package is
installed, `FALSE` otherwise.
