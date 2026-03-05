# bascule

### Installation of the package

You can install bascule from GitHub using `devtools`.

``` r
devtools::install_github("caravagnalab/bascule")
```

Load the package.

``` r
library(bascule) 
```

### Python dependencies installation

When the package is loaded, you might execute the
[`configure_environment()`](https:%3A/caravagnalab.github.io/bascule/reference/configure_environment.md)
function to automatically check whether:

- a version of Anaconda or Miniconda is available, otherwise a Miniconda
  installation will be started,

- if there is no `conda` environment loaded, the package will check if
  the `bascule-env` is present, otherwise it will be created,

- to use an existing `conda` environment, it can be loaded **before**
  loading the package, either through the `reticulate` function
  [`reticulate::use_condaenv()`](https://rstudio.github.io/reticulate/reference/use_python.html)
  or by using the `bascule` function
  [`load_conda_env()`](https:%3A/caravagnalab.github.io/bascule/reference/load_conda_env.md):

``` r
reticulate::use_condaenv("env-name", required=TRUE)
load_conda_env(envname="env-name")
```

- eventually, if the required `Python` dependencies are not installed in
  the loaded environment, they will be installed.

### Functions to manually configure an environment

The function
[`configure_environment()`](https:%3A/caravagnalab.github.io/bascule/reference/configure_environment.md)
can be used interactively to manually configure an existing environment,
or to create one from scratch.

``` r
configure_environment(env_name="bascule-env", use_default=F)
```

The function will first check if a Anaconda or Miniconda installation is
available, otherwise it will prompt a Miniconda installation. The input
name of the environment is either the name of an existing environment or
the name of the environment to be created. The environment will be
loaded or created, and the required `Python` dependencies will be
installed.

### Check the loaded Python version and environment

The package provides also a set of helper functions to check if an
environment is loaded.

- [`have_loaded_env()`](https:%3A/caravagnalab.github.io/bascule/reference/have_loaded_env.md)
  to check if an environment is already loaded,
- [`which_conda_env()`](https:%3A/caravagnalab.github.io/bascule/reference/which_conda_env.md)
  to check which environment is loaded.
- [`have_python_deps()`](https:%3A/caravagnalab.github.io/bascule/reference/have_python_deps.md)
  to check if a `Python` packages list is installed in the specified
  environment.
- [`load_conda_env()`](https:%3A/caravagnalab.github.io/bascule/reference/load_conda_env.md)
  to load the specified environment.
