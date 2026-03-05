# Fit clustering

Fit clustering

## Usage

``` r
fit_clustering(
  x,
  cluster,
  hyperparameters = NULL,
  lr = 0.005,
  optim_gamma = 0.1,
  n_steps = 3000,
  py = NULL,
  enumer = "parallel",
  nonparametric = TRUE,
  autoguide = TRUE,
  CUDA = TRUE,
  compile = FALSE,
  store_parameters = FALSE,
  store_fits = TRUE,
  seed_list = c(10)
)
```

## Arguments

- x:

  Bascule object with signatures deconvolution performed.

- cluster:

  Maximum number of clusters.

- hyperparameters:

  List of hyperparameters passed to the NMF and clustering models.

- lr:

  Learning rate for SVI optimizer.

- optim_gamma:

  Deprecated.

- n_steps:

  Number of steps for the inference.

- py:

  User-installed version of `pybascule` package

- enumer:

  Enumeration used for clustering (either \`parallel\` or
  \`sequential\`).

- nonparametric:

  Deprecated. The model only works in nonparametric way.

- autoguide:

  Logical. If \`TRUE\`, the clustering model will use the Pyro
  autoguide.

- CUDA:

  Logical. If \`TRUE\` and a GPU is available, the models will run on
  GPU.

- compile:

  Deprecated.

- store_parameters:

  Logical. If \`TRUE\`, parameters at every step of inference will be
  stored in the object.

- store_fits:

  Logical. If \`TRUE\`, all tested fits, i.e., for every value of \`K\`,
  will be stored in the object.

- seed_list:

  List of seeds used for every input configuration.

## Value

Bascule object.
