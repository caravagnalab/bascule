# Fit a bascule object

Fit a bascule object

## Usage

``` r
fit(
  counts,
  k_list,
  cluster = NULL,
  reference_cat = list(SBS = COSMIC_filt, DBS = COSMIC_dbs),
  keep_sigs = c("SBS1", "SBS5"),
  hyperparameters = NULL,
  lr = 0.005,
  optim_gamma = 0.1,
  n_steps = 3000,
  py = NULL,
  enumer = "parallel",
  nonparametric = TRUE,
  autoguide = FALSE,
  filter_dn = FALSE,
  min_exposure = 0.2,
  CUDA = TRUE,
  compile = FALSE,
  store_parameters = FALSE,
  store_fits = TRUE,
  seed_list = c(10)
)
```

## Arguments

- counts:

  List of mutation counts matrices from multiple variant types.

- k_list:

  List of number of denovo signatures to test.

- cluster:

  Maximum number of clusters. If \`NULL\`, no clustering will be
  performed.

- reference_cat:

  List of reference catalogues to use for NMF. Names must be the same as
  input counts.

- keep_sigs:

  List of reference signatures to keep even if found with low exposures.

- hyperparameters:

  List of hyperparameters passed to the NMF and clustering models.

- lr:

  Learning rate used for SVI.

- optim_gamma:

  Deprecated

- n_steps:

  Number of iterations for inference.

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

- filter_dn:

  Logical. If \`TRUE\`, all contexts below 0.01 in denovo signatures
  will be set to 0, provided the filtered signatures remain consistent
  with the inferred ones.

- min_exposure:

  Reference signatures with an exposures lower than \`min_exposure\`
  will be dropped.

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
