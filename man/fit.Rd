% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{fit}
\alias{fit}
\title{Fit a bascule object}
\usage{
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
}
\arguments{
\item{counts}{List of mutation counts matrices from multiple variant types.}

\item{k_list}{List of number of denovo signatures to test.}

\item{cluster}{Maximum number of clusters. If `NULL`, no clustering will be performed.}

\item{reference_cat}{List of reference catalogues to use for NMF. Names must be the same as input counts.}

\item{keep_sigs}{List of reference signatures to keep even if found with low exposures.}

\item{hyperparameters}{List of hyperparameters passed to the NMF and clustering models.}

\item{lr}{Learning rate used for SVI.}

\item{optim_gamma}{Deprecated}

\item{n_steps}{Number of iterations for inference.}

\item{py}{User-installed version of \code{pybascule} package}

\item{enumer}{Enumeration used for clustering (either `parallel` or `sequential`).}

\item{nonparametric}{Deprecated. The model only works in nonparametric way.}

\item{autoguide}{Logical. If `TRUE`, the clustering model will use the Pyro autoguide.}

\item{filter_dn}{Logical. If `TRUE`, all contexts below 0.01 in denovo signatures will be set to 0, provided the filtered signatures remain consistent with the inferred ones.}

\item{min_exposure}{Reference signatures with an exposures lower than `min_exposure` will be dropped.}

\item{CUDA}{Logical. If `TRUE` and a GPU is available, the models will run on GPU.}

\item{compile}{Deprecated.}

\item{store_parameters}{Logical. If `TRUE`, parameters at every step of inference will be stored in the object.}

\item{store_fits}{Logical. If `TRUE`, all tested fits, i.e., for every value of `K`, will be stored in the object.}

\item{seed_list}{List of seeds used for every input configuration.}
}
\value{
Bascule object.
}
\description{
Fit a bascule object
}
