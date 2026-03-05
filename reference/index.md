# Package index

## Perform inference

Functions to fit NMF and clustering.

- [`fit()`](https:%3A/caravagnalab.github.io/bascule/reference/fit.md) :
  Fit a bascule object
- [`fit_clustering()`](https:%3A/caravagnalab.github.io/bascule/reference/fit_clustering.md)
  : Fit clustering
- [`refine_denovo_signatures()`](https:%3A/caravagnalab.github.io/bascule/reference/refine_denovo_signatures.md)
  : De novo signatures refinement.
- [`merge_clusters()`](https:%3A/caravagnalab.github.io/bascule/reference/merge_clusters.md)
  : Function to merge similar clusters

## Getters

Utility functions to extract attributes.

- [`get_input()`](https:%3A/caravagnalab.github.io/bascule/reference/get_input.md)
  : Get input data
- [`get_signames()`](https:%3A/caravagnalab.github.io/bascule/reference/get_signames.md)
  : Get signatures names
- [`get_fixed_signames()`](https:%3A/caravagnalab.github.io/bascule/reference/get_fixed_signames.md)
  : Get reference signatures names
- [`get_denovo_signames()`](https:%3A/caravagnalab.github.io/bascule/reference/get_denovo_signames.md)
  : Get denovo signatures names
- [`get_signatures()`](https:%3A/caravagnalab.github.io/bascule/reference/get_signatures.md)
  : Get signatures
- [`get_fixed_signatures()`](https:%3A/caravagnalab.github.io/bascule/reference/get_fixed_signatures.md)
  : Get fixed (reference) signatures
- [`get_denovo_signatures()`](https:%3A/caravagnalab.github.io/bascule/reference/get_denovo_signatures.md)
  : Get de novo (reference) signatures
- [`get_exposure()`](https:%3A/caravagnalab.github.io/bascule/reference/get_exposure.md)
  : Get exposures
- [`get_n_denovo()`](https:%3A/caravagnalab.github.io/bascule/reference/get_n_denovo.md)
  : Get number of de novo signatures
- [`get_n_groups()`](https:%3A/caravagnalab.github.io/bascule/reference/get_n_groups.md)
  : Get number of groups
- [`convert_dn_names()`](https:%3A/caravagnalab.github.io/bascule/reference/convert_dn_names.md)
  : Function to map de novo signatures to input catalogues

## Plotting functions

Plotting functions available in the package

- [`plot_data()`](https:%3A/caravagnalab.github.io/bascule/reference/plot_data.md)
  : Function to visualize the input and reconstructed mutation counts.
- [`plot_signatures()`](https:%3A/caravagnalab.github.io/bascule/reference/plot_signatures.md)
  : Function to visualize the estimated and reference signatures.
- [`plot_exposures()`](https:%3A/caravagnalab.github.io/bascule/reference/plot_exposures.md)
  : Function to visualize the estimated exposures
- [`plot_centroids()`](https:%3A/caravagnalab.github.io/bascule/reference/plot_centroids.md)
  : Plot clustering centroids
- [`plot_fit()`](https:%3A/caravagnalab.github.io/bascule/reference/plot_fit.md)
  : Report of the fit
- [`plot_similarity_reference()`](https:%3A/caravagnalab.github.io/bascule/reference/plot_similarity_reference.md)
  : Function to visualize the similarity of denovo and reference
  signatures.
- [`plot_scores()`](https:%3A/caravagnalab.github.io/bascule/reference/plot_scores.md)
  : Plot model selection scores
- [`plot_gradient_norms()`](https:%3A/caravagnalab.github.io/bascule/reference/plot_gradient_norms.md)
  : Plot parameters gradients norms
- [`plot_posterior_probs()`](https:%3A/caravagnalab.github.io/bascule/reference/plot_posterior_probs.md)
  : Plot posterior probabilities

## S3 objects

Print and plot S3 functions

- [`print(`*`<bascule_obj>`*`)`](https:%3A/caravagnalab.github.io/bascule/reference/print.bascule_obj.md)
  :

  Print for class `'bascule_obj'`.

- [`plot(`*`<bascule_obj>`*`)`](https:%3A/caravagnalab.github.io/bascule/reference/plot.bascule_obj.md)
  :

  Plot for class `'bascule_obj'`.

## Data

Description of datasets available in the package.

- [`COSMIC_sbs`](https:%3A/caravagnalab.github.io/bascule/reference/COSMIC_sbs.md)
  : COSMIC catalogue for SBS (version 3.4, GRCh37)
- [`COSMIC_sbs_filt`](https:%3A/caravagnalab.github.io/bascule/reference/COSMIC_sbs_filt.md)
  : COSMIC catalogue for SBS filtered (version 3.4, GRCh37)
- [`COSMIC_dbs`](https:%3A/caravagnalab.github.io/bascule/reference/COSMIC_dbs.md)
  : COSMIC catalogue for DBS (version 3.4, GRCh37)
- [`COSMIC_indels`](https:%3A/caravagnalab.github.io/bascule/reference/COSMIC_indels.md)
  : COSMIC catalogue for Indels (version 3.4, GRCh37)
- [`COSMIC_cn`](https:%3A/caravagnalab.github.io/bascule/reference/COSMIC_cn.md)
  : COSMIC catalogue for CN (version 3.4, GRCh37)
- [`Degasperi_SBS`](https:%3A/caravagnalab.github.io/bascule/reference/Degasperi_SBS.md)
  : Degasperi SBS catalogue
- [`Degasperi_DBS`](https:%3A/caravagnalab.github.io/bascule/reference/Degasperi_DBS.md)
  : Degasperi DBS catalogue
- [`synthetic_data`](https:%3A/caravagnalab.github.io/bascule/reference/synthetic_data.md)
  : Analysis of a synthetic cohort
- [`breast_data`](https:%3A/caravagnalab.github.io/bascule/reference/breast_data.md)
  : Analysis of breast tumours
- [`skin_metadata`](https:%3A/caravagnalab.github.io/bascule/reference/skin_metadata.md)
  : Skin cohort metadata
- [`skin_fit`](https:%3A/caravagnalab.github.io/bascule/reference/skin_fit.md)
  : Skin cohort BASCULE fit

## Configure conda environment

Helper functions to create and check dependencies.

- [`configure_environment()`](https:%3A/caravagnalab.github.io/bascule/reference/configure_environment.md)
  : Configure the reticulate environment

- [`have_loaded_env()`](https:%3A/caravagnalab.github.io/bascule/reference/have_loaded_env.md)
  :

  Check if there is a loaded `conda` environment.

- [`have_python_deps()`](https:%3A/caravagnalab.github.io/bascule/reference/have_python_deps.md)
  :

  Check if `Python` packages are installed in the environment.

- [`load_conda_env()`](https:%3A/caravagnalab.github.io/bascule/reference/load_conda_env.md)
  :

  Load the input `conda` environment.

- [`which_conda_env()`](https:%3A/caravagnalab.github.io/bascule/reference/which_conda_env.md)
  : Retrieve the name of the currently loaded environment.
