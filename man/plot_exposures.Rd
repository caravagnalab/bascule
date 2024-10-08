% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_exposures.R
\name{plot_exposures}
\alias{plot_exposures}
\title{Function to visualize the estimated exposures}
\usage{
plot_exposures(
  x,
  types = get_types(x),
  samples = get_samples(x),
  clusters = get_cluster_labels(x),
  sample_name = FALSE,
  color_palette = NULL,
  add_centroid = FALSE,
  sort_by = NULL,
  exposure_thr = 0,
  quantile_thr = 0,
  signatures_list = get_signames(x)
)
}
\arguments{
\item{x}{bascule object.}

\item{types}{List of variant types to visualize.}

\item{samples}{List of samples to visualize.}

\item{clusters}{List of clusters to visualize.}

\item{sample_name}{Logical. If `TRUE`, sample names will be reported on the x axis.}

\item{color_palette}{Custom color palette for signatures.}

\item{add_centroid}{Logical. If `TRUE`, also clustering's centroids will be plotted.}

\item{sort_by}{Signature to sort patients' exposures by.}

\item{exposure_thr}{Only signatures with exposures greater than `exposure_thr` in all samples will be highlighted.}

\item{quantile_thr}{add}

\item{signatures_list}{add}
}
\value{
ggplot object.
}
\description{
Function to visualize the estimated exposures
}
