% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_environment.R
\name{load_conda_env}
\alias{load_conda_env}
\title{Load the input \code{conda} environment.}
\usage{
load_conda_env(envname = "bascule-env")
}
\description{
Function to load the input \code{conda} environment. The function
will raise an error if a \code{Python} version has already been attached to the
\code{reticulate} package.
In that case, it will be necessary to restart the \code{R} session and to load the
desired environment **before** calling \code{lineaGT} function interfacing with
\code{Python} - \code{filter_dataset() and fit()}.
}
