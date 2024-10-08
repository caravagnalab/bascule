% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/configure_environment.R
\name{configure_environment}
\alias{configure_environment}
\title{Configure the reticulate environment}
\usage{
configure_environment(envname = "bascule-env", use_default = FALSE)
}
\arguments{
\item{env_name}{name of the \code{conda} environment to use, if available.}
}
\description{
Function to configure the Python dependencies on R.
If a Python environment is not available, the function will check if there is a version of
\code{conda} or \code{miniconda}, otherwise it will install \code{miniconda}, on which
install the Python package \code{pybascule}.
}
