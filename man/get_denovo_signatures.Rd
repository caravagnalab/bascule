% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getters.R
\name{get_denovo_signatures}
\alias{get_denovo_signatures}
\title{Get de novo (reference) signatures}
\usage{
get_denovo_signatures(x, types = get_types(x), matrix = FALSE)
}
\arguments{
\item{x}{bascule object.}

\item{types}{List of variant types to retrieve signames for.}

\item{matrix}{Logical. If `TRUE`, it will return the signatures in wide format.}
}
\value{
De novo signature matrix in long or wide format names with names equal to `types`.
}
\description{
Get de novo (reference) signatures
}
