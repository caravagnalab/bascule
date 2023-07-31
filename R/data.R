#' COSMIC catalogue (version XXXX)
#'
#' @docType data
#'
#' @usage data(COSMIC_catalogue)
#'
#' @format add
#'
#' @keywords datasets

"COSMIC_catalogue"

#' COSMIC catalogue (version XXXX)
#'
#' Filtered version of the COSMIC catalogue. Only validated signatures are present.
#' For each signature (besides SBS5) the contexts with probability mass below 0.02 are set to 0
#' and the probabilities are re-normalized over contexts.
#'
#' @usage data(COSMIC_filt)
#'
#' @format add
#'
#' @keywords datasets

"COSMIC_filt"

#' COSMIC catalogue (version XXXX)
#'
#' Filtered version of the COSMIC catalogue. Only validated signatures are present.
#' For each signature (besides SBS5) the contexts with probability mass below 0.02 are set to 0
#' and the probabilities are re-normalized over contexts.
#' Moreover, all signatures with high cosine similarity are merged.
#'
#' @usage data(COSMIC_filt_merged)
#'
#' @format add
#'
#' @keywords datasets

"COSMIC_filt_merged"

#' Degasperi catalogue (version XXXX)
#'
#' Description
#'
#' @docType data
#'
#' @usage data(Degasperi_catalogue)
#'
#' @format add
#'
#' @keywords datasets

"Degasperi_catalogue"

#' SBS data for PMID: 35949260
#'
#' SBS data released with the paper "Substitution mutational signatures in
#' whole-genome–sequenced cancers in the UK population."
#' Degasperi et al.Science 376, Issue 6591, 2022. PMID: 35949260
#'
#' @docType data
#'
#' @usage data(Degasperi_SBS_GEL_PCAWG_HW)
#'
#' @keywords datasets

"Degasperi_SBS_GEL_PCAWG_HW"

#' Example input SBS data (synthetic)
#'
#' Description
#'
#' @docType data
#'
#' @usage data(example_input)
#'
#' @keywords datasets

"example_input"
