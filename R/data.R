#' COSMIC catalogue for SBS (version 3.4, GRCh37)
#'
#' @docType data
#'
#' @usage data(COSMIC_sbs)
#'
#' @format add
#'
#' @keywords datasets

"COSMIC_sbs"


#' COSMIC catalogue for DBS (version 3.4, GRCh37)
#'
#' @docType data
#'
#' @usage data(COSMIC_dbs)
#'
#' @keywords datasets

"COSMIC_dbs"


#' COSMIC catalogue for Indels (version 3.4, GRCh37)
#'
#' @docType data
#'
#' @usage data(COSMIC_indels)
#'
#' @keywords datasets

"COSMIC_indels"


#' COSMIC catalogue for CN (version 3.4, GRCh37)
#'
#' @docType data
#'
#' @usage data(COSMIC_cn)
#'
#' @keywords datasets

"COSMIC_cn"


#' COSMIC catalogue for SBS filtered (version 3.4, GRCh37)
#'
#' Filtered version of the COSMIC catalogue. Only validated signatures are present.
#' For each signature (besides SBS5) the contexts with probability mass below 0.02 are set to 0
#' and the probabilities are re-normalized over contexts.
#'
#' @usage data(COSMIC_sbs_filt)
#'
#' @format add
#'
#' @keywords datasets

"COSMIC_sbs_filt"


#' Degasperi SBS catalogue
#'
#' SBS signatures extracted and released in "Substitution mutational signatures
#' in whole-genome–sequenced cancers in the UK population."
#' Degasperi et al.Science 376, Issue 6591, 2022. PMID: 35949260
#'
#' @docType data
#'
#' @usage data(Degasperi_SBS)
#'
#' @keywords datasets

"Degasperi_SBS"


#' Degasperi DBS catalogue
#'
#' DBS signatures extracted and released in "Substitution mutational signatures
#' in whole-genome–sequenced cancers in the UK population."
#' Degasperi et al.Science 376, Issue 6591, 2022. PMID: 35949260
#'
#' @docType data
#'
#' @usage data(Degasperi_DBS)
#'
#' @keywords datasets

"Degasperi_DBS"


#' Analysis of a synthetic cohort
#'
#' Input SBS and DBS counts and BASCULE fits for a cohort of 150 samples whose
#' counts and exposures have been generated from the BASCULE generative model.
#'
#' @docType data
#'
#' @usage data(synthetic_data)
#'
#' @keywords datasets

"synthetic_data"


#' Analysis of breast tumours
#'
#' Input SBS and DBS counts and BASCULE fits from 2682 breast tumours
#' released in "Substitution mutational signatures in whole-genome–sequenced
#' cancers in the UK population."
#' Degasperi et al.Science 376, Issue 6591, 2022. PMID: 35949260
#'
#' @docType data
#'
#' @usage data(breast_data)
#'
#' @keywords datasets

"breast_data"


#' Skin cohort BASCULE fit
#'
#' BASCULE fit for skin cohort data used in the survival analysis.
#'
#' @docType data
#'
#' @usage data(skin_fit)
#'
#' @keywords datasets

"skin_fit"


#' Skin cohort metadata
#'
#' Clinical data for the skin cancer cohort used in the survival analysis.
#'
#' @docType data
#'
#' @usage data(skin_metadata)
#'
#' @keywords datasets

"skin_metadata"
