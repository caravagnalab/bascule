require(dplyr)

# # reticulate::use_condaenv('basilica-env')
# # reticulate::py_discover_config('basilica-env')
# # reticulate::py_list_packages('basilica-env')
# #
# # devtools::install()
#
# Degasperi_SBS_GEL_PCAWG_HW = lapply(Degasperi_SBS_GEL_PCAWG_HW, function(x){
#   x %>% select(sample, organ, cohort, usedInCommonExtractionSBS, usedInCommonExtractionDBS, everything())
# })
#
# usethis::use_data(Degasperi_SBS_GEL_PCAWG_HW, overwrite = TRUE)

M = basilica::Degasperi_SBS_GEL_PCAWG_HW[['Myeloid']]
M_matrix = M %>% select(-sample, -organ, -cohort, -usedInCommonExtractionSBS,
                        -usedInCommonExtractionDBS)

M_matrix = basilica::example_input

x = fit(
    x = M_matrix,
    reference_catalogue = basilica::COSMIC_catalogue,
    k = 7,
    lr = 0.01,
    steps = 1000,
    phi = 0.05,
    delta = 0.9,
    groups=NULL,
    input_catalogue=NULL,
    lambda_rate=NULL,
    sigma=FALSE
)

# export(fit)
# export(get_catalogue_signatures)
# export(get_denovo_signatures)
# export(get_exposure)
# export(plot_exposure)
# export(plot_signatures)

x$S9$X1 %>% length()
Degasperi_catalogue = x$S21
rownames(Degasperi_catalogue) = Degasperi_catalogue$mutationClass
Degasperi_catalogue = Degasperi_catalogue[, -1]
Degasperi_catalogue = Degasperi_catalogue %>% t
Degasperi_catalogue = apply(Degasperi_catalogue, c(1,2), as.numeric)

(basilica::COSMIC_catalogue %>% colnames) == (Degasperi_catalogue %>% colnames)

Degasperi_catalogue = Degasperi_catalogue[, (basilica::COSMIC_catalogue %>% colnames)]

usethis::use_data(Degasperi_catalogue)
