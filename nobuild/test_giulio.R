require(dplyr)

reticulate::use_condaenv('basilica-env')
reticulate::py_discover_config('basilica-env')
reticulate::py_list_packages('basilica-env')

devtools::install()

Degasperi_SBS_GEL_PCAWG_HW = lapply(Degasperi_SBS_GEL_PCAWG_HW, function(x){
  x %>% select(sample, organ, cohort, usedInCommonExtractionSBS, usedInCommonExtractionDBS, everything())
})

usethis::use_data(Degasperi_SBS_GEL_PCAWG_HW, overwrite = TRUE)

M =  %>% [['Myeloid']]
M_matrix = M %>% select(-sample, -organ, -cohort, -usedInCommonExtractionSBS,
                        -usedInCommonExtractionDBS)

M_matrix = basilica::example_input

M_fit = basilica::fit(
    x = M_matrix,
    reference_catalogue = basilica::COSMIC_catalogue,
    k = 3,
    lr = 0.01,
    steps = 1000,
    phi = 0.025,
    delta = 0.025,
    groups=NULL,
    input_catalogue=NULL,
    lambda_rate=NULL,
    sigma=FALSE
)
