require(dplyr)

reticulate::use_condaenv('basilica-env')
reticulate::py_discover_config('basilica-env')
reticulate::py_list_packages('basilica-env')

devtools::install()

M = basilica::Degasperi_SBS_GEL_PCAWG_HW[['Myeloid']]
M_matrix = M %>% select(-sample, -organ, -cohort, -usedInCommonExtractionSBS,
                        -usedInCommonExtractionDBS)

M_matrix = readr::read_csv("~/Downloads/data_sigphylo.csv")
M_matrix = read.table("~/Downloads/data_sigphylo.csv", header = T, check.names = F, sep = ',')

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
