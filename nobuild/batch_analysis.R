library(dplyr)

ttypes = basilica::Degasperi_SBS_GEL_PCAWG_HW %>% names

cache = TRUE

runs = easypar::run(
  FUN = function(cohort){

    cli::cli_h1(cohort)

    input = basilica::Degasperi_SBS_GEL_PCAWG_HW[[cohort]] %>%
      dplyr::select(
        -sample,
        -organ,
        -cohort,
        -usedInCommonExtractionSBS,
        -usedInCommonExtractionDBS
        )

    input %>% print()

    dir.create(cohort)

    # Outputs
    output_fc = paste0(cohort, '/fit_COSMIC.rds')
    output_fd = paste0(cohort, '/fit_degasperi.rds')

    if(!file.exists(output_fc))
    {
      fit_COSMIC = basilica::fit(
        x = input,
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

      saveRDS(fit_COSMIC, file = output_fc)
    }
    else
      {
        cli::boxx(paste(output_fc, " -- CACHED!"), col = 'white', background_col = 'red') %>% cat
      }

    if(!file.exists(output_fd))
    {
      fit_degasperi = basilica::fit(
        x = input,
        reference_catalogue = basilica::Degasperi_catalogue,
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

      saveRDS(fit_degasperi, file = output_fd)
    }
    else
    {
      cli::boxx(paste(output_fd, " -- CACHED!"), col = 'white', background_col = 'red') %>% cat
    }

    return(0)
  },
  PARAMS = lapply(list, ttypes),
  parallel = FALSE,
  filter_errors = FALSE
)

saveRDS(runs, file = 'easypar-log.rds')
