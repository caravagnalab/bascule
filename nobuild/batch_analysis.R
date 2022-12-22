library(dplyr)
library(ggplot2)

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
      fit_COSMIC = fit(
        x = input,
        reference_catalogue = basilica::COSMIC_catalogue,
        input_catalogue = basilica::COSMIC_catalogue["SBS1", ],
        cohort = cohort,
        k = 1:5,
        lr = 0.01,
        steps = 1000,
        max_iterations = 30,
        phi = 0.15,
        blacklist = 'freq',
        delta = 0.7,
        groups=NULL,
        lambda_rate=NULL,
        sigma=FALSE,
        CUDA = FALSE,
        compile = TRUE
      )

      saveRDS(fit_COSMIC, file = output_fc)

      ggsave(
        filename = paste0(cohort, '/cosmic_signatures.pdf'),
        fit_COSMIC %>% plot_signatures(),
        width = 15,
        height = fit_COSMIC$n_catalogue + fit_COSMIC$n_denovo
      )

      ggsave(
        filename = paste0(cohort, '/cosmic_exposure.pdf'),
        fit_COSMIC %>% plot_exposure(),
        width = fit_COSMIC$n_samples * 0.1,
        height = 6,
        limitsize = FALSE
      )

      ggsave(
        filename = paste0(cohort, '/cosmic_exposure.pdf'),
        fit_COSMIC %>% plot_exposure(),
        width = fit_COSMIC$n_samples * 0.1,
        height = 6,
        limitsize = FALSE
      )

      ggsave(
        filename = paste0(cohort, '/cosmic_prevalence.pdf'),
        fit_COSMIC %>% plot_cohort_prevalence(),
        width = 8,
        height = 6,
        limitsize = FALSE
      )

      ggsave(
        filename = paste0(cohort, '/cosmic_reference.pdf'),
        fit_COSMIC %>% plot_similarity_reference(),
        width = 20,
        height = 20,
        limitsize = FALSE
      )
    }
    else
      {
        cli::boxx(paste(output_fc, " -- CACHED!"), col = 'white', background_col = 'red') %>% cat
      }

    if(!file.exists(output_fd))
    {
      fit_degasperi = fit(
        x = input,
        reference_catalogue = basilica::Degasperi_catalogue,
        cohort = cohort,
        k = 1:5,
        lr = 0.01,
        steps = 1000,
        max_iterations = 30,
        phi = 0.15,
        blacklist = 'freq',
        delta = 0.7,
        groups=NULL,
        input_catalogue=NULL,
        lambda_rate=NULL,
        sigma=FALSE
      )

      saveRDS(fit_degasperi, file = output_fd)

      ggsave(
        filename = paste0(cohort, '/degasperi_signatures.pdf'),
        fit_degasperi %>% plot_signatures(),
        width = 15,
        height = fit_COSMIC$n_catalogue + fit_COSMIC$n_denovo
      )

      ggsave(
        filename = paste0(cohort, '/degasperi_exposure.pdf'),
        fit_degasperi %>% plot_exposure(),
        width = fit_COSMIC$n_samples * 0.1,
        height = 6,
        limitsize = FALSE
      )

      ggsave(
        filename = paste0(cohort, '/degasperi_exposure.pdf'),
        fit_degasperi %>% plot_exposure(),
        width = fit_COSMIC$n_samples * 0.1,
        height = 6,
        limitsize = FALSE
      )

      ggsave(
        filename = paste0(cohort, '/degasperi_prevalence.pdf'),
        fit_degasperi %>% plot_cohort_prevalence(),
        width = 11,
        height = 6,
        limitsize = FALSE
      )

      ggsave(
        filename = paste0(cohort, '/degasperi_reference.pdf'),
        fit_degasperi %>% plot_similarity_reference(),
        width = 20,
        height = 30,
        limitsize = FALSE
      )
    }
    else
    {
      cli::boxx(paste(output_fd, " -- CACHED!"), col = 'white', background_col = 'red') %>% cat
    }

    fit_COSMIC = readRDS(file = output_fc)
    fit_degasperi = readRDS(file = output_fd)

    comparisons = plot_compare_fits(fit_COSMIC, fit_degasperi)

    pdf(paste0(cohort, '/comparison.pdf'), width = 20, height = 10)
    lapply(comparisons, print)
    dev.off()

    return(0)
  },
  PARAMS = lapply(ttypes, list),
  parallel = FALSE,
  filter_errors = FALSE
)

saveRDS(runs, file = 'easypar-log.rds')
