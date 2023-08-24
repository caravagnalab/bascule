# fit <- function(x,
#                 k,
#                 py = NULL,
#                 reference_catalogue = COSMIC_filtered,
#                 input_catalogue = COSMIC_filtered["SBS1", ],
#                 hyperparameters = NULL,
#                 filtered_cat = FALSE,
#                 cohort = "MyCohort",
#                 lr = 0.05,
#                 steps = 500,
#                 max_iterations = 20,
#                 blacklist = NULL,
#                 phi = 0.05,
#                 delta = 0.9,
#                 filt_pi = 0.1,
#                 groups = NULL,
#                 lambda_rate = NULL,
#                 sigma = FALSE,
#                 CUDA = FALSE,
#                 compile = FALSE,
#                 enforce_sparsity = FALSE,
#                 store_parameters = FALSE,
#                 regularizer="cosine",
#                 reg_weight = 1,
#                 reg_bic = TRUE,
#                 cosine_by_subs = FALSE,
#                 stage="",
#                 verbose = TRUE
# )
#
# {
#
#   sig_col = function(x) crayon::blue(x)
#
#   fit = list()
#   class(fit) = 'basilica_obj'
#
#   # cli::cli_h1("MUSICA - MUtational Signature Inference with a CAtalogue ")
#   cli::cli_h1("Basilica - Bayesian signature learning with a catalogue")
#   cat("\n")
#
#   # First, sanitize inputs
#   sanitized_inputs = sanitize_inputs(
#     x = x,
#     reference_catalogue = reference_catalogue,
#     k = k,
#     lr = lr,
#     steps = steps,
#     phi = phi,
#     delta = delta,
#     groups = NULL,
#     input_catalogue = input_catalogue,
#     lambda_rate = lambda_rate,
#     sigma = sigma,
#     blacklist = blacklist
#   )
#
#   x = sanitized_inputs$x
#   reference_catalogue = sanitized_inputs$reference_catalogue
#   input_catalogue = sanitized_inputs$input_catalogue
#
#   # Report messages for the inputs
#
#   cli::cli_alert("                       Input samples : {.field n = {nrow(x)}}")
#   cli::cli_alert("                     Output clusters : {.field k = {k}}.")
#   cli::cli_alert("                        Blacklist by : {.field {blacklist}}.")
#   cli::cli_alert("                  Maximum iterations : {.field {max_iterations}}.")
#
#   # Report messages for the reference
#   cli::cli_alert(
#     "Reference Catalogue Signatures (RCSs): {.field {rownames(reference_catalogue) %>% sig_col}} ({.field n = {nrow(reference_catalogue)}})"
#   )
#
#   # Report messages for the input
#   if (!is.null(input_catalogue)) {
#     cli::cli_alert(
#       "    Input Catalogue Signatures (ICSs): {.field {head(rownames(input_catalogue)) %>% sig_col}, ...} ({.field n = {nrow(input_catalogue)}})"
#     )
#     input_catalogue %>% dplyr::as_tibble() %>% print()
#   } else {
#     cli::cli_alert_warning("    Input Catalogue Signatures (ICSs): none (no supervision).")
#   }
#
#   cli::cli_h2("Iterative Basilica inference algorithm")
#
#   TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")
#
#   counter <- 1
#   black_list <- c()
#
#   # Iterative objects
#   BIC_trajectories = ICS_trajectories = DNS_trajectories =
#     blacklist_trajectories = NULL
#
#   # repeat untill convergence
#   repeat {
#     cli::cli_h2("Basilica step {.field {counter}}")
#
#     ICS_trajectories = append(ICS_trajectories, list(input_catalogue))
#     blacklist_trajectories = append(blacklist_trajectories, list(black_list))
#
#     n_ICSs = ifelse(is.null(input_catalogue), 0, nrow(input_catalogue))
#     cli::cli_h3("Bayesian NMF via SVI [{.field {steps}} steps, ICSs {.field k = {n_ICSs}}]")
#
#     TIME_it = as.POSIXct(Sys.time(), format = "%H:%M:%S")
#
#     ################### Bayesian NMF
#     k_aux = k
#     if(0 %in% k & n_ICSs == 0) {
#       k = k[k != 0]
#       cli::cli_alert_info("ICSs {.field k = {n_ICSs}}, removing {.field k = {0}} for this run, using {.field k = {k}} .")
#
#       if(length(k) == 0) cli::cli_abort("Cannot proceed: ICSs {.field k = {n_ICSs}} and {.field k = {0}}.")
#     }
#
#     obj = pyfit(
#       x = x,
#       py = py,
#       k_list = k,
#       lr = lr,
#       n_steps = steps,
#       groups = groups,
#       input_catalogue = input_catalogue,
#       # lambda_rate = lambda_rate,
#       # sigma = sigma,
#       hyperparameters = hyperparameters,
#       CUDA = CUDA,
#       compile = compile,
#       enforce_sparsity = enforce_sparsity,
#       store_parameters = store_parameters,
#       regularizer = regularizer,
#       reg_weight = reg_weight,
#       reg_bic = reg_bic,
#       stage = stage,
#       verbose = verbose
#     )
#
#     if (filtered_cat && !is.null(obj$denovo_signatures) && nrow(obj$denovo_signatures)>0)
#       obj$denovo_signatures = renormalize_denovo_thr(obj$denovo_signatures)
#
#     k = k_aux
#
#     TIME_it = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"),
#                        TIME_it, units = "mins") %>% round(2)
#
#     loss = obj$losses[steps] %>% round(2)
#     bic = obj$bic %>% round(2)
#     BIC_trajectories = c(BIC_trajectories, bic)
#
#     best_k = obj$exposure %>% ncol()
#     fix_k = input_catalogue %>% nrow()
#     if(fix_k %>% is.null) fix_k = 0
#
#     best_k = best_k - fix_k
#
#     cli::cli_alert_success(
#       paste0(
#         "NMF completed in {.field {TIME_it}} minutes [BIC {.field {bic}}, Loss {.field {loss}}, best DNSs {.field k = {best_k}} given ICSs {.field k = {n_ICSs}}]. "
#       )
#     )
#
#     ################### Filter out small-exposure signatures
#     cli::cli_h3("Checking ICSs exposure, \u03b1 > \u03A6 ({.field \u03A6 = {phi}}).")
#
#     if (is.null(input_catalogue) || nrow(input_catalogue) == 0)
#       cli::cli_alert_danger("No ICSs in this step were used.")
#
#
#     if(!is.null(blacklist)) {
#
#       # drop non-significant fixed signatures --------------------------------
#       if(!is.null(blacklist) && blacklist == "TMB")
#
#         if(!is.null(blacklist)) {
#           if(blacklist == "TMB")
#             a <- filter.fixed(
#               M = x,
#               alpha = obj$exposure,
#               beta_fixed = input_catalogue,
#               phi = phi)
#
#           if(!is.null(blacklist) && blacklist == "freq")
#             a = filter.fixed_minfreq(
#               alpha = obj$exposure,
#               beta_fixed = input_catalogue,
#               phi = phi)
#
#         } else {
#           a = filter.fixed_nofilter( # fake function
#             alpha = obj$exposure,
#             beta_fixed = input_catalogue)
#         }
#
#     } else {
#       a = filter.fixed_nofilter( # fake function
#         alpha = obj$exposure,
#         beta_fixed = input_catalogue)
#     }
#
#
#     remained_fixed <- a$remained_fixed # data.frame / NULL
#
#     # Here we might have dropped something
#     if (!is.null(a$dropped_fixed)) {
#       n_dropped = a$dropped_fixed %>% nrow()
#       w_dropped = a$dropped_fixed %>% rownames
#
#       cli::cli_alert_warning(
#         paste0(
#           "{.field {n_dropped}} ICSs will be dropped ({.field {w_dropped %>% sig_col}}), and blacklisted."
#         )
#       )
#
#       black_list <- union(black_list, w_dropped)  # character vector
#
#       cli::cli_alert_warning(paste0(
#         "The updated blacklist is: {.field {black_list %>% crayon::red()}}."
#       ))
#     } else {
#
#       if (!is.null(input_catalogue)) {
#         n_cat = input_catalogue %>% nrow()
#
#         # Here we did NOT drop, and had some input signatures - we report it
#         if (n_cat > 0)
#           cli::cli_alert_success("No ICSs will be dropped.")
#       }
#     }
#
#     # TEST
#     # cat('class(a$dropped_fixed):', class(a$dropped_fixed), '\n')
#     # cat('a$dropped_fixed:')
#     # print(a$dropped_fixed)
#     # cat('class(black_list):', class(black_list), '\n')
#     # cat('black_list:', black_list, '\n')
#     # TEST
#
#     n_denovo = obj$denovo_signatures %>% nrow
#     n_catalogue = reference_catalogue %>% nrow
#
#     cli::cli_h3(
#       "Comparing {.field {n_denovo}} DNSs to {.field {n_catalogue}} RCSs, imposing cosine similarlity above {.field \u0394 = {delta}}."
#     )
#
#     # detect denovo signatures which are similar to reference signatures -------
#     if (cosine_by_subs)
#       substitutions = get_contexts(obj) %>% dplyr::pull(subs) %>% unique() else
#         substitutions = NULL
#
#     obj$exposure = filter.denovo.phi(exposures = obj$exposure,
#                                      denovo = obj$denovo_signatures %>% rownames(),
#                                      phi = phi) # filter denovo if low exposure in all patients
#
#     # check if reference are linear comb of denovo sigs co-occurring (exp>thr)
#     b_denovo <- filter.denovo.QP(
#       reference = rbind(input_catalogue, reference_catalogue),
#       # beta_fixed = rbind(input_catalogue, reference_catalogue),
#       beta_denovo = obj$denovo_signatures,
#       thr_exposure = phi,
#       exposures = obj$exposure,
#       black_list = black_list,
#       delta = delta,
#       filt_pi = filt_pi,
#       substitutions = substitutions)
#
#     if (!is.null(b_denovo$new_fixed))
#       print(b_denovo)
#
#     denovo_filt = b_denovo$reduced_denovo
#
#
#     if (!is.null(rbind(input_catalogue, b_denovo$new_fixed))) {
#       ref = setdiff(reference_catalogue,
#                     rbind(input_catalogue, b_denovo$new_fixed))
#     } else{
#       ref = reference_catalogue
#     }
#
#     if (!is.null(denovo_filt) && nrow(denovo_filt) > 0)
#       # check if denovo are linear comb of reference sigs
#       b_reference <- filter.denovo.QP(
#         reference = ref,
#         # beta_fixed = rbind(input_catalogue, b_denovo$new_fixed),
#         beta_denovo = denovo_filt,
#         black_list = black_list,
#         delta = delta,
#         filt_pi = filt_pi,
#         substitutions = substitutions) else b_reference = list("new_fixed"=NULL, "reduced_denovo"=NULL)
#
#     new_fixed <- rbind(b_denovo$new_fixed, b_reference$new_fixed) %>% unique()
#
#     # if (!is.null(b_denovo$new_fixed))
#
#     if (!is.null(new_fixed) && nrow(new_fixed) > 0) {
#       n_denovo_denovo = b_reference$reduced_denovo %>% nrow()
#       # cli::cli_alert_warning(
#       #   "{.field {nrow(new_fixed)}}/{.field {n_denovo}} DNSs were found in the reference catalogue, and will become part of the ICSs!"
#       # )
#     } else {
#       # cli::cli_alert_success("No DNSs were found in the reference catalogue!")
#     }
#
#     reduced_denovo <- b_reference$reduced_denovo  # data.frame / NULL (remaining denovo signatures)
#
#     DNS_trajectories = append(DNS_trajectories, list(reduced_denovo))
#
#     #TEST---------------------------------------------------------
#     # cat("        fixed          :", rownames(input_catalogue), '\n')
#     # cat("        remained fixed :", rownames(remained_fixed), '\n')
#     # cat("        black list     :", black_list, "\n\n")
#     # cat("        denovo         :", rownames(obj$denovo_signatures), '\n')
#     # cat("        new fixed      :", rownames(new_fixed), '\n')
#     # cat("        reduced denovo :", rownames(reduced_denovo), '\n')
#     #TEST---------------------------------------------------------
#
#     if (is.null(input_catalogue)) {
#       col_names <- colnames(x)
#       input_catalogue = data.frame(matrix(nrow = 0, ncol = length(col_names)))
#       colnames(input_catalogue) = col_names
#     }
#
#     if (is.null(new_fixed)) {
#       col_names <- colnames(x)
#       new_fixed = data.frame(matrix(nrow = 0, ncol = length(col_names)))
#       colnames(new_fixed) = col_names
#     }
#
#     if (is.null(remained_fixed)) {
#       col_names <- colnames(x)
#       remained_fixed = data.frame(matrix(nrow = 0, ncol = length(col_names)))
#       colnames(remained_fixed) = col_names
#     }
#
#     # Convergency test:
#     # - we keep finding exactly the ICSs
#     # - there are no DNSs that seem to actually be part of RCSs
#     if (nrow(dplyr::setdiff(input_catalogue, remained_fixed)) == 0 &
#         nrow(new_fixed) == 0) {
#       cat('\n')
#       cli::cli_alert_success("Converged - ICSs is stable and all DNSs are genuine.")
#
#       break
#     }
#
#     if (nrow(remained_fixed) == 0 & nrow(new_fixed) == 0) {
#       input_catalogue <- NULL
#     } else {
#       input_catalogue <- rbind(remained_fixed, new_fixed)
#
#       cli::cli_alert_warning(
#         paste0(
#           "The updated ICSs contains now: {.field {input_catalogue %>% rownames %>% sig_col}}."
#         )
#       )
#     }
#
#     counter <- counter + 1
#     if (counter > max_iterations) {
#
#       cat('\n')
#       cli::cli_alert_danger("Converged forced after {.value {crayon::red(counter)}} iterations.")
#
#       break
#     }
#     # cat('    ------------------------------------\n')
#   }
#   # /end repeat untill convergence
#
#
#   if (nrow(input_catalogue) == 0) {
#     obj$catalogue_signatures <- NULL
#   } else {
#     obj$catalogue_signatures <- input_catalogue
#   }
#
#   # output ---> dtype: list
#   #-------------------------------------:
#   # exposure              --> data.frame
#   # denovo_signatures     --> data.frame
#   # bic                   --> numeric
#   # losses                --> numeric
#   # catalogue_signatures  --> data.frame
#
#   TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"),
#                   TIME, units = "mins")
#   TIME = TIME %>% round(2)
#
#   cli::cli_h3("Basilica completed in {.field {TIME}} minutes and {.field {counter}} iterations.")
#
#   if (length(BIC_trajectories) > 1) {
#     cli::cli_h3("Fit statistics")
#
#     # Print the BIC trajectory in +/- %
#     BIC_trajectories_delta = sapply(1:(length(BIC_trajectories) - 1),
#                                     function(i) {
#                                       d = BIC_trajectories[i + 1] - BIC_trajectories[i]
#                                       d_abs = abs(d)
#                                       d_pro = (d_abs / BIC_trajectories[i])
#                                       if (d > 0)
#                                         d_pro * 100
#                                       else-1 * d_pro  * 100
#                                     })
#     BIC_trajectories_delta = BIC_trajectories_delta %>% round(2)
#
#     BIC_trajectories_delta = sapply(BIC_trajectories_delta,
#                                     function(x) {
#                                       if (x > 0)
#                                         paste0(x, '%') %>% crayon::green()
#                                       else
#                                         paste0(x, '%') %>% crayon::red()
#                                     })
#
#     paste0(
#       "> BIC ",
#       BIC_trajectories[1],
#       " : ",
#       paste(BIC_trajectories_delta, collapse = " \u2192 ")
#     ) %>% cat()
#   }
#
#
#   # Store in the output obj all relevant information
#   fit$cohort = cohort
#
#   fit$n_samples = x %>% nrow()
#
#   fit$time = TIME
#
#   #### CHECK ########
#   # cat_expo = obj$denovo_signatures %>% rownames()
#   # cat_expo = setdiff(obj$exposure %>% colnames, cat_expo)
#   # fit$n_catalogue = obj$exposure[, cat_expo] %>% ncol()
#
#   fit$n_denovo = ifelse(
#     obj$denovo_signatures %>% is.null,
#     0,
#     obj$denovo_signatures %>% nrow
#   )
#
#   if (is.matrix(obj$alpha)) {
#     obj$alpha = obj$alpha %>% as.data.frame()
#     colnames(obj$alpha) = c(obj$catalogue_signatures %>% rownames(),
#                             obj$denovo_signatures %>% rownames())
#   }
#
#
#   fit$input = list(
#     counts = x,
#     reference_catalogue = reference_catalogue,
#     input_catalogue = ICS_trajectories[[1]]
#   )
#
#   fit$params = list(
#     k = k,
#     lr = lr,
#     steps = steps,
#     phi = phi,
#     delta = delta,
#     groups = groups,
#     lambda_rate = lambda_rate,
#     sigma = sigma
#   )
#
#   fit$iterations = list(
#     BIC = BIC_trajectories,
#     ICS = ICS_trajectories,
#     DNS = DNS_trajectories,
#     blacklist = blacklist_trajectories
#   )
#
#   fit$fit = obj
#
#   fit$k_list = k
#
#
#   return(fit)
# }
#
#
# sanitize_inputs = function(x,
#                            reference_catalogue,
#                            k,
#                            lr,
#                            steps,
#                            phi,
#                            delta,
#                            groups = NULL,
#                            input_catalogue = NULL,
#                            lambda_rate = NULL,
#                            sigma = FALSE,
#                            blacklist)
#
# {
#   # Input counts
#   if (!is.data.frame(x))
#     cli::cli_abort("The count matrix should be a dataframe!")
#   if (nrow(x) == 0)
#     cli::cli_abort("The count matrix has no rows!")
#
#   if (rownames(x) %>% is.null())
#   {
#     cm = paste0("Sample_", 1:nrow(x))
#     rownames(x) = cm
#
#     cli::cli_alert_warning("Sample names were missing, will use {.field {head(cm)}, ...}")
#   }
#
#   # Input reference_catalogue
#   if (!is.data.frame(reference_catalogue))
#     cli::cli_abort("The reference catalogue should be a dataframe!")
#   if (nrow(reference_catalogue) == 0)
#     cli::cli_abort("The reference catalogue has no rows!")
#
#   if (rownames(reference_catalogue) %>% is.null)
#     cli::cli_abort("The reference catalogue has no rownames (signature names)!")
#
#   # What is there
#   if (!all(colnames(x) %in% colnames(reference_catalogue)))
#     cli::cli_abort("Some columns in the input counts miss from the reference catalogue")
#
#   # ordering
#   if (!all(colnames(x) == colnames(reference_catalogue)))
#   {
#     cli::cli_alert_warning("Reference signature and catalgous have different column orders, will re-order")
#
#     reference_catalogue = reference_catalogue[names(x)]
#   }
#
#   # Input catalogue
#   if (!is.null(input_catalogue))
#   {
#     if (!is.data.frame(input_catalogue))
#       cli::cli_abort("The input catalogue should be a dataframe!")
#     if (nrow(input_catalogue) == 0)
#       cli::cli_abort("The input catalogue has no rows!")
#
#     if (rownames(input_catalogue) %>% is.null)
#       cli::cli_abort("The input catalogue has no rownames (signature names)!")
#
#     # Ref inluded in input
#     if (!all(rownames(input_catalogue) %in% rownames(reference_catalogue)))
#       cli::cli_abort("The input catalogue has signatures that are not in there reference!")
#
#     # What is there
#     if (!all(colnames(x) %in% colnames(input_catalogue)))
#       cli::cli_abort("Some columns in the input counts miss from the input catalogue")
#
#     # ordering
#     if (!all(colnames(x) == colnames(input_catalogue)))
#     {
#       cli::cli_alert_warning("Reference catalgous and input counts have different column orders, will re-order")
#
#       input_catalogue = input_catalogue[names(x)]
#     }
#
#     input_catalogue = input_catalogue %>% as.data.frame()
#   }
#
#   # Numerics
#   if (!(is.numeric(k) |
#         (!is.list(k) & all(sapply(k, is.numeric)))))
#     cli::cli_abort("k must be a list of integers, or a single integer.")
#
#   if (!(is.numeric(lr) &
#         lr > 0))
#     cli::cli_abort("Invalid learning rate.")
#   if (!(is.numeric(steps) &
#         steps > 0))
#     cli::cli_abort("Invalid number of steps.")
#
#   # TODO - complete
#   #   phi,
#   #   delta,
#   #   groups = NULL,
#   #   input_catalogue = NULL,
#   #   lambda_rate = NULL,
#   #   sigma = FALSE
#   # blacklist
#   # max_iterations (to be added yet...)
#
#
#   return(
#     list(
#       x = x %>% as.data.frame(),
#       reference_catalogue = reference_catalogue %>% as.data.frame(),
#       input_catalogue = input_catalogue
#     )
#   )
# }
#
#
#
#
#
#
# # get_reconstructed_data = function(x){
# #
# #   signatures = NULL
# #   if("catalogue_signatures" %in% names(x$fit)){ signatures = rbind(signatures,x$fit$catalogue_signatures) }
# #   if("denovo_signatures" %in% names(x$fit)){ signatures = rbind(signatures,x$fit$denovo_signatures)}
# #
# #   signatures = signatures[colnames(x$fit$exposure),]
# #
# #   as.matrix(x$fit$exposure[rownames(x$input$counts),]*rowSums(x$input$counts)) %*%  as.matrix(signatures) %>%
# #     as.data.frame()
# #
# # }
#
# get_similarity_scores = function(reconstr_data,real_data){
#
#   reconstr_data = reconstr_data[rownames(real_data),]
#   cosines = lapply(1:nrow(reconstr_data),function(j){
#
#     tibble(id = rownames(real_data)[j], cos =  cosine.vector(reconstr_data[j,],real_data[j,]))
#
#   }) %>% bind_rows()
#
#   cosines
#
# }
#
#
#
# plot_gof = function(scores){
#
#   if(!"group" %in% (scores %>% colnames())){  scores = scores %>% mutate(group = "1")}
#
#   ggplot(scores, aes(x = group,y =  cos,fill = group)) + geom_boxplot() + ylim(0,1) +
#     CNAqc:::my_ggplot_theme() + labs(title = "GOF of data counts", y = "GOF")  + ggsci::scale_fill_lancet()
#
# }
#
#
# plot_cohort_prevalence <- function(x) {
#
#   alpha <- get_exposure(x, long = TRUE) %>%
#     dplyr::group_by(Signature) %>%
#     dplyr::summarise(Exposure = sum(Exposure)) %>%
#     dplyr::mutate(Exposure = Exposure/sum(Exposure)) %>%
#     dplyr::arrange(dplyr::desc(Exposure))
#
#   alpha$Signature = factor(alpha$Signature,
#                            levels = alpha$Signature)
#
#
#   plt = ggplot2::ggplot(
#     data = alpha,
#     ggplot2::aes(x="Overall", y=Exposure, fill=Signature)
#   ) +
#     ggplot2::geom_bar(stat = "identity") +
#     my_ggplot_theme() +
#     ggplot2::scale_y_continuous(labels=scales::percent) +
#     ggplot2::scale_fill_manual(values = get_signature_colors(x)) +
#     ggplot2::theme(
#       axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
#     ) +
#     ggplot2::labs(
#       title = paste0(x$cohort, ' (n = ', x$n_samples, ')'),
#       y = "Overall exposure",
#       x = ''
#     ) +
#     ggplot2::guides(
#       fill = ggplot2::guide_legend(
#         nrow = ifelse(x$n_catalogue + x$n_denovo > 8, 2, 1))
#     )
#
#   alpha <- get_exposure(x, long = TRUE)
#
#   plt2 = ggplot(alpha,
#                 aes(x = Signature, y = Exposure, color = Signature))+
#     geom_jitter(size = .5, alpha = .6) +
#     geom_boxplot(color = 'black', width = .3) +
#     # scale_fill_manual(values = get_signature_colors(x)) +
#     scale_color_manual(values = get_signature_colors(x)) +
#     my_ggplot_theme() +
#     geom_hline(
#       yintercept = x$params$phi,
#       color = 'indianred3',
#       linetype = 'dashed'
#     ) +
#     labs(title = "Per sample exposure") +
#     ggplot2::theme(
#       axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
#     )
#
#   alpha <- get_exposure(x, long = TRUE)
#
#   seq_cuts = seq(0, 1, by = 0.01)
#
#   regression_data = lapply(seq_cuts,
#                            function(s)
#                            {
#                              alpha %>%
#                                filter(Exposure >= s) %>%
#                                group_by(Signature) %>%
#                                summarise(n = n(), cut = s)
#                            }
#   ) %>%
#     Reduce(f = bind_rows)
#
#   pl3 = regression_data %>%
#     ggplot() +
#     geom_line(
#       aes(x = cut, y = n, color = Signature)
#     ) +
#     geom_point(
#       data = regression_data %>% group_by(Signature) %>% filter(n == min(n)) %>% filter(row_number() == n()),
#       aes(x = cut, y = n, color = Signature)
#     ) +
#     scale_color_manual(values = get_signature_colors(x)) +
#     my_ggplot_theme()
#
#   pl4 = pl3 +
#     scale_y_log10()
#
#   plt = ggpubr::ggarrange(
#     plt,
#     plt2,
#     ncol = 2,
#     widths = c(1, 2),
#     common.legend = TRUE,
#     legend = 'bottom'
#   )
#
#   plt34 = ggpubr::ggarrange(pl3, pl4, ncol = 2)
#
#   plt = ggpubr::ggarrange(plt, plt34, ncol = 1,
#                           common.legend = TRUE,
#                           legend = 'bottom')
#
#
#   return(plt)
# }
#
#
#
#
# plot_compare_fits <- function(x, y, similarity_cutoff = 0.4) {
#   # Similarity to the reference
#   cosine_matrix <- cosine.matrix(x %>% get_signatures(),
#                                  y %>% get_signatures())
#
#   # Nice colors
#   color_gradient = (RColorBrewer::brewer.pal(10, 'Spectral')) %>% rev
#   color_gradient[1:5] = color_gradient[1:5] %>% ggplot2::alpha(0.7)
#   color_breaks = seq(0, 1, 0.1)
#
#   # Numbers where worth
#   numbers = cosine_matrix %>% round(2)
#   numbers[numbers < similarity_cutoff] = ''
#
#   # The world is a better place now that I can pheatmap -> ggplot
#   ggp = pheatmap::pheatmap(
#     mat = cosine_matrix,
#     color = color_gradient,
#     breaks = color_breaks,
#     border_color = 'white',
#     cellwidth = 25,
#     cellheight = 15,
#     display_numbers = numbers
#   ) %>% ggplotify::as.ggplot()
#
#   # Closest matches: x to y and viceverse
#   x_y = cosine_matrix %>% as.matrix() %>% reshape2::melt()
#
#   x_y_best = x_y %>% group_by(Var1) %>% filter(value == max(value))
#   y_x_best = x_y %>% group_by(Var2) %>% filter(value == max(value))
#
#   extra_plots = NULL
#
#   x_y_best = easypar::run(
#     FUN = function(i){
#       comparative_plot_sbs(
#         x %>% get_signatures(long = TRUE) %>% dplyr::filter(Signature == x_y_best$Var1[i]),
#         y %>% get_signatures(long = TRUE) %>% dplyr::filter(Signature == x_y_best$Var2[i])
#       ) +
#         labs(title = paste0("x (", x_y_best$Var1[i], ')',
#                             "vs y (", x_y_best$Var2[i], ') = ',
#                             x_y_best$value[i] %>% round(3)))
#     },
#     PARAMS = lapply(1:nrow(x_y_best), list),
#     parallel = FALSE
#   )
#
#   y_x_best = easypar::run(
#     FUN = function(i){
#       comparative_plot_sbs(
#         x %>% get_signatures(long = TRUE) %>% dplyr::filter(Signature == y_x_best$Var1[i]),
#         y %>% get_signatures(long = TRUE) %>% dplyr::filter(Signature == y_x_best$Var2[i])
#       ) +
#         labs(title = paste0("x (", y_x_best$Var1[i], ')',
#                             "vs y (", y_x_best$Var2[i], ') = ',
#                             y_x_best$value[i] %>% round(3)))
#     },
#     PARAMS = lapply(1:nrow(y_x_best), list),
#     parallel = FALSE
#   )
#
#   # np = extra_plots %>% length()
#   # nc = np %>% sqrt() %>% ceiling()
#   # nr = nc
#   # if ((nr - 1) * nc >= np)
#   #   nr = nr - 1
#   #
#   # extra_plots = ggpubr::ggarrange(
#   #   plotlist = extra_plots,
#   #   ncol = nc,
#   #   nrow = nr,
#   #   common.legend = TRUE,
#   #   legend = 'bottom'
#   # )
#
#   # ggp = ggpubr::ggarrange(ggp, extra_plots, ncol = 2, widths = c(1, nc))
#
#   return(list(cosine = ggp, x_y = x_y_best, y_x = y_x_best))
# }
#
#
# comparative_plot_sbs = function(x, y)
# {
#   x = x %>% mutate(Object = 'x')
#   y = y %>% dplyr::mutate(Value = -1 * Value) %>% mutate(Object = 'y')
#
#   sigs = bind_rows(x, y)
#
#   # Remove parenthesis
#   sigs$substitution = stringr::str_extract_all(sigs$Feature, "\\[[^\\]\\[]*]") %>% unlist()
#   sigs$substitution = gsub(
#     pattern = '\\[',
#     replacement = '',
#     x = sigs$substitution
#   )
#   sigs$substitution = gsub(
#     pattern = '\\]',
#     replacement = '',
#     x = sigs$substitution
#   )
#
#   sigs = sigs %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(context = paste0(substr(Feature, 1, 1), '_', substr(Feature, 7, 7)), )
#
#   max_range = sigs$Value %>% abs %>% max
#   brange = seq(-max_range, max_range, max_range / 5) %>% round(3)
#
#   ggplot2::ggplot(sigs) +
#     ggplot2::geom_bar(
#       ggplot2::aes(x = context, y = Value, fill = Object),
#       stat = "identity",
#       position = "identity"
#     ) +
#     facet_wrap(~ substitution, nrow = 1) +
#     my_ggplot_theme() +
#     theme(
#       axis.text.x = ggplot2::element_blank(),
#       axis.text.y = ggplot2::element_blank(),
#       axis.title.x = ggplot2::element_blank(),
#       axis.title.y = ggplot2::element_blank()
#     ) +
#     ylim(-max_range, max_range)
# }
#
#
#
#
# plot_simulation_report <- function(cohort) {
#
#   plotlist <- lapply(1:nrow(cohort),FUN =  function(i) plot_simulation_report_aux(cohort[i,]))
#
#   return(plotlist)
#
# }
#
#
# plot_simulation_report_aux <- function(x) {
#
#   library(patchwork)
#
#   fit <- x$fit[[1]]
#
#   plot_signatures_fit <- plot_signatures(fit)
#   sbs_real_fixed <- x$exp_fixed[[1]] %>% as.matrix() %>% reshape2::melt()
#   sbs_real_fixed$Type = "Catalogue"
#   sbs_real_denovo <- x$exp_denovo[[1]] %>% as.matrix() %>% reshape2::melt()
#   sbs_real_denovo$Type = "De novo"
#   sbs_real <- rbind(sbs_real_fixed, sbs_real_denovo)
#   colnames(sbs_real) <- c("Signature", "Feature", "Value", "Type")
#   plot_signatures_simulation <- plot_signatures_sbs(fit, sbs_real) + scale_fill_discrete() +
#     ggplot2::geom_bar(
#       ggplot2::aes(x = crazy_map, y = Value, fill = Signature, color = Type),
#       stat="identity",
#       position="identity") +
#     scale_color_manual("Type",values = c("grey10", "grey60", "white"), breaks = c("Catalogue", "De novo", "")) +
#     ggtitle(glue::glue("Real signatures ( n = {fit$n_samples})"))
#
#   plot_exposure_fit <- plot_exposure(fit)
#   exposure_real  <- x$exp_exposure[[1]] %>% as.matrix() %>% reshape2::melt() %>%
#     mutate(Type = if_else(Var2  %in% unique(sbs_real_denovo$Var2) ,"De novo", "Catalogue"))
#   colnames(exposure_real) <- c("Sample", "Signature", "Exposure", "Type")
#
#   plot_exposure_real <- ggplot2::ggplot(
#     data = exposure_real,
#     ggplot2::aes(x=Sample, y=Exposure, fill=Signature)
#   ) +
#     ggplot2::geom_bar(stat = "identity") +
#     my_ggplot_theme() +
#     ggplot2::scale_y_continuous(labels=scales::percent) +
#     ggplot2::theme(
#       axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
#     ) +
#     ggplot2::labs(
#       title = paste0(fit$cohort, ' (n = ', fit$n_samples, ')')    )
#
#   #fit_similarity <- plot_similarity_reference(fit)
#
#   stats <- evaluate.cohort(x)
#   ps_MAE <- simple_barplot(data.frame(value = stats$mae), "MAE")
#   ps_MSE <- simple_barplot(data.frame(value = stats$mse), "MSE")
#   ps_FA <- simple_barplot(data.frame(value = stats$fixed_acc), "Catalogue accuracy")
#   ps_DA <- simple_barplot(data.frame(value = stats$denovo_ratio), "De novo accuracy")
#   if(!is.null(stats$denovo_match[[1]])) {
#     ps_DS <- ggplot(stats$denovo_match[[1]] %>% as.data.frame(), aes(x = match, y = similarity)) +
#       geom_col(color = "grey30") + my_ggplot_theme() +
#       ylab("") + xlab("") + ggtitle("De novo similarity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   } else {
#     ps_DS <- ggplot()
#   }
#
#   ps_plot <- ps_MAE | ps_MSE | ps_FA | ps_DA | ps_DS
#
#   final_plot <- ggpubr::ggarrange(plot_signatures_fit ,
#                                   plot_signatures_simulation,
#                                   plot_exposure_fit,
#                                   plot_exposure_real,
#                                   ps_plot, ncol = 1)
#
#   return(final_plot)
# }
#
#
# simple_barplot <- function(df,title) {
#   ggplot(df, aes(x = "", y = value)) + geom_col(color = "grey30") + my_ggplot_theme() + ylab("") + xlab("") + ggtitle (title)
# }
#
#
#
#
#
#
# # plot exposure--------------------------------------------------------QC:almost
#
# .plot.synthetic.alpha <- function(exp_alpha, inf_alpha) {
#
#   if (!identical(rownames(exp_alpha), rownames(inf_alpha))) {
#     stop("expected and inferred exposure have no same sample names")
#   }
#
#   #rownames(exp_alpha) <- rownames(inf_alpha)  # just to be consistent, should be fixed later
#   exp_alpha$sample <- rownames(exp_alpha)
#   exp_alpha$sample <- factor(exp_alpha$sample, levels = as.character(1:length(exp_alpha$sample)))
#   exp_alpha_long <- tidyr::gather(exp_alpha, key="signature", value="exposure", c(-sample))
#   exp_alpha_long$type <- rep('expected', each=nrow(exp_alpha_long))
#
#   #inf_alpha <- x$exposure[[1]]
#   inf_alpha$sample <- rownames(inf_alpha)
#   inf_alpha$sample <- factor(inf_alpha$sample, levels = as.character(1:length(inf_alpha$sample)))
#   inf_alpha_long <- tidyr::gather(inf_alpha, key="signature", value="exposure", c(-sample))
#   inf_alpha_long$type <- rep('inferred', each=nrow(inf_alpha_long))
#
#   alpha <- rbind(exp_alpha_long, inf_alpha_long)
#
#   # TEST -----------------------------------------------------------------------
#   cols <- unique(alpha$signature)
#   exp_denovo <- cols[stringr::str_detect(cols, '_D', negate = FALSE)]
#   aa <- setdiff(cols, exp_denovo)
#   fixed <- aa[stringr::str_detect(aa, 'SBS', negate = FALSE)]
#   inf_denovo <- setdiff(aa, fixed)
#   new_cols <- c(fixed, exp_denovo, inf_denovo)
#   alpha$signature <- factor(alpha$signature, levels = new_cols)
#   # TEST -----------------------------------------------------------------------
#
#   plt <- ggplot(data = alpha, aes(x=sample, y=exposure, fill=signature)) +
#     geom_bar(stat = "identity") +
#     facet_grid(type ~ .) +
#     theme_minimal() +
#     ggtitle("Signatures exposure (Expected vs. Inferred)") +
#     scale_fill_brewer(palette='Spectral')
#   #ggsci::scale_fill_lancet()
#   #scale_y_continuous(labels=scales::percent)
#
#   return(plt)
# }
#
#
# # plot signatures cosine matrix ------------------------------------------------
#
# .plot.synthetic.beta.similarity <- function(exp_denovo, inf_denovo) {
#
#   if (is.null(exp_denovo) | is.null(inf_denovo)) {
#     cplot <- ggplot() +
#       theme_void() +
#       geom_text(aes(0,0,label='N/A')) +
#       xlab(NULL) #optional, but safer in case another theme is applied later
#   } else {
#     # expected vs inferred signatures cosine similarity matrix
#     cos <- basilica:::cosine.matrix(inf_denovo, exp_denovo)
#     cos1 <- tibble::rownames_to_column(cos, var = 'inferred_denovo')
#     cos_long <- tidyr::gather(cos1, key="expected_denovo", value="cosine_similarity", c(-inferred_denovo))
#     # plot data
#     cplot <- ggplot(data = cos_long, aes(x=expected_denovo, y=inferred_denovo, fill=cosine_similarity)) +
#       geom_tile() +
#       geom_text(aes(label = round(cosine_similarity, 3))) +
#       #scale_fill_gradient(low = "orange", high = "green") +
#       scale_fill_gradient2(low = "orange",
#                            mid = "white",
#                            high = "purple") +
#       ggtitle("Cosine similarity matrix (expected vs. inferred)") +
#       xlab("Expected") +
#       ylab("Inferred")
#   }
#
#   return(cplot)
# }
#
#
#
#
#
# filter.fixed = function(M,
#                         alpha,
#                         beta_fixed = NULL,
#                         phi = 0.05) {
#   if (!is.data.frame(M)) {
#     warning("invalid count matrix (M) !")
#   }
#   if (!is.data.frame(alpha)) {
#     warning("invalid exposure matrix (alpha) !")
#   }
#   if (!is.numeric(phi)) {
#     warning("invalid parameter phi !")
#   }
#
#   if (is.null(beta_fixed)) {
#     remained_fixed = NULL
#     dropped_fixed = NULL
#   } else if (is.data.frame(beta_fixed)) {
#     theta = matrix(rowSums(M), nrow = 1)
#     alpha0 = theta %*% as.matrix(alpha)
#     contribution = colSums(alpha0) / sum(alpha0)
#
#     atb = alpha0 %>%
#       as_tibble() %>%
#       reshape2::melt(id = NULL) %>%
#       dplyr::rename(Signature = variable, TMB = value)
#
#     ctb = contribution %>%
#       as.list() %>%
#       as_tibble %>%
#       reshape2::melt(id = NULL) %>%
#       dplyr::rename(Signature = variable, proportion = value)
#
#     predicate_collapsed = dplyr::full_join(atb, ctb, by = 'Signature') %>%
#       dplyr::arrange(dplyr::desc(proportion)) %>%
#       mutate(Signature = ifelse(proportion < phi, crayon::red(Signature), Signature))
#
#     predicate_collapsed$TMB = paste0("TMB = ",
#                                      predicate_collapsed$TMB %>% round(0))
#
#     predicate_collapsed$proportion = paste0("\u03c0 = ",
#                                             predicate_collapsed$proportion %>% round(3))
#
#     predicate_collapsed$Signature = sprintf("%20s", predicate_collapsed$Signature)
#     predicate_collapsed$TMB = sprintf("%20s", predicate_collapsed$TMB)
#     predicate_collapsed$proportion = sprintf("%20s", predicate_collapsed$proportion)
#
#     predicate_collapsed = apply(predicate_collapsed, 1, function(x)
#       paste(x, collapse = ' '))
#
#     cli::boxx(
#       predicate_collapsed,
#       header = "TMB filter",
#       float = 'center',
#       footer = paste0("\u03c0 > ", phi)
#     ) %>% cat()
#     cat('\n')
#
#     dropped = which(contribution < phi)
#
#     if (sum(dropped) == 0) {
#       remained_fixed = beta_fixed
#       dropped_fixed = NULL
#     } else {
#       remained_fixed = beta_fixed[-c(dropped),]
#       if (nrow(remained_fixed) == 0) {
#         remained_fixed = NULL
#       }
#
#       if (any(dropped > nrow(beta_fixed))) {
#         cli::boxx("AZAD this is a bug") %>% cat()
#         dropped = dropped[dropped < nrow(beta_fixed)]
#       }
#
#       dropped_fixed = beta_fixed[c(dropped),]
#     }
#   } else {
#     warning("invalid fixed signatures (beta_fixed) !")
#   }
#   return(list(remained_fixed = remained_fixed, dropped_fixed = dropped_fixed))
# }
#
#
#
#
# filter.fixed_minfreq = function(alpha, beta_fixed, phi = 0.15)
# {
#   if (is.null(beta_fixed))
#     return(list(remained_fixed = NULL, dropped_fixed = NULL))
#
#   if (!is.null(alpha)) {
#     alpha_cat = alpha[, rownames(beta_fixed)]
#     predicate = apply(alpha_cat, 2, function(x)
#       sum(x > phi))
#
#     stays = colnames(alpha_cat)[predicate > 0]
#     goes = colnames(alpha_cat)[predicate == 0]
#
#     if (length(stays) == 0)
#       remained_fixed = NULL
#     else
#       remained_fixed = beta_fixed[stays,]
#
#     if (length(goes) == 0)
#       dropped_fixed = NULL
#     else
#       dropped_fixed = beta_fixed[goes,]
#
#     predicate_collapsed = predicate %>% as_tibble()
#     predicate_collapsed$Signature = names(predicate)
#
#     predicate_collapsed = predicate_collapsed %>%
#       group_by(value) %>%
#       mutate(Signature = paste(Signature, collapse = ', ')) %>%
#       distinct() %>%
#       arrange(value) %>%
#       mutate(Signature = ifelse(value == 0, crayon::red(Signature), Signature))
#
#     predicate_collapsed = paste('n =',
#                                 predicate_collapsed$value,
#                                 '[',
#                                 predicate_collapsed$Signature,
#                                 ']')
#
#     cli::boxx(
#       predicate_collapsed,
#       header = "Frequency filter",
#       float = 'center',
#       footer = paste0("\u03C6 > ", phi, " in n samples")
#     )  %>% cat()
#
#     cat('\n')
#
#     return(list(remained_fixed = remained_fixed, dropped_fixed = dropped_fixed))
#   }
# }
#
#
#
#
# filter.denovo.phi = function(exposures, phi, denovo) {
#   exposures.denovo = exposures[,denovo] %>% as.data.frame()
#   colnames(exposures.denovo) = denovo
#   sbs_low_exp = (exposures.denovo < phi) %>% colSums()
#
#   rmv = sbs_low_exp[sbs_low_exp > nrow(exposures.denovo)] %>% names()
#   if (length(rmv) == 0) return(exposures)
#
#   return(
#     exposures %>% dplyr::select(-dplyr::all_of(rmv)) %>% renormalize_denovo_thr(thr=0)
#   )
# }
#
#
#
#
# cosine.vector <- function(a, b, substitutions = NULL) {
#
#   if (is.matrix(a) && nrow(a)>1) a = t(a)
#   if (is.matrix(b) && nrow(b)>1) b = t(b)
#
#   if (is.null(substitutions)) {
#     if (!identical(colnames(a), colnames(b))) {
#       a = a[names(b)]
#     }
#
#     numerator <- sum(a * b)
#     denominator <- sqrt(sum(a^2)) * sqrt(sum(b^2))
#     return(numerator / denominator)
#   }
#
#   cosine.tot = 0
#   keep_subs = length(substitutions)
#   for (ss in substitutions) {
#     keep_cols.tmp = grep(ss, colnames(b), value = T)
#
#     if (all(c(a[,keep_cols.tmp], b[,keep_cols.tmp])==0)) {
#       keep_subs = keep_subs - 1
#       next
#     }
#
#     num.tmp = sum(a[,keep_cols.tmp] * b[,keep_cols.tmp])
#     denomin.tmp = sqrt(sum(a[,keep_cols.tmp]^2)) * sqrt(sum(b[,keep_cols.tmp]^2))
#
#     if (num.tmp != 0 && denomin.tmp != 0)
#       cosine.tot = cosine.tot + (num.tmp / denomin.tmp)
#   }
#
#   return(cosine.tot / keep_subs)
# }
#
#
#
#
#
# # Checks if a signature has at least one patient where it's exposure exceeds phi
#
#
# filter.fixed_nofilter = function(alpha, beta_fixed)
# {
#   if (is.null(beta_fixed))
#     return(list(remained_fixed = NULL, dropped_fixed = NULL))
#
#   if (!is.null(beta_fixed))
#     return(list(remained_fixed = beta_fixed, dropped_fixed = NULL))
#
#   if (!is.null(alpha))
#     return(list(remained_fixed = beta_fixed[colnames(alpha)[colnames(alpha) %in% rownames(beta_fixed)],],dropped_fixed = NULL))
#
# }
#
#
# filter.denovo = function(reference_catalogue,
#                          beta_fixed,
#                          beta_denovo = NULL,
#                          black_list = NULL,
#                          delta = 0.9) {
#   if (!is.data.frame(reference_catalogue)) {
#     warning("Invalid reference catalogue!")
#   }
#   if (!is.numeric(delta)) {
#     warning("Invalid delta argument!")
#   }
#
#   if (is.data.frame(beta_fixed)) {
#     reference = dplyr::setdiff(reference_catalogue, beta_fixed)
#   } else if (is.null(beta_fixed)) {
#     reference = reference_catalogue
#   } else {
#     warning('invalid fixed signatures (beta_fixed) !')
#   }
#
#   if (is.null(beta_denovo)) {
#     return(list(new_fixed = NULL, reduced_denovo = NULL))
#   } else if (is.data.frame(beta_denovo)) {
#     match_list = c()
#     cos_matrix = cosine.matrix(beta_denovo, reference)
#     while (TRUE) {
#       max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
#       if (cos_matrix[max] < delta) {
#         break
#       }
#       row_index = as.numeric(max)[1]
#       col_index = as.numeric(max)[2]
#       match_list[length(match_list) + 1] =
#         colnames(cos_matrix[col_index])
#
#       if (dim(cos_matrix)[1] == 1 | dim(cos_matrix)[2] == 1) {
#         cos_matrix = cos_matrix[-c(row_index),-c(col_index)]
#         break
#       } else {
#         cos_matrix = cos_matrix[-c(row_index),-c(col_index)]
#       }
#     }
#   } else {
#     warning("Invalid beta denovo!")
#   }
#
#   match_list = setdiff(match_list, black_list)
#   if (length(match_list) == 0) {
#     #col_names = colnames(reference_catalogue)
#     #df = data.frame(matrix(nrow=0, ncol = length(col_names)))
#     #colnames(df) = col_names
#     return(list(new_fixed = NULL, reduced_denovo = beta_denovo))
#   } else {
#     if (is.null(dim(cos_matrix))) {
#       return(list(new_fixed = reference[match_list,], reduced_denovo = NULL))
#     } else {
#       reduced_denovo = beta_denovo[rownames(cos_matrix),]
#       return(list(new_fixed = reference[match_list,], reduced_denovo = reduced_denovo))
#     }
#   }
# }
#
#
# adjust.denovo.fixed =
#   function(alpha,
#            fixed_signatures,
#            denovo_signatures,
#            limit = 0.9) {
#     if (is.null(fixed_signatures) | is.null(denovo_signatures)) {
#       return(list(exposure = alpha, denovo_signatures = denovo_signatures))
#     }
#
#     cos_matrix =
#       cosine.matrix(denovo_signatures, fixed_signatures)
#     while (TRUE) {
#       max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
#       if (cos_matrix[max] < limit) {
#         break
#       } else {
#         row_index = as.numeric(max)[1]
#         col_index = as.numeric(max)[2]
#         denovo_name = rownames(cos_matrix[col_index])[row_index]
#         fixed_name = colnames(cos_matrix[col_index])
#
#         # remove denovo signature
#         denovo_signatures =
#           denovo_signatures[!(rownames(denovo_signatures) %in% denovo_name),]
#
#         # adjust fixed signature exposure
#         alpha[fixed_name] =
#           alpha[, fixed_name] + alpha[, denovo_name]
#         # remove denovo signature exposure
#         alpha = alpha[,!names(alpha) %in% denovo_name]
#
#         if (dim(cos_matrix)[1] == 1 | dim(cos_matrix)[2] == 1) {
#           break
#         } else {
#           cos_matrix = cos_matrix[-c(row_index),-c(col_index)]
#         }
#       }
#     }
#     return(list(exposure = alpha, denovo_signatures = denovo_signatures))
#   }
#
#
# adjust.denovo.denovo =
#   function(alpha, denovo_signatures, limit = 0.9) {
#     if (is.null(denovo_signatures)) {
#       return(list(exposure = alpha, denovo_signatures = denovo_signatures))
#     } else if (nrow(denovo_signatures) == 1) {
#       return(list(exposure = alpha, denovo_signatures = denovo_signatures))
#     }
#     cos_matrix =
#       cosine.matrix(denovo_signatures, denovo_signatures)
#     for (i in 1:nrow(cos_matrix)) {
#       cos_matrix[i, i] = 0
#     }
#     while (TRUE) {
#       max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
#       if (cos_matrix[max][1] < limit) {
#         break
#       } else {
#         row_index = as.numeric(max)[1]
#         col_index = as.numeric(max)[2]
#         denovo_one = rownames(cos_matrix[col_index])[row_index]
#         denovo_two = colnames(cos_matrix[col_index])
#
#         denovo_signatures[paste(denovo_one, denovo_two, sep = ''),] =
#           denovo_signatures[denovo_one,] + denovo_signatures[denovo_two,]
#
#         denovo_signatures =
#           denovo_signatures[!(rownames(denovo_signatures) %in% c(denovo_one, denovo_two)),]
#
#         # adjust signature in exposure
#         alpha[paste(denovo_one, denovo_two, sep = '')] =
#           alpha[, denovo_one] + alpha[, denovo_two]
#         # remove signature from exposure
#         alpha =
#           alpha[,!names(alpha) %in% c(denovo_one, denovo_two)]
#
#         if (dim(cos_matrix)[1] == 2 | dim(cos_matrix)[2] == 2) {
#           break
#         } else {
#           cos_matrix =
#             cos_matrix[-c(row_index, col_index),-c(col_index, row_index)]
#         }
#       }
#     }
#     return(list(exposure = alpha, denovo_signatures = denovo_signatures))
#   }
#
#
# # filter based on linear projection with constraints
#
# filter.denovo.QP =
#   function(reference,
#            # beta_fixed,
#            beta_denovo = NULL,
#            black_list = NULL,
#            delta = 0.9,
#            filt_pi = 0.05,
#            thr_exposure = 0.05,
#            exposures = NULL,
#            substitutions = NULL) {
#     ## if denovo = TRUE -> check also if the denovo have exposure > thr in same samples
#
#     if (!is.data.frame(reference)) warning("Invalid reference catalogue!")
#     if (!is.numeric(delta)) warning("Invalid delta argument!")
#
#     # (Reference - Beta Fixed)
#     # if (is.data.frame(beta_fixed)) {
#     #   reference = dplyr::setdiff(reference_catalogue, beta_fixed)
#     # } else if (is.null(beta_fixed)) {
#     #   reference = reference_catalogue
#     # } else {
#     #   warning('invalid fixed signatures (beta_fixed) !')
#     # }
#
#     if (nrow(reference)==0) return(list(new_fixed = NULL, reduced_denovo = beta_denovo))
#
#     # BETA DENOVO
#     if (is.null(beta_denovo)) return(list(new_fixed = NULL, reduced_denovo = NULL))
#
#     if (is.data.frame(beta_denovo)) {
#       if (is.null(exposures)) {
#         a = beta_denovo
#         b = reference
#       } else {
#         a = reference
#         b = beta_denovo
#       }
#       ### names of catalogue signatures to include + names de novo to remove
#       res_optimization =
#         solve.quadratic.optimization(a,
#                                      b,
#                                      delta = delta,
#                                      filt_pi = filt_pi,
#                                      thr_exposure = thr_exposure,
#                                      exposures = exposures,
#                                      substitutions = substitutions)
#       match_list = res_optimization$catalogue_to_include
#     } else warning("Invalid beta denovo!")
#
#     match_list = setdiff(match_list, black_list)
#     if (length(match_list) == 0) {
#       # col_names = colnames(reference_catalogue)
#       # df = data.frame(matrix(nrow=0, ncol = length(col_names)))
#       # colnames(df) = col_names
#       return(list(new_fixed = NULL, reduced_denovo = beta_denovo))
#     } else {
#       if (is.null(res_optimization$denovo_to_include)) {
#         return(list(new_fixed = reference[match_list,], reduced_denovo = NULL))
#       } else {
#         reduced_denovo = beta_denovo[res_optimization$denovo_to_include,]
#         return(list(new_fixed = reference[match_list,], reduced_denovo = reduced_denovo))
#       }
#     }
#   }

