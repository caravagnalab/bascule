# # Sanitize input (check all patients are defined for every data type)
# sanitizer_long = function(x)
# {
#   shared_samples = lapply(
#     x,
#     function(x){
#       x$sample %>% unique()
#     }
#   ) %>%
#     Reduce(f = intersect)
#
#   test_in = sapply(
#     x, function(x) all(x$sample %in% shared_samples)
#   ) %>%
#     Reduce(f = all)
#
#   stopifnot(test_in)
# }
#
# # Internal, everything in long tibble format
# #
# # - input: list with data (sample | type | feature | value)
# # - nmf: list with Pyro fit, exposures and signatures tibbles (sample | type | signature | value)
# # - clustering: list with a tibble  (sample | cluster)
# # - stat: named vector with number of patients, number of data types, number of signatures per data type,
# #         number of clusters, fit scores, etc.
# # - alternative: alternative model fits
#
# # Get type of data used for signatures: e.g., "SBS", "DBS"
# get_types = function(x) { x$input %>% names() }
#
# # Get samples names
# get_samples = function(x) { x$input[[1]]$sample %>% unique() }
#
# # Get cluster labels ("C1", "C2")
# get_cluster_labels = function(x)
# {
#   if(is.null(x$clustering)) return(NULL)
#
#   x$clustering$cluster %>% unique()
# }
#
# # Get clustering assignments (tibble)
# get_cluster_assignments = function(x,
#                                    samples = get_samples(x),
#                                    clusters = get_cluster_labels(x)
# )
# {
#   if(is.null(x$clustering)) return(NULL)
#
#   x$clustering %>%
#     dplyr::filter(
#       sample %in% samples,
#       cluster %in% clusters
#     )
# }
#
# # Get input, can subset by types, samples and clustering assignment (if available), and
# # can transform the data in matrix format
# get_input = function(x,
#                      types = get_types(x),
#                      samples = get_samples(x),
#                      clusters = get_cluster_labels(x),
#                      matrix = FALSE
# )
# {
#   out = lapply(
#     types,
#     function(t){
#       w = x$input[[t]] %>%
#         dplyr::filter(sample %in% !!samples)
#
#       if(!is.null(clusters))
#       {
#         which_selection = x %>%
#           get_cluster_assignments(samples = samples, clusters = clusters) %>%
#           dplyr::pull(sample)
#
#         w = w %>% filter(sample %in% which_selection)
#       }
#
#       return(w)
#     })
#
#   if(matrix)
#   {
#     out = lapply(out,
#                  function(x){
#                    x %>%
#                      dplyr::select(-type) %>%
#                      tidyr::pivot_wider(
#                        names_from = feature,
#                        values_from = value
#                      )
#                  }
#     )
#   }
#
#   return(out)
# }
#
# # Get exposure, it can subset by types, samples and clusters. It can return
# # a list of matrices.
# get_exposure = function(x,
#                         types = get_types(x),
#                         samples = get_samples(x),
#                         clusters = get_cluster_labels(x),
#                         matrix = FALSE
# )
# {
#   out = lapply(
#     types,
#     function(t){
#       w = x$nmf[[t]]$exposure %>%
#         dplyr::filter(sample %in% !!samples)
#
#       if(!is.null(clusters))
#       {
#         which_selection = x %>%
#           get_cluster_assignments(samples = samples, clusters = clusters) %>%
#           dplyr::pull(sample)
#
#         w = w %>% filter(sample %in% which_selection)
#       }
#
#       return(w)
#     })
#
#   if(matrix)
#   {
#     out = lapply(out,
#                  function(x){
#                    x %>%
#                      dplyr::select(-type) %>%
#                      tidyr::pivot_wider(
#                        names_from = feature,
#                        values_from = value
#                      )
#                  }
#     )
#   }
# }
#
# # Get signatures, it can subset by types and return
# # a list of matrices.
# get_signatures = function(
#     x,
#     types = get_types(x),
#     matrix = FALSE
# ){
#   out = lapply(
#     types,
#     function(t){
#       w = x$nmf[[t]]$signatures %>%
#         dplyr::filter(sample %in% !!samples)
#
#       if(!is.null(clusters))
#       {
#         which_selection = x %>%
#           get_cluster_assignments(samples = samples, clusters = clusters) %>%
#           dplyr::pull(sample)
#
#         w = w %>% filter(sample %in% which_selection)
#       }
#
#       return(w)
#     })
#
#   if(matrix)
#   {
#     out = lapply(out,
#                  function(x){
#                    x %>%
#                      dplyr::select(-type) %>%
#                      tidyr::pivot_wider(
#                        names_from = feature,
#                        values_from = value
#                      )
#                  }
#     )
#   }
# }
#
# # NMF
# pyro_nmf = function(x, which){
#
#   fit = .... # pyro reticulate call
#
#   nmf = list(
#     pyro = fit,
#     exposure = ..., # extract in long format the exposure tibble
#     beta = ... # extract in long format the signatures tibble
#   )
#
#   x$nmf[[which]] = nmf # store the output
#
#   return(x)
# }
