#' Print for class \code{'bascule_obj'}.
#'
#' @param x An obj of class \code{'bascule_obj'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @exportS3Method print bascule_obj
#' @export print.bascule_obj

print.bascule_obj = function(x, ...)
{
  stopifnot(inherits(x, "bascule_obj"))

  cli::cli_rule(
    paste(
      crayon::bgYellow(crayon::black(paste0("[ Bascule ] ", x$cohort))),
      '{.field {x$n_samples}} samples with {.field {x$input$counts %>% sum}} total mutations.'
    )
  )

  ncol = 5

  lapply(get_types(x), function(t) {
    cli::cli_h3(paste0(t, " catalogue signatures (5 columns shown max)") %>% crayon::blue())

    writeLines(paste0("\t", capture.output(get_fixed_signatures(x, matrix=T)[[t]][, 1:ncol] %>% print())))


    if (!is.null(get_denovo_signatures(x))) {
      cli::cli_h3("De novo signatures (5 columns shown max)"%>% crayon::yellow())

      writeLines(paste0("\t", capture.output(
        get_denovo_signatures(x, matrix=T)[[t]][, 1:ncol] %>% print()
      )))

      # cat('\n')
      # cli::cli_h3("Most-frequent signatures")
      # bar_print_console(x, t=t)
      # cat('\n')
    }
  })
}

# bar_print_console = function(x, t, top = 10) {
#
#   sorted = get_exposure(x, matrix=T)[[t]] %>% colSums() %>% sort(decreasing = TRUE)
#   sorted_p = sorted/sorted[1]
#
#   max_cols = 30
#
#   boxes = (sorted_p * max_cols) %>% round
#
#   is_denovo = function(n) {
#     n %in% (get_denovo_signatures(x, matrix=T)[[t]] %>% rownames())
#   }
#
#   sapply(boxes %>% seq, function(i) {
#
#     if(is_denovo(names(boxes)[i]))
#       cat(
#         "\t\t",
#         sprintf( "%8s", names(boxes)[i]) %>% crayon::yellow(),
#         paste(rep("\u25A0", boxes[i]), collapse = ''),
#         '\n'
#       )
#     else
#       cat(
#         "\t\t",
#         sprintf( "%8s", names(boxes)[i]) %>% crayon::blue(),
#         paste(rep("\u25A0", boxes[i]), collapse = ''),
#         '\n'
#       )
#   })
# }


#' Plot for class \code{'bascule_obj'}.
#'
#' @description
#'
#' The default plot is the CNA segments in wide format
#'
#' @param x An obj of class \code{'bascule_obj'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @export plot.bascule_obj

plot.bascule_obj = function(x, ...)
{
  stopifnot(inherits(x, "cnaqc"))
}
