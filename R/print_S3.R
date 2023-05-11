#' Print for class \code{'basilica_obj'}.
#'
#' @param x An obj of class \code{'basilica_obj'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @exportS3Method print basilica_obj
#' @export print.basilica_obj

print.basilica_obj = function(x, ...)
{
  stopifnot(inherits(x, "basilica_obj"))

  cli::cli_rule(
    paste(
      crayon::bgYellow(crayon::black(paste0("[ Basilica ] ", x$cohort))),
      '{.field {x$n_samples}} samples with {.field {x$input$counts %>% sum}} total mutations.'
    )
  )

  # cli::cli_alert_info(paste0(" CNA segments: ", x$n_cna_clonal, " clonal, ", x$n_cna_sbclonal, " subclonal."))

  # cli::cli_alert_info(paste0("Mutation mapping (head): ", paste0(head(x$n_karyotype), ' (',
  # names(head(x$n_karyotype)), ')', collapse = '; ')))

  # cli::cli_alert_info(paste0("Mutation mapping (up to top 5): "))

  cli::cli_h3("Catalogue signatures (5 columns shown max)" %>% crayon::blue())

  ncol = 5
  writeLines(paste0("\t", capture.output(
    x$fit$catalogue_signatures[, 1:ncol] %>% print()
    )))


  cli::cli_h3("De novo signatures (5 columns shown max)"%>% crayon::yellow())

  writeLines(paste0("\t", capture.output(
    x$fit$denovo_signatures[, 1:ncol] %>% print()
  )))

  cat('\n')
  cli::cli_h3("Most-frequent signatures")
  bar_print_console(x)
  cat('\n')
}


bar_print_console = function(x, top = 10) {

  sorted = x$fit$exposure %>% colSums() %>% sort(decreasing = TRUE)
  sorted_p = sorted/sorted[1]

  max_cols = 30

  boxes = (sorted_p * max_cols) %>% round

  is_denovo = function(n){
    n %in% (x$fit$denovo_signatures %>% rownames())
  }

  sapply(boxes %>% seq, function(i){

    if(is_denovo(names(boxes)[i]))
      cat(
        "\t\t",
        sprintf( "%8s", names(boxes)[i]) %>% crayon::yellow(),
        paste(rep("\u25A0", boxes[i]), collapse = ''),
        '\n'
      )
    else
      cat(
        "\t\t",
        sprintf( "%8s", names(boxes)[i]) %>% crayon::blue(),
        paste(rep("\u25A0", boxes[i]), collapse = ''),
        '\n'
      )

  })



}




#' Plot for class \code{'basilica_obj'}.
#'
#' @description
#'
#' The default plot is the CNA segments in wide format
#'
#' @param x An obj of class \code{'basilica_obj'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' }
plot.basilica_obj = function(x, ...)
{
  stopifnot(inherits(x, "cnaqc"))


}
