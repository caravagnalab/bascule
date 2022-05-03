

#' Title
#'
#' @param x basilica object
#'
#' @return plot exposure matrix
#' @export
#'
#' @examples
plot_exposure <- function(x) {
  alpha <- x$exposure
  alpha$Branch <- c(1:nrow(alpha))
  alpha_long <- gather(alpha,
                       key="Signature",
                       value="Exposure",
                       c(-Branch)
  )
  
  ggplot(data = alpha_long, aes(x=Branch, y=Exposure, fill=Signature)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    #ggtitle(paste(title)) +
    scale_y_continuous(labels=scales::percent)
}
