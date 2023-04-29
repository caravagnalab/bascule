#' plot exposure matrix
#'
#' @description creates bar plot of relative exposure matrix, where x-axis are samples and y-axis are their relative contribution.
#' @param x basilica object
#' @param sample_name 
#' @param levels
#' @param cls
#'  @param flip_coord
#' @return plot
#' @export plot_exposure
#'
#' @examples

# plot_exposure <- function(x, labels = NULL,sort_by = NULL, thr=0.1){
# 
# 
#   alpha <- get_exposure(x, long = TRUE)
# 
#   # plt <- .plot_exposure(x = alpha)
#   if ("groups" %in% colnames(alpha))
#     samples_order = alpha %>% dplyr::arrange(groups) %>% dplyr::pull(Sample) %>% unique() else
#     samples_order = alpha$Sample %>% unique %>% gtools::mixedsort()
# 
#   if(!is.null(sort_by))
#   {
#     samples_order = alpha %>%
#       dplyr::filter(Signature == sort_by) %>%
#       dplyr::arrange(dplyr::desc(Exposure)) %>%
#       dplyr::pull(Sample)
#   }
# 
# 
#  alpha = alpha %>%
#     dplyr::mutate(Signature=ifelse(Exposure < thr, "Other", Signature))
# 
#   caption = paste0("Sorted by ", sort_by)
#   if(is.null(sort_by)) caption = "Sorted by sample"
#   
#   if(!is.null(labels)){
#     
#     alpha = alpha %>% dplyr::select(-groups) %>% full_join(labels, by = "Sample")
#     
# }
#   
# 
#   keep = alpha$Signature %>% unique()
# 
#   other_col = list("Other"="gainsboro") %>% unlist()
# 
#   plt = ggplot2::ggplot(
#     data = alpha,
#     ggplot2::aes(x=factor(Sample,levels = samples_order), y=Exposure, fill=Signature)
#   ) +
#     ggplot2::geom_bar(stat = "identity") +
#     my_ggplot_theme() +
#     ggplot2::scale_y_continuous(labels=scales::percent) +
#     ggplot2::scale_fill_manual(values = c(get_signature_colors(x),other_col),
#                                breaks = keep) +
#     ggplot2::theme(
#       axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
#     ) +
#     ggplot2::labs(
#       title = paste0(x$cohort, ' (n = ', x$n_samples, ')'),
#       caption = caption,
#       x = "Sample"
#     )
#   
#   # ggplot2::guides(
#   #   fill = ggplot2::guide_legend(
#   #     nrow = ifelse(x$n_catalogue + x$n_denovo > 8, 2, 1))
#   #   )
#   
#  if ("groups" %in% names(x$fit))
#     plt = plt + ggplot2::facet_grid(~groups)
# 
# 
#   return(plt)
# }


plot_exposure = function(x,sample_name = T,levels= NULL, cls = NULL, flip_coord = F){
  
  b = x$fit$exposure
  
  if(is.null(cls)){ cls = ggsci::pal_simpsons()(ncol(b))
  names(cls) = colnames(b)
  }
  
  if(is.null(levels)){ levels =   colnames(b) }
  
  p = ggplot(data = b %>% as.data.frame() %>% mutate(sample = rownames(b)) %>%
               reshape2::melt() %>% dplyr::rename(Signature = variable),
             aes(x = sample, y  = value,
                 fill = factor(Signature,levels = levels))) +
    geom_bar(stat = "identity")  + ggplot2::scale_fill_manual(values = cls) + labs(title = "Expsosure", x = "") +
    theme(axis.text.x = element_text(angle = 90)) +  
    guides(fill=guide_legend(title="Signatures"))
  
  if (!sample_name) {
    p =  p +  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(x = "")
    
  }
  
  if(flip_coord){
    
    p =  p + coord_flip()
  }
  
  p
}

