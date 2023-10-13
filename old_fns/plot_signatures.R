#' #' Plot signatures
#' #'
#' #' @description creates bar plot of inferred signature profiles,
#' #' where x-axis are 96 substitution bases and y-axis are their relative contribution.
#' #'
#' #' @param what add
#' #' @param cls add
#' #' @param x basilica object
#' #'
#' #' @return plot
#' #' @export plot_signatures
#'
#' plot_signatures = function(x, what = "SBS", context = T, cls = NULL,
#'                            signames = NULL, catalogue = NULL) {
#'
#'   if (is.null(signames))
#'     signames = rownames(x %>% get_signatures)
#'   if (!is.null(catalogue))
#'     signames = unique(c(signames, rownames(catalogue)))
#'
#'   a = NULL
#'
#'   if("catalogue_signatures" %in% names(x$fit)) a = rbind(x$fit$catalogue_signatures, a)
#'   if("denovo_signatures" %in% names(x$fit)) a = rbind(x$fit$denovo_signatures, a)
#'
#'   if (!is.null(catalogue))
#'     a = rbind(a, catalogue[setdiff(rownames(catalogue), rownames(a)),])
#'
#'   if(is.null(cls) && !have_color_palette(x)) cls = gen_palette(nrow(a)) %>% setNames(rownames(a))
#'   else if (is.null(cls) && have_color_palette(x)) cls = get_color_palette(x)
#'
#'   return(
#'     plot_signatures_aux(catalogue=a, what="SBS", context=context, cls=cls, signames=signames)
#'   )
#' }
#'



# plot_signatures_aux = function(signatures, what="SBS",
#                                context=T, cls=NULL,
#                                signames=NULL) {
#   if (is.null(signames)) signames = rownames(signatures)
#   if (is.null(cls)) cls = gen_palette(n=nrow(signatures))
#
#   p = signatures %>%
#     ggplot() +
#     geom_bar(aes(value, x=context, fill=sigs), stat="identity") +
#     ggh4x::facet_nested(sigs ~ type + factor(variant, levels=signatures$variant %>% unique()), scales="free") +
#     theme_bw() +
#     scale_fill_manual(values=cls) +
#     theme(axis.text.x=element_text(angle=90)) +
#     guides(fill="none") +
#     labs(y="", x="", title="Signatures")
#
#   if(!context)
#     p = p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) + labs(x="")
#
#   return(p)
# }
