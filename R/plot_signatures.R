#' plot signatures
#'
#' @description creates bar plot of inferred signature profiles,
#' where x-axis are 96 substitution bases and y-axis are their relative contribution.
#'
#' @param what
#' @param cls
#' @param x basilica object
#'
#' @return plot
#' @export plot_signatures
#' @examples


# plot_signatures <- function(x, Type = c("De novo", "Catalogue"), highlight = 0.05)
# {
#   sigs = x %>% get_signatures(long = TRUE) %>%
#     dplyr::filter(Type %in% !!Type)
#
#   return(plot_signatures_sbs(x, sigs, highlight))
# }
#
# plot_signatures_sbs = function(x, sigs, highlight = 0.05)
# {
#   # Remove parenthesis
#   sigs$substitution = stringr::str_extract_all(sigs$Feature, "\\[[^\\]\\[]*]") %>% unlist()
#   sigs$substitution = gsub(pattern = '\\[', replacement = '', x = sigs$substitution)
#   sigs$substitution = gsub(pattern = '\\]', replacement = '', x = sigs$substitution)
#
#   # Detect highlights
#   sigs = sigs %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(
#       context = paste0(substr(Feature, 1, 1), '_', substr(Feature, 7, 7)),
#     ) %>%
#     dplyr::group_by(Signature) %>%
#     dplyr::mutate(
#       cutoff = !!highlight * sum(Value),
#       highlight = ifelse(
#         (Value > !!highlight * sum(Value)),
#         paste0("probability mass above ", !!highlight * 100, '%'),
#         NA
#       )
#     )
#
#   # Extract a crazy map to colour the reference nucleotide
#   sigs$ref = substr(sigs$substitution, 1, 1)
#
#   sigs = sigs %>%
#     rowwise() %>%
#     mutate(crazy_map =
#              gsub(
#                '_',
#                paste0("<span style='color:indianred3'>", ref, "</span>"),
#                context
#              ))
#
#   # Nice plot
#   plt = ggplot2::ggplot(sigs) +
#     ggplot2::geom_hline(
#       data = sigs %>% dplyr::distinct(Signature, cutoff),
#       ggplot2::aes(yintercept = cutoff),
#       size = .2,
#       linetype = 'dashed',
#       color = 'black'
#     )  +
#     ggplot2::geom_bar(
#       ggplot2::aes(x = crazy_map, y = Value, fill = Signature, color = highlight),
#       stat="identity",
#       position="identity") +
#     ggplot2::facet_grid(as.character(Signature)~substitution, scales="free", drop = FALSE, space="free_x") +
#     my_ggplot_theme() +
#     ggplot2::theme(
#       axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
#     ) +
#     ggplot2::labs(
#       y = "Frequency",
#       title = paste0(x$cohort, ' (n = ', x$n_samples, ')')
#     ) +
#     ggplot2::scale_fill_manual(values = get_signature_colors(x) %>% ggplot2::alpha(0.7)) +
#     ggplot2::scale_color_manual(values = "black", na.value = NA, na.translate = F) +
#     ggplot2::guides(fill = ggplot2::guide_legend(
#       paste0("p >", highlight * 100, '%')
#     )) +
#     # ggplot2::guides(
#     #   fill = ggplot2::guide_legend(
#     #     nrow = ifelse(x$n_catalogue + x$n_denovo > 8, 2, 1)),
#     #   color = ggplot2::guide_legend("", override.aes = ggplot2::aes(fill = NA))
#     # ) +
#     ggplot2::theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 0))
#
#   return(plt)
# }


plot_signatures = function(x,what = "SBS", context = T, cls = NULL){

  a = NULL

  if("catalogue_signatures" %in% names(x$fit)){ a = rbind(x$fit$catalogue_signatures,a)}
  if("denovo_signatures" %in% names(x$fit)){ a = rbind(x$fit$denovo_signatures,a)}


  if(is.null(cls)){ cls = ggsci::pal_simpsons()(nrow(a))
  names(cls) = rownames(a)
  }

  a =  a %>% dplyr::mutate(sbs = rownames(a)) %>% as_tibble() %>%
    reshape2::melt()


  if(what == 'SBS'){
    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate(
        substitution = paste0(substr(start = 3, stop = 3, Var2),">",substr(start = 5, stop = 5, Var2)),
        context = paste0(
          substr(start = 1, stop = 1, Var2),
          '_',
          substr(start = 7, stop = 7, Var2)
        )
      )  }

  if(what == "DBS"){

    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate(
        substitution = paste0(substr(start = 1, stop = 2, Var2),">NN"),
        context = substr(start = 4, stop = 5, Var2)
      )
  }

  if(what == "ID"){

    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate( Var2 = as.character(Var2),
              substitution = substr(start = 1, stop = nchar(Var2) - 2, Var2),
              context = substr(start = nchar(Var2), stop = nchar(Var2), Var2)
      )
  }

  if(what == "CNV"){

    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>% rowwise() %>%
      mutate( Var2 = as.character(Var2),
              substitution =
                paste0(str_split(Var2,pattern = ":")[[1]][1],":",str_split(Var2,pattern = ":")[[1]][2]),
              context =  paste0(str_split(Var2,pattern = ":")[[1]][3])
      )
  }

  # library(CNAqc)

  p = a  %>%
    ggplot() +
    geom_bar(aes(value, x = context, fill = Var1), stat = 'identity') +
    facet_grid(Var1 ~ factor(substitution,levels = a$substitution %>% unique()), scales = 'free') +
    # CNAqc:::my_ggplot_theme() +
    theme_bw() +
    scale_fill_manual(values = cls) +
    theme(axis.text.x = element_text(angle = 90)) +
    guides(fill = 'none')  +
    labs(y = "", title = "Signatures")

  if(!context) {
    p = p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(x = "")
  }

  p
}

plot_exposure = function(x,sample_name = T,levels= NULL, flip_coord = F){

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


