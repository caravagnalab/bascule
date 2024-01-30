plot_signatures = function(x, types=get_types(x), context=T, cls=NULL,
                           signames=get_signames(x)) {

  signatures = lapply(types, function(tid)
    get_signatures(x, types=tid)[[tid]] %>%
      reformat_contexts(what=tid) %>%
      dplyr::mutate(type=tid)) %>%
    do.call(rbind, .) %>% dplyr::filter(sigs %in% unlist(signames))

  if (is.null(cls)) cls = gen_palette(x, types=types)

  lapply(types, function(t) plot_signatures_aux(signatures %>% dplyr::filter(type==t),
                                                what=t, context=context, signames=signames,
                                                cls=cls)) %>%
    patchwork::wrap_plots(ncol=length(types))
}



plot_signatures_aux = function(signatures, what="SBS",
                               context=T, cls=NULL,
                               signames=NULL) {
  if (is.null(signames)) signames = signatures$sigs %>% unique()
  if (is.null(cls)) cls = gen_palette(n=length(unlist(signames)))

  p = signatures %>%
    ggplot() +
    geom_bar(aes(value, x=context, fill=sigs), stat="identity") +
    ggh4x::facet_nested(sigs ~ type + factor(variant, levels=signatures$variant %>% unique()), scales="free") +
    theme_bw() +
    scale_fill_manual(values=cls) +
    theme(axis.text.x=element_text(angle=90), legend.position="bottom") +
    guides(fill="none") +
    labs(y="", x="")

  if(!context)
    p = p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) + labs(x="")

  return(p)
}

