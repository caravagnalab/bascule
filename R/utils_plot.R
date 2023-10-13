reformat_contexts = function(signatures, what) {
  if(what == "SBS") {
    signatures = signatures %>%
      dplyr::mutate(variant=paste0(substr(start=3, stop=3, features), ">", substr(start=5, stop=5, features)),
                    context=paste0(substr(start=1, stop=1, features), '_', substr(start=7, stop=7, features)))
  }

  if(what == "DBS") {
    signatures = signatures %>%
      dplyr::mutate(variant=paste0(substr(start=1, stop=2, features),">NN"),
                    context=substr(start=4, stop=5, features))
  }

  if(what == "ID") {
    signatures = signatures %>%
      dplyr::mutate(features=as.character(features),
                    variant=substr(start=1, stop=nchar(features) - 2, features),
                    context=substr(start=nchar(features), stop=nchar(features), features))
  }

  if(what == "CNV") {
    signatures = signatures %>%
      dplyr::mutate(features=as.character(features),
                    variant=paste0(str_split(features, pattern=":")[[1]][1], ":",
                                   str_split(features, pattern=":")[[1]][2]),
                    context= paste0(str_split(features,pattern=":")[[1]][3]))
  }

  return(signatures)
}



