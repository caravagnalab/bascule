have_groups = function(x) {
  return("clustering" %in% names(x))
}


have_alternatives = function(x) {
  list(
    "nmf"=sapply(get_types(x), function(tid)
      !purrr::is_empty(x$nmf[[tid]]$pyro$alternatives)),
    "clustering"=!purrr::is_empty(x$clustering$pyro$alternatives)
    )
}
