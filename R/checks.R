have_groups = function(x) {
  return("clustering" %in% names(x))
}


have_alternatives = function(x) {
  any(sapply(get_alternatives(x, what="nmf"), function(x) !purrr::is_empty(x$fits)))
}
