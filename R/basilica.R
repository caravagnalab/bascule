



basilica <- function(
    x,
    reference_catalogue,
    k_list,
    lr,
    n_steps,
    phi,
    delta,
    groups=NULL,
    input_catalogue=NULL
    ) {

  #obj <- list()
  #obj$beta_fixed <- input_catalogue

  counter <- 1
  while (TRUE) {
    print(paste('iter:', counter))

    obj <- fit(
      x=x,
      k_list=k_list,
      lr=lr,
      n_steps=n_steps,
      groups=groups,
      input_catalogue=input_catalogue
      )

    a <- filter_fixed(x, obj$alpha, obj$beta_fixed, phi)
    print(paste('a:', a))
    b <- filter_denovo(obj$beta_denovo, reference_catalogue, delta)
    print(paste('b:', b))

    if (is.null(a) & is.null(b)) {
      break
    } else if (is.null(a)) {
      input_catalogue <- b
    } else if (is.null(b)) {
      input_catalogue <- a
    } else {
      input_catalogue <- rbind(a, b)
    }

    if (nrow(a)==nrow(input_catalogue) & nrow(b)==0) {
      break
    }
    counter <- counter + 1
  }
  return(obj)
}













