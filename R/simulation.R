

generate_data <- function(ref_path,
                          target_complexity,
                          input_complexity,
                          num_samples) {

  pybasilica <- reticulate::import("pybasilica")

  x <- pybasilica$input_generator(
    cosmic_path = ref_path,
    target_complexity = target_complexity,
    input_complexity = input_complexity,
    num_samples = num_samples
  )

  obj = list()

  # counts (data.frame)
  obj$x <- reticulate::py_to_r(x$M)
  # expected exposure (data.frame)
  obj$exp_exposure <- reticulate::py_to_r(x$alpha)

  # expected catalogue signatures (data.frame / NULL)
  if (is.null(x$beta_fixed)) {
    obj$exp_fixed <- x$beta_fixed
  } else {
    obj$exp_fixed <- reticulate::py_to_r(x$beta_fixed)
  }
  # expected denovo signatures (data.frame / NULL)
  if (is.null(x$beta_denovo)) {
    obj$exp_denovo <- x$beta_denovo
  } else {
    obj$exp_denovo <- reticulate::py_to_r(x$beta_denovo)
  }
  # input catalogue (data.frame / NULL)
  if (is.null(x$beta_input)) {
    obj$input_catalogue <- x$beta_input
  } else {
    obj$input_catalogue <- reticulate::py_to_r(x$beta_input)
  }

  # reference catalogue (data.frame)
  obj$ref_catalogue <- reticulate::py_to_r(x$cosmic_df)

  obj$targetX <- target_complexity  # character
  obj$inputX <- input_complexity    # character
  obj$num_samples <- num_samples    # numeric

  return(obj)
}

#-------------------------------------------------------------------------------
sim_fit <- function(
    ref_path,
    target_complexity,
    input_complexity,
    num_samples,
    k,
    lr,
    steps,
    phi,
    delta,
    iter
) {

  x <- generate_data(
    ref_path,
    target_complexity,
    input_complexity,
    num_samples
  )

  pybasilica <- reticulate::import("pybasilica")
  fit <- pybasilica$pyfit(
    M=x$x,
    groups=NULL,
    input_catalogue=x$input_catalogue,
    reference_catalogue=x$ref_catalogue,
    k=k,
    lr=lr,
    steps=steps,
    phi=phi,
    delta=delta
  )

  x$inf_exposure = fit[[1]] # data.frame
  x$inf_fixed = fit[[2]]    # data.frame
  x$inf_denovo = fit[[3]]   # data.frame

  x$k = k         # integer
  x$lr = lr       # numeric
  x$steps = steps # numeric
  x$phi = phi     # numeric
  x$delta = delta # numeric
  x$iter = iter

  return(x)
}

#-------------------------------------------------------------------------------
sim_fit_batch <- function(ref_path,
                          targetX_list,
                          inputX_list,
                          num_samples_list,
                          num_iters,
                          k,
                          lr,
                          steps,
                          phi,
                          delta) {

  data <- init_tibble()

  for (i in targetX_list) {
    for (j in inputX_list) {
      for (n in num_samples_list) {
        for (iter in num_iters) {
          b <- sim_fit(
            ref_path = ref_path,
            target_complexity = i,
            input_complexity = j,
            num_samples = n,
            k,
            lr,
            steps,
            phi,
            delta,
            iter
          )

          #---------------------------------------------
          data <- data %>% tibble::add_row(
            x = list(b$x),
            Input_Catalogue = list(b$input_catalogue),
            Ref_Catalogue = list(b$ref_catalogue),

            Exp_Exposure = list(b$exp_exposure),
            Exp_Fixed = list(b$exp_fixed),
            Exp_Denovo = list(b$exp_denovo),

            Inf_Exposure = list(b$inf_exposure),
            Inf_Fixed = list(b$inf_fixed),
            Inf_Denovo = list(b$inf_denovo),

            TargetX = b$targetX,
            InputX = b$inputX,
            Num_Samples = b$num_samples,
            IterNum = b$iter,

            K = list(b$k),
            Lr = b$lr,
            Steps = b$steps,
            Phi = b$phi,
            Delta = b$delta
            #---------------------------------------------
          )
        }
      }
    }
  }
  return(data)
}


init_tibble <- function() {

  obj <- tibble(
    x = list(),
    Input_Catalogue = list(),
    Ref_Catalogue = list(),

    Exp_Exposure = list(),
    Exp_Fixed = list(),
    Exp_Denovo = list(),

    Inf_Exposure = list(),
    Inf_Fixed = list(),
    Inf_Denovo = list(),

    TargetX = character(),
    InputX = character(),
    Num_Samples = numeric(),
    IterNum = numeric(),

    K = list(),
    Lr = numeric(),
    Steps = numeric(),
    Phi = numeric(),
    Delta = numeric()
  )

  return(obj)
}

#-------------------------------------------------------------[Waiting for Test]
catalogue_perf <- function(input_catalogue, expected_fixed, inferred_fixed) {
  # all args are data.frame
  inp <- rownames(input_catalogue)
  exp <- rownames(expected_fixed)
  inf <- rownames(inferred_fixed)

  TP <- length(intersect(inf, exp))
  FP <- length(setdiff(inf, exp))
  TN <- length(setdiff(setdiff(inp, exp), inf))
  FN <- length(setdiff(exp, inf))
  Accuracy <- (TP + TN) / (TP + TN + FP + FN)
  return(Accuracy)
}

#-------------------------------------------------------------------------------
reconstruction_matrix <- function(m, alpha, fixed, denovo) {
  # all args are data.frame
  theta <- diag(rowSums(m)) # matrix
  alpha <- theta %*% as.matrix(alpha) # matrix
  beta <- as.matrix(rbind(fixed, denovo)) # matrix

  mr_matrix <- alpha %*% as.matrix(beta)
  mr <- round(as.data.frame(mr_matrix))
}
#----------------------------
