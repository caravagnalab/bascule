
#' @import tibble
#' @import dplyr
#'
#'
generate_data <- function(ref_catalogue,
                          target_complexity,
                          input_complexity,
                          num_samples) {

  pybasilica <- reticulate::import("pybasilica")

  x <- pybasilica$input_generator(
    full_cosmic_df = ref_catalogue,
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
    ref_catalogue,
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
    ref_catalogue,
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
sim_fit_batch <- function(ref_catalogue = COSMIC_catalogue,
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

  counter <- 1

  for (i in targetX_list) {
    for (j in inputX_list) {
      for (n in num_samples_list) {
        for (iter in num_iters) {

          print(paste("iter:", counter))

          try(
            {
              b <- sim_fit(
                ref_catalogue = ref_catalogue,
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
              )
              #---------------------------------------------
            }
          )

          counter <- counter + 1
        }
      }
    }
  }
  return(data)
}


init_tibble <- function() {

  obj <- tibble::tibble(
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

#-------------------------------------------------------------------------------
# split main catalogue to reference and denovo
reference_denovo <- function(ref_path, num_ref, seed) {

  ref_org <- read.table(ref_path, sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
  all_sigs <- rownames(ref_org)

  set.seed(seed = seed)
  ref_list <- sample(all_sigs, num_ref)
  denovo_list <- setdiff(all_sigs, ref_list)

  ref <- ref_org[ref_list, ]
  denovo <- ref_org[denovo_list, ]

  obj <- list(ref=ref, denovo=denovo)
  return(obj)
}


#-------------------------------------------------------------------------------
generate_exposure <- function(signatures, groups, seed) {

  set.seed(seed = seed)
  df_list <- list()

  for (group in unique(groups)) {

    sigNums <- sample(2:length(signatures), 1)
    sigNames <- sample(signatures, sigNums)
    num_samples <- length(groups[groups==group])

    print(paste("group", group, "has", sigNums, "signatures, and", num_samples, "samples"))

    x <- matrix( runif(num_samples * sigNums, 0, 1), ncol = sigNums )
    alpha <- x / rowSums(x)
    alpha <- as.data.frame(alpha)
    colnames(alpha) <- sigNames
    alpha$group <- rep(group, num_samples)
    print(alpha)

    df_list[length(df_list)+1] <- list(alpha)
  }

  data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list) # merge all different group exposure matrices

  # sort columns
  column_names <- colnames(data)
  column_names <- column_names[order(column_names)]
  column_names <- append(setdiff(column_names, "group"), "group")
  #column_names[length(column_names)+1] <- "group"
  data <- data[, column_names]

  data[is.na(data)] <- 0    # convert 'NA' to zero
  data[order(data$group), ] # sort rows by group column

  return(data)

}




#-------------------------------------------------------------------------------
# generate theta vector from counts catalogue
generate_theta <- function(x) {
  theta <- rowSums(x)
  return(theta)
}

#-------------------------------------------------------------------------------
generate_signatures <- function(reference_catalogue, denovo_catalogue, complexity, num_samples, seed) {

  set.seed(seed = seed)

  if (complexity=='low') {
    fixed_num <- sample(3:5, 1)
    denovo_num <- sample(0:2, 1)
  }
  else if (complexity=='medium') {
    fixed_num <- sample(0:2, 1)
    denovo_num <- sample(3:5, 1)
  }
  else if (complexity=='high') {
    fixed_num <- sample(3:5, 1)
    denovo_num <- sample(3:5, 1)
  }
  else {
    stop("wrong complexity!")
  }

  in_reference_list <- rownames(reference_catalogue)
  out_reference_list <- rownames(denovo_catalogue)
  mutation_features <- colnames(reference_catalogue)

  # catalogue signatures -------------------------------------------------------
  if (fixed_num > 0) {
    fixed_list <- sample(in_reference_list, fixed_num)
    fixed_df <- reference_catalogue[fixed_list, ]
  }
  else {
    fixed_df <- NULL
  }

  # denovo signatures ----------------------------------------------------------
  if (denovo_num > 0) {
    denovo_list <- sample(out_reference_list, denovo_num)
    denovo_df <- denovo_catalogue[denovo_list, ]
    rownames(denovo_df) <- paste0(rownames(denovo_df), "_D")
  }
  else {
    denovo_df <- NULL
  }

  if (is.null(fixed_df)) {
    beta <- denovo_df
  }
  else if (is.null(denovo_df)) {
    beta <- fixed_df
  }
  else {
    beta <- rbind(fixed_df, denovo_df)
  }

  return(beta)
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

'
ref_path <- "/home/azad/Documents/thesis/pybasilica/pybasilica/data/cosmic/cosmic_catalogue.csv"
reference_catalogue <- read.table(ref_path, sep = ",", row.names = 1, header = TRUE, check.names = FALSE)

signatures <- rownames(beta)
groups <- c(1,1,3,1,2,1,3,2,2)
seed = 123

a <- reference_denovo(ref_path, num_ref=50, seed=123)
beta <- generate_signatures(a$ref, a$denovo, complexity="medium", num_samples=5, seed=123)
alpha <- generate_exposure(signatures=signatures , groups=groups, seed=123)
'





