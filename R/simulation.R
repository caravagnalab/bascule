

#----------------------------------------------------------------------QC:PASSED
# split reference catalogue to 2 sub catalogue:
# reference catalogue (SBS1 included) + denovo catalogue
split.reference <- function(reference, ratio, seed=NULL) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  num_ref <- round(ratio * nrow(reference))

  SBS1 <- reference['SBS1', ]   # save SBS1 (data.frame)
  reference <- reference[!(rownames(reference) %in% c("SBS1")), ] # excludes SBS1

  # suffle the reference catalogue
  shuffled_reference = reference[sample(1:nrow(reference)), ]

  ref <- shuffled_reference[1:(num_ref-1), ]
  ref <- ref[order(rownames(ref)), ]
  ref <- rbind(SBS1, ref) # includes SBS1

  denovo <- shuffled_reference[num_ref:nrow(shuffled_reference), ]
  denovo <- denovo[order(rownames(denovo)), ]

  obj <- list(reference=ref, denovo=denovo)
  return(obj)
}

#----------------------------------------------------------------------QC:PASSED

generate.theta <- function(mut_range, num_samples, seed=NULL) {

  if (!(is.integer(mut_range))) {
    stop("not valid range argument in generate_theta function!")
  }

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  theta = sample(mut_range, num_samples)  # integer
  return(theta)
}

#----------------------------------------------------------------------QC:PASSED
# generate signatures which includes:
#   * fixed signatures (SBS1 included)
#   * denovo signatures
# similarity between catalogue signatures are less than threshold
# similarity between denovo signatures are less than threshold
# but similarity between catalogue and denovo signatures are not taken to the account (future work)
generate.signatures <- function(
    reference_catalogue,
    denovo_catalogue,
    reference_cosine, # cosine similarity matrix of reference signatures (SBS1 excluded)
    denovo_cosine,    # cosine similarity matrix of denovo signatures
    complexity,
    similarity_limit,
    seed=NULL
    ) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  if (is.numeric(complexity) & length(complexity)==2) {
    fixed_num <- complexity[1]
    denovo_num <- complexity[2]
  } else if (is.character(complexity)) {
    if (complexity=='low') {
      fixed_num <- sample(3:5, 1)
      denovo_num <- sample(0:2, 1)
    }
    else if (complexity=='medium') {
      fixed_num <- sample(1:2, 1)
      denovo_num <- sample(3:5, 1)
    }
    else if (complexity=='high') {
      fixed_num <- sample(3:5, 1)
      denovo_num <- sample(3:5, 1)
    }
    else {
      stop("complexity argument should be selected from {'low', 'medium', 'high'}")
    }
  } else {
    stop("wrong complexity argument!")
  }

  SBS1 <- reference_catalogue['SBS1', ] # save SBS1 (data.frame)
  reference <- reference_catalogue[!(rownames(reference_catalogue) %in% c("SBS1")), ] # excludes SBS1

  # catalogue signatures -------------------------------------------------------

  if (fixed_num > 1) {

    while (TRUE) {
      shuffled_reference = reference[sample(1:nrow(reference)), ]
      signatures <- rownames(shuffled_reference[1:(fixed_num-1), ])
      cos_matrix <- reference_cosine[c("SBS1", signatures), c("SBS1", signatures)]
      for (i in 1:nrow(cos_matrix)) {
        cos_matrix[i, i] <- 0
      }
      max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
      if (cos_matrix[max][1] < similarity_limit) {
        fixed_df <- rbind(SBS1, reference[signatures, ])
        break
      }
    }
  }
  else {
    fixed_df <- SBS1
  }

  # denovo signatures ----------------------------------------------------------

  if (denovo_num > 1) {

    while (TRUE) {
      shuffled_denovo = denovo_catalogue[sample(1:nrow(denovo_catalogue)), ]
      signatures <- rownames(shuffled_denovo[1:denovo_num, ])
      cos_matrix <- denovo_cosine[signatures, signatures]
      for (i in 1:nrow(cos_matrix)) {
        cos_matrix[i, i] <- 0
      }
      max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
      if (cos_matrix[max][1] < similarity_limit) {
        denovo_df <- denovo_catalogue[signatures, ]
        rownames(denovo_df) <- paste0(rownames(denovo_df), "_D")
        break
      }
    }
  }
  else if (denovo_num==1) {
    shuffled_denovo = denovo_catalogue[sample(1:nrow(denovo_catalogue)), ]
    denovo_df <- shuffled_denovo[1, ]
    rownames(denovo_df) <- paste0(rownames(denovo_df), "_D")
  }
  else {
    denovo_df <- NULL
  }

  obj <- list(fixed = fixed_df, denovo = denovo_df)
  return(obj)
}

#----------------------------------------------------------------------QC:PASSED

generate.input <- function(
    reference_catalogue,
    beta_fixed,
    complexity,
    seed=NULL
    ) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  k_fixed <- nrow(beta_fixed)

  if (is.null(complexity)) {
    return(NULL)
  } else if (is.numeric(complexity) & length(complexity)==2) {
    n_overlap <- complexity[1]
    n_extra <- complexity[2]
    if (n_overlap > k_fixed) {
      stop("overlap signatures are more than fixed signatures!")
    }
  } else if (is.character(complexity)) {
    if (complexity=='low') {
      n_overlap <- sample(1:k_fixed, 1)
      n_extra <- 0
    }
    else if (complexity=='medium') {
      n_overlap <- sample(1:k_fixed, 1)
      n_extra <- sample(1:k_fixed, 1)
    }
    else if (complexity=='high') {
      n_overlap <- 0
      n_extra <- sample(1:(k_fixed), 1)
    }
    else {
      stop("complexity argument should be selected from {'low', 'medium', 'high'}")
    }
  } else {
    stop("wrong complexity argument!")
  }

  extra_ref <- reference_catalogue[setdiff(rownames(reference_catalogue), rownames(beta_fixed)), ]

  if (n_overlap > 0) {
    overlap <- sample(rownames(beta_fixed))[1:n_overlap]
  } else {
    overlap <- NULL
  }

  if (n_extra > 0) {
    extra <- sample(rownames(extra_ref))[1:n_extra]
  } else {
    extra <- NULL
  }

  df <- reference_catalogue[c(overlap, extra), ]

  return(df)
}

#----------------------------------------------------------------------QC:PASSED

generate.exposure <- function(beta, groups, seed=NULL) {

  signatures <- rownames(beta)
  if (!('SBS1' %in% signatures)) {
    stop('Wrong signatures! SBS1 not included!')
  }

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  if (length(signatures) < 2) {
    stop("not valid! there are not enough signatures!")
  }

  df_list <- list()

  signatures <- signatures[! signatures %in% c('SBS1')] # excludes SBS1

  for (group in unique(groups)) {

    if (length(unique(groups))==1) {
      sigNums <- length(signatures)
      sigNames <- c('SBS1', signatures)
    } else {
      sigNums <- sample(1:length(signatures), 1)
      sigNames <- c('SBS1', sample(signatures, sigNums))
    }
    #sigNums <- sample(1:length(signatures), 1)
    #sigNames <- c('SBS1', sample(signatures, sigNums))

    num_samples <- length(groups[groups==group])

    #print(paste("group", group, "has", sigNums+1, "signatures, and", num_samples, "samples"))

    x <- matrix( runif(num_samples * (sigNums+1), 0, 1), ncol = sigNums+1 )
    alpha <- x / rowSums(x)
    alpha <- as.data.frame(alpha)
    colnames(alpha) <- sigNames
    alpha$group <- rep(group, num_samples)
    #print(alpha)

    df_list[length(df_list)+1] <- list(alpha)
  }

  # merge all different group exposure matrices
  data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

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

#----------------------------------------------------------------------QC:PASSED
# may use theta*alpha*beta or generate.counts (should be discussed to )
generate.counts <- function(alpha, beta, theta, seed=NULL) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  if (!is.null(alpha$group)) {
    alpha <- subset(alpha, select = -c(group))
  }

  alpha <- alpha[, order(colnames(alpha))]
  beta <- beta[order(rownames(beta)), ]
  if (!identical(colnames(alpha), rownames(beta))) {
    print(colnames(alpha))
    print(rownames(beta))
    stop("alpha and beta are NOT valid!")
  }

  num_samples <- nrow(alpha)

  M <- matrix(rep(0, num_samples*96) , ncol = 96)

  # iterate over samples
  for (sample in 1:num_samples) {

    p <- as.numeric(alpha[sample, ]) # select sample i

    # iterate over number of mutations in sample i
    for (j in 1:theta[sample]) {

      # sample signature profile index from categorical data
      signature_idx <- extraDistr::rcat(1, p)
      signature <- beta[signature_idx, ]

      # sample mutation feature index for corresponding signature from categorical data
      mutation_idx <- extraDistr::rcat(1, as.numeric(signature))

      # add +1 to the mutation feature in position j in branch i
      M[sample, mutation_idx] <- M[sample, mutation_idx] + 1

    }
  }
  M <- as.data.frame(M)
  colnames(M) <- colnames(beta)
  rownames(M) <- rownames(alpha)
  return(M)
}

#----------------------------------------------------------------------QC:PASSED

generate.data <- function(
    reference_catalogue,
    denovo_catalogue,
    reference_cosine,
    denovo_cosine,
    targetX,
    inputX,
    similarity_limit,
    groups,
    mut_range,
    seed=NULL
    ) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  # SIGNATURES --------------------------------
  signatures <- generate.signatures(
    reference_catalogue = reference_catalogue,
    denovo_catalogue  = denovo_catalogue,
    reference_cosine = reference_cosine,
    denovo_cosine = denovo_cosine,
    complexity = targetX,
    similarity_limit = similarity_limit,
    seed = seed
  )
  beta <- rbind(signatures$fixed, signatures$denovo)

  # INPUT ------------------------------------
  input <- generate.input(
    reference_catalogue = reference_catalogue,
    beta_fixed = signatures$fixed,
    complexity = inputX,
    seed = seed
  )

  # EXPOSURE ---------------------------------
  alpha <- generate.exposure(beta=beta, groups=groups, seed=seed) # include group column

  # THETA ------------------------------------
  num_samples <- length(groups)
  theta <- generate.theta(mut_range=mut_range, num_samples=num_samples, seed=seed)

  # COUNT MATRIX -----------------------------
  #m <- generate.counts(alpha=alpha, beta=beta, theta=theta, seed=seed)

  # removing group column
  if (!is.null(alpha$group)) {
    alpha <- subset(alpha, select = -c(group))
  }
  M <- as.data.frame(round(as.matrix(alpha*theta) %*% as.matrix(beta), digits = 0))
  rownames(M) <- rownames(alpha)
  colnames(M) <- colnames(beta)

  # MODIFY COMPLEXITY VALUES -----------------
  if (is.numeric(targetX)) {
    targetX <- paste(targetX, collapse = "|")
  }
  if (is.numeric(inputX)) {
    inputX <- paste(inputX, collapse = "|")
  }

  # CREATE TIBBLE ----------------------------
  obj <- tibble::tibble(
    x = list(M),
    input_cat = list(input),
    ref_cat = list(reference_catalogue),
    exp_exposure = list(alpha),
    exp_fixed = list(signatures$fixed),
    exp_denovo = list(signatures$denovo),
    targetX = targetX,
    inputX = inputX,
  )
  return(obj)
}

#----------------------------------------------------------------------QC:PASSED

run.data <- function(
    data,
    k,
    lr = 0.01,
    steps = 500,
    max_iterations = 20,
    blacklist = NULL,
    phi = 0.05,
    delta = 0.9,
    filt_pi =0.1,
    groups = NULL,
    lambda_rate = NULL,
    sigma = FALSE,
    CUDA = FALSE,
    compile = TRUE,
    enforce_sparsity = FALSE,
    input=TRUE
) {

  x <- data$x[[1]]
  reference <- data$ref_cat[[1]]
  if (input) {
    input <- data$input_cat[[1]]
  } else {
    input = NULL
  }

  obj <- basilica::fit(
    x=x,
    reference_catalogue = reference,
    k = k,
    cohort = "Simulation",
    lr = lr,
    steps = steps,
    max_iterations = max_iterations,
    blacklist = blacklist,
    phi = phi,
    delta = delta,
    filt_pi =filt_pi,
    groups = groups,
    lambda_rate = lambda_rate,
    sigma = sigma,
    CUDA = CUDA,
    compile = compile,
    enforce_sparsity = enforce_sparsity
  )

  simulation.fit.obj <- tibble::add_column(
    data,
    #inf_exposure = list(obj$exposure),
    #inf_denovo = list(obj$denovo_signatures),
    #inf_fixed = list(obj$catalogue_signatures),
    #bic = obj$bic,
    #losses = list(obj$losses),
    fit = list(obj)
  )
  return(simulation.fit.obj)
}

#----------------------------------------------------------------------QC:PASSED

generate.cohort <- function(
    reference_signatures,
    ratio,
    targetX,
    inputX,
    similarity_limit,
    groups,
    mut_range,
    seed = NULL,
    num_data
    ) {

  reference_denovo <- split.reference(reference =reference_signatures, ratio=ratio, seed=seed)

  reference_catalogue <- reference_denovo$reference
  denovo_catalogue <- reference_denovo$denovo

  reference_cosine <- basilica:::cosine.matrix(reference_catalogue, reference_catalogue)
  denovo_cosine <- basilica:::cosine.matrix(denovo_catalogue, denovo_catalogue)

  data <- NULL

  if (is.null(seed)) {
    for (i in 1:num_data) {
      xx <- basilica:::generate.data(
        reference_catalogue = reference_catalogue,
        denovo_catalogue = denovo_catalogue,
        reference_cosine = reference_cosine,
        denovo_cosine = denovo_cosine,
        targetX = targetX,
        inputX = inputX,
        similarity_limit = similarity_limit,
        groups = groups,
        mut_range = mut_range,
        seed = seed
      )
      data <- rbind(data, xx)
    }
  } else {
    for (i in 1:num_data) {
      xx <- basilica:::generate.data(
        reference_catalogue = reference_catalogue,
        denovo_catalogue = denovo_catalogue,
        reference_cosine = reference_cosine,
        denovo_cosine = denovo_cosine,
        targetX = targetX,
        inputX = inputX,
        similarity_limit = similarity_limit,
        groups = groups,
        mut_range = mut_range,
        seed = seed
      )
      data <- rbind(data, xx)
      seed <- seed + 1
    }
  }

  return(data)
}

#----------------------------------------------------------------------QC:PASSED

run.cohort <- function(
    cohort,
    k = 0:5,
    lr = 0.01,
    steps = 500,
    max_iterations = 20,
    blacklist = NULL,
    phi = 0.05,
    delta = 0.9,
    filt_pi =0.1,
    groups = NULL,
    lambda_rate = NULL,
    sigma = FALSE,
    CUDA = FALSE,
    compile = TRUE,
    enforce_sparsity = FALSE,
    input=TRUE
    ) {

  results <- NULL
  for (i in 1:nrow(cohort)) {

    cat('============================================\n') # TEST
    cat('                 Data No.', i, '\n') # TEST
    cat('============================================\n') # TEST

    xx <- basilica:::run.data(
      data = cohort[i, ],
      k = k,
      lr = lr,
      steps = steps,
      max_iterations = max_iterations,
      blacklist = blacklist,
      phi = phi,
      delta = delta,
      filt_pi = filt_pi,
      groups = groups,
      lambda_rate = lambda_rate,
      sigma = sigma,
      CUDA = CUDA,
      compile = compile,
      enforce_sparsity = enforce_sparsity,
      input=input
    )
    results <- rbind(results, xx)
  }
  return(results)
}

#===============================================================================
#=========================== EVALUATION ========================================
#===============================================================================


#----------------------------------------------------------------------QC:PASSED

fixed.accuracy <- function(reference, input, expected_fixed, inferred_fixed) {
  ref_list <- rownames(reference)
  if (is.null(expected_fixed)) {exp_list <- c()} else {exp_list <- rownames(expected_fixed)}
  if (is.null(inferred_fixed)) {inf_list <- c()} else {inf_list <- rownames(inferred_fixed)}
  if (is.null(input)) {input_list <- c()} else {input_list <- rownames(input)}

  TP <- length(intersect(inf_list, exp_list))
  FP <- length(setdiff(inf_list, exp_list))
  #TN <- length( setdiff( setdiff(ref_list, exp_list), inf_list) )
  TN <- length( setdiff(setdiff(input_list, exp_list), inf_list)  )
  FN <- length(setdiff(exp_list, inf_list))

  accuracy <- list(TP=TP, FP=FP, TN=TN, FN=FN)

  #accuracy <- (TP + TN) / (TP + TN + FP + FN)
  return(accuracy)
}

#----------------------------------------------------------------------QC:PASSED

reconstruct.count <- function(m, alpha, beta) {
  # all args are data.frame
  theta <- diag(rowSums(m))               # matrix
  alpha <- theta %*% as.matrix(alpha)     # matrix
  beta <- as.matrix(beta)                 # matrix

  mr_matrix <- alpha %*% beta
  mr <- round(as.data.frame(mr_matrix))
  rownames(mr) <- rownames(m)
  return(mr)
}

#----------------------------------------------------------------------QC:PASSED

compute.mae <- function(m , mr) {
  mae <- sum(abs(m - mr)) / (dim(m)[1] * dim(m)[2])
  return(mae)
}

#----------------------------------------------------------------------QC:PASSED

compute.mse <- function(m , mr) {
  mse <- sum((m - mr)^2) / (dim(m)[1] * dim(m)[2])
  return(mse)
}

#----------------------------------------------------------------------QC:PASSED

denovo.similarity <- function(expected_denovo, inferred_denovo) {

  if (length(expected_denovo)==0 | length(inferred_denovo)==0) {
    return(NULL)
  } else {
    df <- data.frame(matrix(nrow = nrow(inferred_denovo), ncol = nrow(expected_denovo)))
    colnames(df) <- rownames(expected_denovo)
    rownames(df) <- rownames(inferred_denovo)

    for (i in 1:nrow(inferred_denovo)) {
      inferred <- inferred_denovo[i,]
      inferred_name <- rownames(inferred)
      for (j in 1:nrow(expected_denovo)) {
        target <- expected_denovo[j, ]
        target_name <- rownames(target)
        score <- cosine.vector(inferred, target)
        df[inferred_name, target_name] <- score
      }
    }

    #------------------------------
    #match_list <- list()
    match_df <- data.frame(matrix(nrow = nrow(inferred_denovo), ncol = 2))
    colnames(match_df) <- c("match", "similarity")
    rownames(match_df) <- rownames(inferred_denovo)

    similarity <- 0
    iter <- min(nrow(inferred_denovo), nrow(expected_denovo))
    for (i in 1:iter) {

      max = which(df == max(df), arr.ind = TRUE)
      similarity <- similarity + df[max]

      row <- row.names(df[max[,1],])
      column <- names(df[max[,2]])

      #match_list[row] <- column
      match_df[row, 'match'] <- column
      match_df[row, 'similarity'] <- df[max]

      df[row, column] <- 0
    }

  return( list( similarity_average=(similarity / iter), match_df=match_df ) )
  }
}

#----------------------------------------------------------------------QC:PASSED

denovo.ratio <- function(expected_denovo, inferred_denovo) {

  if (is.null(expected_denovo)) {n_exp <- 0} else {n_exp <- nrow(expected_denovo)}
  if (is.null(inferred_denovo)) {n_inf <- 0} else {n_inf <- nrow(inferred_denovo)}

  denovo_ratio <- (n_inf + 1) / (n_exp + 1)

  return(denovo_ratio)
}

#----------------------------------------------------------------------QC:PASSED

#' @import dplyr
evaluate.data <- function(x) {
  #--------------------------
  reference <- x$ref_cat[[1]]
  input <- x$input_cat[[1]]
  expected_fixed <- x$exp_fixed[[1]]
  #inferred_fixed <- x$inf_fixed[[1]]
  inferred_fixed <- x$fit[[1]]$fit$catalogue_signatures
  a <- basilica:::fixed.accuracy(reference, input, expected_fixed, inferred_fixed)
  TP <- a$TP
  FP <- a$FP
  TN <- a$TN
  FN <- a$FN
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  #--------------------------
  m <- x$x[[1]]
  #alpha <- x$inf_exposure[[1]]
  alpha <- x$fit[[1]]$fit$exposure
  beta <- rbind(x$fit[[1]]$fit$catalogue_signatures, x$fit[[1]]$fit$denovo_signatures)
  mr <- basilica:::reconstruct.count(m, alpha, beta)
  mae <- basilica:::compute.mae(m, mr)
  mse <- basilica:::compute.mse(m, mr)
  #--------------------------
  b <- basilica:::denovo.similarity(x$exp_denovo[[1]], x$fit[[1]]$fit$denovo_signatures)
  denovo_similarity <- b$similarity_average  # numeric
  denovo_match <- b$match_df                 # data.frame
  #--------------------------
  denovo_ratio <- basilica:::denovo.ratio(x$exp_denovo[[1]], x$fit[[1]]$fit$denovo_signatures)
  #--------------------------

  # CREATE TIBBLE ----------------------------
  obj <- tibble::tibble(

    targetX = x$targetX,
    inputX = x$inputX,
    num_samples = nrow(m),

    mae = mae,
    mse = mse,
    fixed_acc = accuracy,
    denovo_ratio = denovo_ratio,
    denovo_sim = list(denovo_similarity),
    denovo_match = list(denovo_match),
  )

  return(obj)
}

#----------------------------------------------------------------------QC:PASSED

evaluate.cohort <- function(x) {
  res <- NULL
  counter <- 1 # TEST
  for (i in 1:nrow(x)) {
    res <- rbind(res, evaluate.data(x = x[i, ]))
    #print(paste('counter:', counter)) # TEST
    counter <- counter + 1 # TEST
  }
  return(res)
}




