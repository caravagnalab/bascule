
visData <- function(x) {

  df <- tibble::tibble(
    TargetX = character(),
    InputX = character(),
    N_Samples = numeric(),
    Fixed_TP = numeric(),
    Fixed_FP = numeric(),
    Fixed_TN = numeric(),
    Fixed_FN = numeric(),
    Fixed_TPR = numeric(),
    Fixed_Prec = numeric(),
    Fixed_Rec = numeric(),
    Fixed_Acc = numeric(),
    Denovo_Ratio = numeric(),
    Fitness = numeric(),
    MAE = numeric(),
    MSE = numeric()
  )

  for (i in 1:nrow(x)) {
    input_catalogue <- x[i, ]$Input_Catalogue[[1]]
    inferred_alpha <- x[i, ]$Inf_Exposure[[1]]
    inferred_fixed <- x[i, ]$Inf_Fixed[[1]]
    inferred_denovo <- x[i, ]$Inf_Denovo[[1]]
    expected_fixed <- x[i, ]$Exp_Fixed[[1]]
    expected_denovo <- x[i, ]$Exp_Denovo[[1]]
    m <- x[i, ]$x[[1]]
    mr <- reconstruction_matrix(m, inferred_alpha, inferred_fixed, inferred_denovo)

    # fixed signatures
    if (is.null(input_catalogue)) {inp <- c()} else {inp <- rownames(input_catalogue)}
    if (is.null(expected_fixed)) {exp <- c()} else {exp <- rownames(expected_fixed)}
    if (is.null(inferred_fixed)) {inf <- c()} else {inf <- rownames(inferred_fixed)}
    fixed_TP <- length(intersect(inf, exp))
    fixed_FP <- length(setdiff(inf, exp))
    fixed_TN <- length(setdiff(setdiff(inp, exp), inf))
    fixed_FN <- length(setdiff(exp, inf))
    if (length(exp)==0) {fixed_TPR <- fixed_TP+1} else {fixed_TPR <- fixed_TP / length(exp)}
    fixed_Prec <- fixed_TP / (fixed_TP + fixed_FP)
    fixed_Rec <- fixed_TP / (fixed_TP + fixed_FN)
    fixed_Acc <- (fixed_TP + fixed_TN) / (fixed_TP + fixed_TN + fixed_FP + fixed_FN)

    # denovo signatures
    if (is.null(expected_denovo)) {n_exp_denovo <- 0} else {n_exp_denovo <- nrow(expected_denovo)}
    if (is.null(inferred_denovo)) {n_inf_denovo <- 0} else {n_inf_denovo <- nrow(inferred_denovo)}
    denovo_Ratio <- (n_inf_denovo + 1) / (n_exp_denovo + 1)

    # goodness of fitness
    fitness <- fitness.quality(m, mr)
    mae <- compute.mae(m, mr)
    mse <- compute.mse(m, mr)


    df <- df %>% tibble::add_row(
      TargetX = x[i, ]$TargetX,
      InputX = x[i, ]$InputX,
      N_Samples = x[i, ]$Num_Samples,
      Fixed_TP = fixed_TP,
      Fixed_FP = fixed_FP,
      Fixed_TN = fixed_TN,
      Fixed_FN = fixed_FN,
      Fixed_TPR = fixed_TPR,
      Fixed_Prec = fixed_Prec,
      Fixed_Rec = fixed_Rec,
      Fixed_Acc = fixed_Acc,
      Denovo_Ratio = denovo_Ratio,
      Fitness = fitness,
      MAE = mae,
      MSE = mse
    )
  }

  df$TargetX <- factor(df$TargetX, levels = c("low", "medium", "high"))
  df$InputX <- factor(df$InputX, levels = c("low", "medium", "high"))

  return(df)
}


#-------------------------------------------------------------------------------

reconstruction_matrix <- function(m, alpha, fixed, denovo) {
  # all args are data.frame
  theta <- diag(rowSums(m)) # matrix
  alpha <- theta %*% as.matrix(alpha) # matrix
  beta <- as.matrix(rbind(fixed, denovo)) # matrix

  mr_matrix <- alpha %*% as.matrix(beta)
  mr <- round(as.data.frame(mr_matrix))
  return(mr)
}

#-------------------------------------------------------------------------------

fitness.quality <- function(m, mr) {
  total <- 0
  for (i in 1:nrow(m)) {
    cos <- cosine_sim(m[i,], mr[i, ])
    total <- total + cos
  }
  return(total / nrow(m))
}

#-------------------------------------------------------------------------------

compute.mae <- function(m , mr) {
  mae <- sum(abs(m - mr)) / (dim(m)[1] * dim(m)[2])
  return(mae)
}

#-------------------------------------------------------------------------------

compute.mse <- function(m , mr) {
  mse <- sum((m - mr)^2) / (dim(m)[1] * dim(m)[2])
  return(mse)
}

#-------------------------------------------------------------------------------

denovo.quality <- function(exp, inf) {

  if (length(exp)==0 | length(inf)==0) {
    return(list(matrix=0, ratio=0, cosine=0))
  } else {
    df <- data.frame(matrix(nrow = nrow(inf), ncol = nrow(exp)))
    colnames(df) <- rownames(exp)
    rownames(df) <- rownames(inf)

    for (i in 1:nrow(inf)) {
      inferred <- inf[i,]
      inferred_name <- rownames(inferred)
      maxScore <- 0
      bestMatch <- NULL
      for (j in 1:nrow(exp)) {
        target <- exp[j, ]
        target_name <- rownames(target)
        score <- cosine_sim(inferred, target)
        df[inferred_name, target_name] <- score
        '
        if (score > maxScore) {
          maxScore <- score
          bestMatch <- target_name
        }
        '
      }
    }
    '
    match_list <- colnames(df)[max.col(df, ties.method = "first")]
    s <- 0
    for (i in 1:length(match_list)) {
      s <- s + df[rownames(df)[i], match_list[i]]
    }
    cosine_score <- (s / length(match_list))
    '
  }

  a <- nrow(inf) / nrow(exp)
  b <- mean(apply(df, 1, max))
  return(list(matrix=df, ratio=a, cosine=b))
}

#-------------------------------------------------------------------------------


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

fill_tibble <- function(x) {

  data <- init_tibble()
  for (row in x) {
    data <- data %>% tibble::add_row(
      x = list(row$x),
      Input_Catalogue = list(row$input_catalogue),
      Ref_Catalogue = list(row$ref_catalogue),

      Exp_Exposure = list(row$exp_exposure),
      Exp_Fixed = list(row$exp_fixed),
      Exp_Denovo = list(row$exp_denovo),

      Inf_Exposure = list(row$inf_exposure),
      Inf_Fixed = list(row$inf_fixed),
      Inf_Denovo = list(row$inf_denovo),

      TargetX = row$targetX,
      InputX = row$inputX,
      Num_Samples = row$num_samples,
      IterNum = row$iter,

      K = list(row$k),
      Lr = row$lr,
      Steps = row$steps,
      Phi = row$phi,
      Delta = row$delta
    )
  }
  return(data)
}







