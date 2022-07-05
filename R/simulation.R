

#----------------------------------------------------------------------QC:PASSED
# split reference catalogue to 2 sub catalogue:
# reference catalogue (SBS1 included) + denovo catalogue
split_reference <- function(ref_path, ratio, seed=NULL) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  # read csv file as data.frame
  reference <- read.table(ref_path, sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
  num_ref <- round(ratio * nrow(reference))

  # save SBS1
  SBS1 <- reference['SBS1', ]
  # excludes SBS1 from reference catalogue
  reference <- reference[!(rownames(reference) %in% c("SBS1")), ]

  # suffle the reference catalogue
  shuffled_reference = reference[sample(1:nrow(reference)), ]

  ref <- shuffled_reference[1:(num_ref-1), ]
  ref <- ref[order(rownames(ref)), ]
  ref <- rbind(SBS1, ref) # includes SBS1

  denovo <- shuffled_reference[num_ref:nrow(shuffled_reference), ]
  denovo <- denovo[order(rownames(denovo)), ]

  obj <- list(ref=ref, denovo=denovo)
  return(obj)
}

#----------------------------------------------------------------------QC:PASSED

similarity_tables <- function(reference, denovo) {
  #ref <- reference[!(rownames(reference) %in% c("SBS1")), ] # excludes SBS1
  cosine_reference <- basilica::cosine_matrix(reference, reference)
  cosine_denovo <- basilica::cosine_matrix(denovo, denovo)

  obj <- list(cosine_reference=cosine_reference, cosine_denovo=cosine_denovo)
  return(obj)
}

#----------------------------------------------------------------------QC:PASSED
# generate signatures which includes:
# fixed signatures (SBS1 included) + denovo signatures
generate_signatures <- function(
    reference_catalogue,
    denovo_catalogue,
    reference_cosine, # cosine similarity matrix of reference signatures (SBS1 excluded)
    denovo_cosine, # cosine similarity matrix of denovo signatures
    complexity,
    limit,
    seed=NULL
    ) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

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
    stop("wrong complexity argument!")
  }

  SBS1 <- reference_catalogue['SBS1', ]
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
      #print(paste("fixed cos_matrix[max]", cos_matrix[max][1]))
      if (cos_matrix[max][1] < limit) {
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
      #print(paste("denovo cos_matrix[max]", cos_matrix[max][1]))
      if (cos_matrix[max][1] < limit) {
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

  #if (is.null(denovo_df)) {
  #  beta <- fixed_df
  #}
  #else {
  #  beta <- rbind(fixed_df, denovo_df)
  #}

  obj <- list(fixed = fixed_df, denovo = denovo_df)
  return(obj)
}

#----------------------------------------------------------------------QC:PASSED

generate_input <- function(
    reference_catalogue,
    beta_fixed,
    complexity,
    seed=NULL
    ) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  k_fixed <- nrow(beta_fixed)

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
    stop("wrong complexity argument!")
  }

  extra_ref <- reference_catalogue[setdiff(rownames(reference_catalogue), rownames(beta_fixed)), ]

  if (n_overlap > 0) {
    overlap <- sample(rownames(beta_fixed))[1:n_overlap]
    #shuffled_fixed = beta_fixed[sample(1:nrow(beta_fixed)), ]
    #overlap <- shuffled_fixed[1:n_overlap, ]
  } else {
    overlap <- NULL
  }

  if (n_extra > 0) {
    extra <- sample(rownames(extra_ref))[1:n_extra]
    #shuffled_extra = extra_ref[sample(1:nrow(extra_ref)), ]
    #extra <- shuffled_extra[1:n_overlap, ]
  } else {
    extra <- NULL
  }

  df <- reference_catalogue[c(overlap, extra), ]

  return(df)
}

#----------------------------------------------------------------------QC:PASSED

generate_exposure <- function(beta, groups, seed=NULL) {

  signatures <- rownames(beta)
  if (!('SBS1' %in% signatures)) {
    stop('Wrong signatures! SBS1 not included!')
  }

  if (!is.null(seed)) {
    set.seed(seed = seed)
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

generate_theta <- function(range, num_samples, seed=NULL) {

  if (!(is.integer(range))) {
    stop("not valid range argument in generate_theta function!")
  }

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  theta = sample(range, num_samples)
  return(theta)
}

#--------------------------------------- may be not needed (should be discussed)

generate_counts <- function(alpha, beta, theta, seed=NULL) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  alpha <- subset(alpha, select = -c(group))
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

#-------------------------------------------------------------------------------

generate.one.synthetic <- function(
    reference_catalogue,
    denovo_catalogue,
    reference_cosine,
    denovo_cosine,
    targetX,
    inputX,
    limit,
    groups,
    range,
    seed=NULL
    ) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  b <- generate_signatures(
    reference_catalogue = reference_catalogue,
    denovo_catalogue = denovo_catalogue,
    reference_cosine = reference_cosine, # cosine similarity matrix of reference signatures (SBS1 excluded)
    denovo_cosine = denovo_cosine, # cosine similarity matrix of denovo signatures
    complexity = targetX,
    limit = limit,
    seed = seed
  )

  beta <- rbind(b$fixed, b$denovo)

  input <- generate_input(
    reference_catalogue=reference_catalogue,
    beta_fixed=b$fixed,
    complexity=inputX,
    seed=seed
  )

  alpha <- generate_exposure(beta=beta, groups=groups, seed=seed) # include group column

  num_samples <- length(groups)
  theta <- generate_theta(range=range, num_samples=num_samples, seed=seed)

  # generate count matrix
  # M <- generate_counts(alpha, beta, theta)
  alpha <- subset(alpha, select = -c(group))  # removing group column
  M <- as.data.frame(round(as.matrix(alpha*theta) %*% as.matrix(beta), digits = 0))
  rownames(M) <- rownames(alpha)
  colnames(M) <- colnames(beta)

  #c('x', 'input_cat', 'ref_cat', 'exp_exposure', 'exp_fixed', 'exp_denovo', 'targetX', 'inputX')
  obj <- tibble::tibble(
    x = list(M),
    input_cat = list(input),
    ref_cat = list(reference_catalogue),
    exp_exposure = list(alpha),
    exp_fixed = list(b$fixed),
    exp_denovo = list(b$denovo),
    targetX = targetX,
    inputX = inputX,
  )

  #obj <- list(m=M, alpha=alpha, beta=beta, theta=theta)
  return(obj)
}

#-------------------------------------------------------------------------------

generate.full.synthetic <- function(
    ref_path,
    ratio,
    num_iter,
    targetX,
    inputX,
    limit,
    groups,
    range,
    seed=NULL
) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  a <- basilica::split_reference(ref_path = ref_path, ratio = ratio, seed = seed)
  ref_cat <- a$ref
  denovo_cat <- a$denovo

  c <- basilica::similarity_tables(ref_cat, denovo_cat)
  ref_cosine <- c$cosine_reference
  denovo_cosine <- c$cosine_denovo

  data <- NULL
  for (i in 1:num_iter) {
    xx <- basilica::generate.one.synthetic(
      reference_catalogue = ref_cat,
      denovo_catalogue = denovo_cat,
      reference_cosine=ref_cosine,
      denovo_cosine=denovo_cosine,
      targetX=targetX,
      inputX=inputX,
      limit=limit,
      groups=groups,
      range=range,
      seed=seed
    )
    data <- rbind(data, xx)
  }
  return(data)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

get.one.sample <- function(
    synthetic,
    k,
    lr,
    steps,
    phi,
    delta
    ) {

  x <- synthetic$x[[1]]
  ref <- synthetic$ref_cat[[1]]
  input <- synthetic$input_cat[[1]]

  obj <- basilica:::fit(
    x=x,
    reference_catalogue=ref,
    k=k,
    lr=lr,
    steps=steps,
    phi=phi,
    delta=delta,
    groups=NULL,
    input_catalogue=input
  )

  #obj$exposure <- list(obj$exposure)
  #obj$denovo_signatures <- list(obj$denovo_signatures)
  #obj$catalogue_signatures <- list(obj$catalogue_signatures)
  #results <- c(synthetic, obj)

  obj <- tibble::add_column(
    synthetic,
    inf_exposure = list(obj$exposure),
    inf_denovo = list(obj$denovo_signatures),
    inf_fixed = list(obj$catalogue_signatures)
    )

  return(obj)
}

#-------------------------------------------------------------------------------

get.all.samples <- function(
    synthetic,
    k,
    lr,
    steps,
    phi,
    delta) {

  df <- NULL
  for (i in 1:nrow(synthetic)) {
    syn <- synthetic[i, ]
    output <- get.one.sample(
      syn,
      k,
      lr,
      steps,
      phi,
      delta
    )
    df <- rbind(df, output)
  }

  return(df)
}


# plot exposure-----------------------------------------------------------------

plot.alpha <- function(exp_alpha, inf_alpha) {
  #exp_alpha <- x$exp_exposure[[1]]
  rownames(exp_alpha) <- rownames(inf_alpha)  # just to be consistent, should be fixed later
  exp_alpha$sample <- rownames(exp_alpha)
  exp_alpha_long <- tidyr::gather(exp_alpha, key="signature", value="exposure", c(-sample))
  exp_alpha_long$type <- rep('expected', each=nrow(exp_alpha_long))

  #inf_alpha <- x$exposure[[1]]
  inf_alpha$sample <- rownames(inf_alpha)
  inf_alpha_long <- tidyr::gather(inf_alpha, key="signature", value="exposure", c(-sample))
  inf_alpha_long$type <- rep('inferred', each=nrow(inf_alpha_long))

  alpha <- rbind(exp_alpha_long, inf_alpha_long)

  plt <- ggplot(data = alpha, aes(x=sample, y=exposure, fill=signature)) +
    geom_bar(stat = "identity") +
    facet_grid(type ~ .) +
    theme_minimal() +
    ggtitle("Signatures exposure (Expected vs. Inferred)")
  #scale_y_continuous(labels=scales::percent)

  #glist <-
  #glist[[1]] <- plt

  return(list(plt))
}


# plot signatures --------------------------------------------------------------

plot.beta <- function(beta) {

  if (is.null(beta)) {
    p <- ggplot() +
      theme_void() +
      geom_text(aes(0,0,label='N/A')) +
      xlab(NULL) #optional, but safer in case another theme is applied later
  } else {
    # separate context and alteration
    x <- data.table::as.data.table(reshape2::melt(as.matrix(beta),varnames=c("signature","cat")))
    x[, Context := paste0(substr(cat,1,1), ".", substr(cat, 7, 7)) ]
    x[, alt := paste0(substr(cat,3,3),">",substr(cat,5,5)) ]

    # make the ggplot2 object
    glist <- list()
    for(i in 1:nrow(beta)) {

      plt <- ggplot(x[signature==rownames(beta)[i]]) +
        geom_bar(aes(x=Context,y=value,fill=alt),stat="identity",position="identity") +
        facet_wrap(~alt,nrow=1,scales="free_x") +
        theme(axis.text.x=element_text(angle=90,hjust=1),panel.background=element_blank(),axis.line=element_line(colour="black")) +
        ggtitle(rownames(beta)[i]) + theme(legend.position="none") + ylab("Frequency of mutations")

      #if(!xlabels) {
      plt <- plt + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
      #}

      glist[[i]] <- plt
    }

    # make the final plot
    #gridExtra::grid.arrange(grobs=glist,ncol=ceiling(nrow(beta)/3))

    p <- ggpubr::ggarrange(plotlist=glist, ncol = 1)
  }

  return(list(p))
}

# plot signatures cosine matrix ------------------------------------------------


plot.beta.cosine <- function(exp_denovo, inf_denovo) {

  if (is.null(exp_denovo) | is.null(inf_denovo)) {
    cplot <- ggplot() +
      theme_void() +
      geom_text(aes(0,0,label='N/A')) +
      xlab(NULL) #optional, but safer in case another theme is applied later
  } else {
    # expected vs inferred signatures cosine similarity matrix
    cos <- cosine_matrix(inf_denovo, exp_denovo)
    cos1 <- tibble::rownames_to_column(cos, var = 'inferred_denovo')
    cos_long <- tidyr::gather(cos1, key="expected_denovo", value="cosine_similarity", c(-inferred_denovo))
    # plot data
    cplot <- ggplot(cos_long, aes(expected_denovo, inferred_denovo)) +
      geom_tile(aes(fill = cosine_similarity)) +
      geom_text(aes(label = round(cosine_similarity, 3))) +
      scale_fill_gradient(low = "white", high = "darkgreen") +
      ggtitle("Cosine similarity matrix (expected vs. inferred)") +
      xlab("Expected") +
      ylab("Inferred")
  }

  return(list(cplot))
}

#-------------------------------------------------------------------------------

multi.plot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# ------------------------------------------------------------------------------

final.plot <- function(exp_alpha, inf_alpha, exp_denovo, inf_denovo) {

  alpha <- basilica:::plot.alpha(exp_alpha, inf_alpha)
  beta.cosine <- plot.beta.cosine(exp_denovo, inf_denovo)
  p1 <- list(ggpubr::ggarrange(plotlist=c(alpha, beta.cosine), ncol = 2))

  exp.beta <- plot.beta(exp_denovo)
  #ggpubr::annotate_figure(exp.beta, top = ggpubr::text_grob("Expected Signatures", color = "red", face = "bold", size = 14))
  inf.beta <- plot.beta(inf_denovo)
  #ggpubr::annotate_figure(inf.beta, top = ggpubr::text_grob("Inferred Signatures", color = "red", face = "bold", size = 14))
  p2 <- list(ggpubr::ggarrange(plotlist=c(exp.beta, inf.beta), ncol = 2))
  p <- c(p1, p2)

  return(basilica:::multi.plot(plotlist = p))
}


#-------------------------------------------------------------------------------

#' @import dplyr
fill_tibble <- function(x) {

  data <- tibble::tibble(
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

#-------------------------------------------------------------------------------

visData <- function(x) {

  df <- tibble::tibble(
    targetX = character(),
    inputX = character(),
    num_samples = numeric(),

    fixed_acc = numeric(),
    denovo_ratio = numeric(),
    denovo_sim = numeric(),
    mae = numeric(),
  )

  for (i in 1:nrow(x)) {

    print(paste('row:', i))

    input_cat <- x[i, ]$input_cat[[1]]

    exp_exposure <- x[i, ]$exp_exposure[[1]]
    exp_fixed <- x[i, ]$exp_fixed[[1]]
    exp_denovo <- x[i, ]$exp_denovo[[1]]

    inf_exposure <- x[i, ]$inf_exposure[[1]]
    inf_fixed <- x[i, ]$inf_fixed[[1]]
    inf_denovo <- x[i, ]$inf_denovo[[1]]

    m <- x[i, ]$x[[1]]
    mr <- reconstruction_matrix(m, inf_exposure, inf_fixed, inf_denovo)

    # fixed signatures
    if (is.null(input_cat)) {inp <- c()} else {inp <- rownames(input_cat)}
    if (is.null(exp_fixed)) {exp <- c()} else {exp <- rownames(exp_fixed)}
    if (is.null(inf_fixed)) {inf <- c()} else {inf <- rownames(inf_fixed)}
    fixed_TP <- length(intersect(inf, exp))
    fixed_FP <- length(setdiff(inf, exp))
    fixed_TN <- length(setdiff(setdiff(inp, exp), inf))
    fixed_FN <- length(setdiff(exp, inf))
    #if (length(exp)==0) {fixed_TPR <- fixed_TP+1} else {fixed_TPR <- fixed_TP / length(exp)}
    #fixed_prec <- fixed_TP / (fixed_TP + fixed_FP)
    #fixed_rec <- fixed_TP / (fixed_TP + fixed_FN)
    fixed_acc <- (fixed_TP + fixed_TN) / (fixed_TP + fixed_TN + fixed_FP + fixed_FN)

    # denovo signatures
    if (is.null(exp_denovo)) {n_exp_denovo <- 0} else {n_exp_denovo <- nrow(exp_denovo)}
    if (is.null(inf_denovo)) {n_inf_denovo <- 0} else {n_inf_denovo <- nrow(inf_denovo)}
    denovo_ratio <- (n_inf_denovo + 1) / (n_exp_denovo + 1)
    denovo_sim <- denovo.sim(exp = exp_denovo, inf = inf_denovo)

    # goodness of fitness
    #fitness <- fitness.quality(m, mr)
    mae <- compute.mae(m, mr)
    #mse <- compute.mse(m, mr)

    # fill visualization tibble
    df <- df %>% tibble::add_row(
      targetX =  x[i, ]$targetX,
      inputX = x[i, ]$inputX,
      num_samples = nrow(x[i, ]$exp_exposure[[1]]),

      fixed_acc = fixed_acc,
      denovo_ratio = denovo_ratio,
      denovo_sim = denovo_sim,
      mae = mae,
    )
  }

  #df$targetX <- factor(df$TargetX, levels = c("low", "medium", "high"))
  #df$inputX <- factor(df$InputX, levels = c("low", "medium", "high"))

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

denovo.sim <- function(exp, inf) {

  if (length(exp)==0 | length(inf)==0) {
    return(NULL)
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
      }
    }

    #------------------------------
    similarity <- 0
    for (i in 1:min(nrow(inf), nrow(exp))) {
      max = which(df == max(df), arr.ind = TRUE)
      similarity <- similarity + df[max]

      row_index <- as.numeric(max)[1]
      col_index <- as.numeric(max)[2]

      df <- df[-c(row_index), -c(col_index)]
    }
    #-------------------------------

  }

  return(similarity / nrow(inf))
}

#-------------------------------------------------------------------------------






