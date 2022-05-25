
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

