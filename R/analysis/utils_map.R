

# input: data.frame of cosine similarity matrix
# output: return list of:
#   - maximum value of cosine similarity
#   - respective signatures
#   - cosine similarity matrix where its max value is set to zero
maxfinder <- function(x) {
  
  # check if data.frame has both dimnesion > 0
  if ( (nrow(x)==0) | (ncol(x)==0) ) {
    return(list(matrix=NULL, first=NULL, second=NULL, value=NULL, end=FALSE))
  }
  
  # check if all values of data.frame are zero
  if (all(x == 0)) {
    return(list(matrix=NULL, first=NULL, second=NULL, value=NULL, end=TRUE))
  }
  
  if ( (nrow(x)==1) & (ncol(x)==1) ) {
    a <- rownames(x)
    b <- colnames(x)
    c <- x[1,1]
    x[1,1] <- 0
    return(list(matrix=x, first=a, second=b, value=c, end=FALSE))
  }
  else {
    row <- which(x == max(x), arr.ind = T)[1]
    col <- which(x == max(x), arr.ind = T)[2]
    a <- rownames(x[row, ])
    b <- colnames(x[col])
    c <- max(x)
    x[row, ] <- 0
    x[, col] <- 0
    return(list(matrix=x, first=a, second=b, value=c, end=FALSE))
  }
}


#-------------------------------------------------------------------------------
mapper <- function(x) {
  ind <- min(nrow(x), ncol(x))
  
  first <- c()
  second <- c()
  values <- c()
  
  for (i in 1:(ind-1)) {
    
    res <- maxfinder(x)
    
    if ( is.null(res[[1]]) ) {
      print("invalid input!")
      return(0)
    }
    else if ( res[[5]] ) {
      break
    }
    else {
      x <- res[[1]]
      first[i] <- res[[2]]
      second[i] <- res[[3]]
      values[i] <- res[[4]]
    }
  }
  
  df <- data.frame(first, second, values)
  return(df)
}

#-------------------------------------------------------------------------------
# basilica_exposure # basilica brca exposure
# serena_exposure   # serena brca exposure


# reference --> reference signature (e.g., COSMIC)
# fixed ------> basilica fixed signatures (wide)
# denovo -----> basilica denovo signatures (wide)
# common -----> serena common signatures (wide)
# rare -------> serena rare signatures (wide)
# threshold --> real number 
## (if cosine similarity between a basilica signature and serena signature is
## higher than threshold, we map them to eachother)
map.data <- function(
    reference, 
    fixed, 
    denovo, 
    common, 
    rare, 
    threshold
) {
  
  cmatrix1 <- basilica:::cosine.matrix(denovo, common)
  denovo_common <- mapper(cmatrix1)
  denovo_common <- subset(denovo_common, values >= threshold )
  
  cmatrix2 <- basilica:::cosine.matrix(denovo, rare)
  denovo_rare <- mapper(cmatrix2)
  denovo_rare <- subset(denovo_rare, values >= threshold )
  
  cmatrix3 <- basilica:::cosine.matrix(denovo, reference)
  denovo_reference <- mapper(cmatrix3)
  denovo_reference <- subset(denovo_reference, values >= threshold )
  
  cmatrix4 <- basilica:::cosine.matrix(fixed, common)
  fixed_common <- mapper(cmatrix4)
  fixed_common <- subset(fixed_common, values >= threshold )
  
  cmatrix5 <- basilica:::cosine.matrix(fixed, rare)
  fixed_rare <- mapper(cmatrix5)
  fixed_rare <- subset(fixed_rare, values >= threshold )
  
  cmatrix6 <- basilica:::cosine.matrix(fixed, reference)
  fixed_reference <- mapper(cmatrix6)
  fixed_reference <- subset(fixed_reference, values >= threshold )
  
  cmatrix7 <- basilica:::cosine.matrix(common, reference)
  common_reference <- mapper(cmatrix7)
  common_reference <- subset(common_reference, values >= threshold )
  
  cmatrix8 <- basilica:::cosine.matrix(rare, reference)
  rare_reference <- mapper(cmatrix8)
  rare_reference <- subset(rare_reference, values >= threshold )
  
  cmatrix9 <- basilica:::cosine.matrix(rbind(fixed, denovo), rbind(common, rare))
  basilica_serena <- mapper(cmatrix9)
  basilica_serena <- subset(basilica_serena, values >= threshold )
  colnames(basilica_serena) <- c("Basilica", "Serena", "value")
  
  # aggregated fixed and denovo signatures (basilica)
  all_basilica <- rbind(fixed, denovo)
  # list of signatures name detected by basilica but un-explained by serena
  ne_basilica <- rownames(all_basilica[!(rownames(all_basilica) %in% basilica_serena$Basilica), ])
  
  # aggregated common and rare signatures (serena)
  all_serena <- rbind(common, rare)
  # list of signatures name detected by serena but un-explained by basilica
  ne_serena <- rownames(all_serena[!(rownames(all_serena) %in% basilica_serena$Serena), ])
  
  # -------------------------------- [ OUTPUT ] --------------------------------
  
  data <- list(
    denovo_common = denovo_common, 
    denovo_rare = denovo_rare, 
    denovo_reference = denovo_reference, 
    fixed_common = fixed_common, 
    fixed_rare = fixed_rare, 
    fixed_reference = fixed_reference, 
    common_reference = common_reference, 
    rare_reference = rare_reference, 
    basilica_serena = basilica_serena, 
    basilica_singles = ne_basilica, 
    serena_singles = ne_serena
  )
  return(data)
}



