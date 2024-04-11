

solve.quadratic.optimization = function(a,
                                        b,
                                        filt_pi = 0.05,
                                        delta = 0.9,
                                        thr_exposure = 0.05,
                                        exposures = NULL,
                                        return_weights = FALSE) {
  
  df = data.frame(matrix(0, nrow(a), nrow(b)))
  rownames(df) = rownames(a)
  colnames(df) = rownames(b)
  
  cmp = nrow(a)
  b_m = t(as.matrix(b))
  Rinv = solve(chol(t(b_m) %*% b_m))
  
  res = lapply(
    1:nrow(a),
    FUN = function(i) {
      optim_res =
        solve.quadratic.optimization.aux(
          v = as.matrix(a)[i, ] %>% t(),
          Z = b_m,
          Rinv = Rinv,
          filt_pi = filt_pi,
          delta = delta,
          exposures = exposures,
          thr_exposure = thr_exposure,
          return_weights = return_weights)
      return(optim_res)
    }
  ) %>% setNames(rownames(a))
  
  return(res)
}


solve.quadratic.optimization.aux = function(v,
                                            Z,
                                            Rinv,
                                            filt_pi = 0.05,
                                            delta = 0.9,
                                            exposures = NULL,
                                            thr_exposure = 0.05,
                                            return_weights = FALSE) {
  d = v %*% Z
  b = c(1, rep(0, length(d)))
  C = cbind(rep(1, length(d)), diag(length(d)))
  pis =
    quadprog::solve.QP(
      Dmat = Rinv,
      factorized = TRUE,
      dvec = d,
      Amat = C,
      bvec = b,
      meq = 1
    )$solution
  pis[pis < filt_pi] = 0
  reconstructed_vector = Z %*% pis
  
  vec1 = matrix(v)
  rownames(vec1) = rownames(reconstructed_vector)
  
  cos_sim = lsa::cosine(as.numeric(v), as.numeric(reconstructed_vector)) %>%
    as.numeric()
  # cos_sim = cosine.vector(vec1, reconstructed_vector)
  
  if (is.null(exposures)) {
    if (!is.na(cos_sim) && cos_sim > delta) {
      if (return_weights) return(pis[pis > 0] %>% setNames(colnames(Z)[pis > 0]))
      return(colnames(Z)[pis > 0])
    }
    
  }
  
  # # exposure not null -> makes sense if we test more than one dn vs cosmic
  # } else {
  #   if (!is.na(cos_sim) && cos_sim > delta) {
  #     # if none or one linear combination
  #     if (sum(pis>0) <= 1) return(colnames(Z)[pis > 0])
  #
  #     keep.tmp = colnames(Z)[pis > 0]
  #     exposures.tmp = exposures[, keep.tmp]
  #     exposures.tmp$n_denovo = apply(exposures.tmp > thr_exposure, 1,
  #                                    function(x) length(unique(x))==1)
  #
  #     if (sum(exposures.tmp$n_denovo) >= nrow(exposures)*0.9)
  #       # return the reference that can be explained as a linear comb of denovo
  #       return(colnames(Z)[pis > 0])
  #   }
  # }
  return(NULL)
}


