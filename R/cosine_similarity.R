
cosine_sim <- function(a, b) {
  numerator <- sum(a * b)
  denominator <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  return(numerator / denominator)
}
