## This function is to simulate NULL (where SNP is contributing to no traits)

require(mvtnorm)
set.seed(1101)
sim_null <- function(R, N, times = 10000) {
  #Y <- rmvnorm(N, sigma = R)
  Z <- matrix(NA, ncol = ncol(R), nrow = times)
  for (i in 1:times) {
    x <- rnorm(N, sd = 1)
    z <- gwas(x, Y)
    Z[i,] <- z
  }
  return(Z)
}






