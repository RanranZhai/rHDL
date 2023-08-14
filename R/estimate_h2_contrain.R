

## Likelehood function for h11 and h22
llfun <- function(param, N, M, Nref = 1000, lam, bstar, lim = exp(-18)){
  h2 = param[1]
  int = param[2]
  lamh2 = h2/M*lam^2 - h2*lam/Nref + int*lam/N
  lamh2 = ifelse(lamh2<lim, lim,lamh2)
  ll = sum(log(lamh2)) + sum(bstar^2/(lamh2))
  return(ll)
}




estimate_h2 <- function(Z, N = 1e4, eigen = .9999) {
  Nref <- N
  M <- length(Z)

  bhat1 <- Z/sqrt(N+Z**2) ## used new bhat in reverse HDL method
  z11 <- (bhat1)^2

  reg <- lm(z11 ~ pp)
  h11.ols <- c(summary(reg)$coef[1:2,1:2]*c(N, M))

  ## add weight
  h11v <- (h11.ols[2]*pp/M + 1/N)^2
  reg <- lm(z11 ~ pp, weight = 1/h11v)
  h11.wls <- c(summary(reg)$coef[1:2,1:2]*c(N,M))

  v <- which(cumsum(aa$d)/sum(aa$d) >= eigen)[1]
  V <- aa$u[,1:v]
  lam <- aa$d[1:v]
  bstar1 <- crossprod(V, bhat1)

  opt <- optim(c(h11.wls[2],1), llfun, N = N, Nref = Nref, lam = lam, bstar = bstar1, M = M,
               lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(Inf,Inf), hessian = T)

  h2 <- opt$par[1]
  int <- opt$par[2]
  val <- opt$value
  conv <- opt$convergence
  hessian <- as.numeric(opt$hessian)
  se <- sqrt(1/(hessian[1]))

  return(c(h2, se, int, val, conv, hessian))
}
















