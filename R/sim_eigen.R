## This function is to simulate different scenarios to determine best eigen cut when running
set.seed(1101)
simZ <- function(R, N) {
  #Y <- rmvnorm(N, sigma = R)
  M = nrow(R)
  h2s <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  pi0s <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  cat('Simulating Z scores')
  Z <- c()
  for (h2 in h2s) {
    zzz <- c()
    for(pi0 in pi0s) {
      pi1=pi2=(1-pi0)/2
      N_eff = round(M*(1-pi0))
      zz <- c()
      for (i in 1:100) {
        sigma = rchisq(1,1)
        b_eff = rnorm(N_eff, 0, sigma)
        b <- rep(0, times = M)
        idx_eff = sample(1:M, N_eff)
        b[idx_eff] <- b_eff
        Vg = var(Y %*% b)
        Ve = Vg/h2*(1 - h2)
        e = rnorm(N, sd = sqrt(Ve))
        x = Y %*% b + e
        z <- gwas(x, Y)
        zz <- rbind(zz, z)
        i = i+1
        #cat(i, ' ')
      }
      zzz <- rbind(zzz, zz)
    }
    Z <- rbind(Z, zzz)
    #cat(h2, ' ')
  }
  cat(' DONE\n')
  return(Z)
}



estZ <- function(Z, N, R, eigen) {
  cat('Estimating with eigen cut ', eigen)
  options(warn = -1)
  hdl <- t(apply(Z, 1, estimate_h2, N = N, eigen = eigen))
  cat(' DONE\n')
  return(hdl)
}


plot_eigen <- function(res, eigens, h2s = seq(.1,.9,.2), pi0s = seq(.1,.9,.2)) {
  dd <- data.frame(True_h2 = rep(rep(h2s, each = 100*length(pi0s)), times = length(eigens)),
                   True_pi0 = rep(rep(rep(pi0s, each = 100), times = length(h2s)), times = length(eigens)),
                   eigen = rep(eigens, each = 2500))
  res <- cbind(res,dd)
  m <- ceiling(length(eigens)/3)
  par(mfrow = c(m,3))
  for (eigen in eigens) {
    tmp <- res[res$eigen == eigen, ]
    boxplot(tmp[,1] ~ tmp$True_h2, main = paste0('eigen cut = ', eigen), xlab = 'True value', ylab = 'Estimate')
    abline(h = h2s, lty = 2, col = 2)
  }
}


#Z <- simZ(R, N)



calMSE <- function(res, h2s = seq(.1,.9,.2)) {
  nrep <- nrow(res)/length(h2s)
  if (nrep %% 1 == 0) {
    err <- 0
    for (i in 1:length(h2s)) {
      a <- (i-1)*nrep+1
      b <- i*nrep
      tmp <- sum(((res[a:b,1])-h2s[i])**2)
      err <- err+tmp
      return(err/nrow(res))
    }
  } else {
    cat('Row of data is not a multiple of object length')
  }
}



det_eigen <- function(Z, N, R, eigens = seq(.6,1,.1), plot = TRUE) {
  #Z <- simZ(R, N)
  MSE <- HDL <- c()
  for (eigen in eigens) {
    res <- estZ(Z, N, R, eigen = eigen)
    HDL <- rbind(HDL, res)
    mse <- calMSE(res)
    cat('Mean squared error is ', mse, '\n')
    MSE <- c(MSE, mse)
  }
  op <- data.frame(eigen = eigens, MSE = MSE)
  cat(paste0('***********',op$eigen[which.min(op$MSE)], ' is the most unbiased eigen cut***********\n'))
  if (plot == TRUE) {
    plot_eigen(HDL, eigens)
  }
  return(op$eigen[which.min(op$MSE)])
}

























