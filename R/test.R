## this script is for test run


R0 <- readRDS('~/Documents/Projects/HDL-TGCA/deCODE4907_lowMAF.R.rds')

idx <- sample(1:4000, 200)
R <- R0[idx,idx]


det_eigen(35559, R)


best_eigen <- det_eigen(35559, R)


R <- readRDS('~/Documents/Projects/Likelihood-TGCA/Gudjonsson1416_lowMAF.R.rds')

aa <- svd(R)
pp <- c()
for (i in 1:ncol(R)) {
  tmp <- sum((R[i,])^2)
  pp <- c(pp, tmp)
}


Z <- simZ(R, 5368)

z1 <- Z[2500,]
b1 <- z1/sqrt(5368+z1^2)


par(mfrow = c(2,2))
hist(z1)
hist(zscore(z1))
hist(b1)
hist(zscore(b1))

res <- t(apply(Z,1,estimate_h2, N=5368, eigen = .75))

res1 <- t(apply(Z,1,estimate_norm_h2, N=5368, eigen = .75))



test <- readRDS('~/Documents/Projects/Likelihood-TGCA/new_bhat/deCODE1416_h2_Pval_res.rds')

test$Chr <- unlist(strsplit(test$chr, 'chr'))[seq(2, 2*nrow(test),2)]
test$Chr[test$Chr == 'X'] <- 23
test$Chr <- as.numeric(test$Chr)

test$pos <-as.numeric(test$pos)


test$logP <- -log10(test$Pval)

runrun(test, 'Chr', 'pos', 'h2', 'logP', file = '~/test_mr.png')




gud_tscp <- readRDS('~/Documents/Projects/Likelihood-TGCA/Reactome/Gudjonsson_Reactome_Transcription_res_h2_HDL_no_constrain.rds')


runrun(gud_tscp, 'Chr', 'pos', 'h2', 'logP', file = '~/test_mr.png')





