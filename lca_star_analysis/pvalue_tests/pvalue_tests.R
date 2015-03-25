# Multinomial supremacy test
Tn <- function(vote_vec){
  X_k = max(vote_vec)
  if (sum(vote_vec == X_k) > 1) {
    return(0)
  }
  max_ind <- match(X_k, vote_vec)
  M <- max(vote_vec[-max_ind])
  
  
  # fit liklihood ratio test from Nettleton
  log1 <- 0
  if (M > 0) {
    log1 <- M * log((2*M)/(M+X_k))
  }
  log2 <- X_k * log((2*X_k)/(M+X_k))
  res <- 2*(log1 + log2)
  print(paste("log1:", log1))
  print(paste("log2:", log2))
  print(paste("res:", res))
  return(res)
}

set.seed(12345)
tests <- rmultinom(100, 30, prob=c(0.1,0.3,0.5))
Tn_test = 1:dim(tests)[2]
unif_pval = 1:dim(tests)[2]
sup_pval = 1:dim(tests)[2]
for (i in 1:dim(tests)[2]) {
  Tn_test[i] = Tn(tests[,i])
  unif_pval[i] = 1-pchisq(Tn_test[i],1)
  sup_pval[i] = chisq.test(tests[,i], p=rep(1/length(tests[,i]), length(tests[,i])))$p.value
  print(tests[,i])
  print(dchisq(Tn_test[i],1))
  print(chisq.test(tests[,i], p=rep(1/length(tests[,i]), length(tests[,i])))$p.value)
}

library(ggplot2)
qplot(unif_pval,sup_pval, size=2, xlim=c(0,1), ylim=c(0,1)) + theme_bw()


hist(pchisq(Tn_test[Tn_test>0],1))

chisq.test(tests[,1], p=rep(1/length(tests[,1]), length(tests[,1])))
