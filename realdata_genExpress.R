### sparse quantile
setwd('/Users/thetm/Dropbox/ongoing_works/ON Going/sparse QUANTILE regres/Rcodes/')
my.quantile.loss = function(u,tau) (u>0)*tau*u + (u<=0)*(1-tau)*abs(u)
source('HSBQR.R')
library(rqPen); library(hqreg)
Iters = 30000
burnin = 500
tau = 1  # in the prior

# random data generation
library(abess)
data("trim32"); rownames(trim32) <- NULL
n = nrow(trim32)
p = 500
trim32$y <- scale(trim32$y)
trim32[,1:p] <- scale(trim32[,1:p])

mytau = 0.9

lasso = mala = lmc = horSH = list()
ac = c()
for (ss in 1:100) {
  # generate data 
  test = sample(1:n, 36)
  Y = trim32$y[-test]
  X = as.matrix(trim32[-test,1:p]) ; tX = t(X)
  Ytest =  trim32$y[test]
  Xtest = as.matrix(trim32[test,1:p] )
  
  # lasso 
  cv.lasso.hqreg = cv.hqreg(X, Y, FUN = 'hqreg_raw', nfolds = 5,intercept = FALSE, method = 'quantile',tau = mytau,alpha = 1)
  predict.lasso = predict(cv.lasso.hqreg, Xtest, lambda = "lambda.min")
  lasso[[ss]] = c(mean( my.quantile.loss(Ytest - predict.lasso,tau = mytau)) ,0 )
  
  while(sum(is.na(coefficients(cv.lasso.hqreg)))==500 ){
    cv.lasso.hqreg = cv.hqreg(X, Y, FUN = 'hqreg_raw', nfolds = 5,intercept = FALSE, method = 'quantile',tau = mytau,alpha = 1)
    predict.lasso = predict(cv.lasso.hqreg, Xtest, lambda = "lambda.min")
    lasso[[ss]] = c(mean( my.quantile.loss(Ytest - predict.lasso,tau = mytau)) ,0 )
  }
  
  ### MALA
  Bm_hinge = matrix( 0 ,nrow = p)
  h = 1/(p)^2.1 # 2.4
  a = 0  
  M = coefficients(cv.lasso.hqreg)
  for(s in 1:Iters){
    YXm = Y-X%*%M
    tam = M + h*tX%*%( YXm > 0 )*mytau + h*tX%*%( YXm <= 0 )*(mytau-1) -
      h*sum(4*M/(tau^2 + M^2) ) +sqrt(2*h)*rnorm(p)
    YXtam = Y-X%*%tam
    pro.tam = - sum(YXtam*( YXtam > 0)*mytau + YXtam*( YXtam <= 0)*(mytau-1) ) -sum(2*log(tau^2 + tam^2))
    pro.M = - sum(YXm*( YXm > 0)*mytau + YXm*( YXm <= 0)*(mytau-1) ) -sum(2*log(tau^2 + M^2))
    
    tran.m = -sum((M-tam -h*tX%*%( YXtam > 0 )*mytau + h*tX%*%( YXtam <= 0 )*(mytau-1)  -
                     h*sum(2*log(tau^2 + tam^2)) )^2)/(4*h)
    tran.tam = -sum((tam-M - h*tX%*%( YXm > 0 )*mytau + h*tX%*%( YXm <= 0 )*(mytau-1)  -
                       h*sum(2*log(tau^2 + M^2)) )^2)/(4*h)
    pro.trans = pro.tam+tran.m-pro.M-tran.tam
    if(log(runif(1)) <= pro.trans){
      M = tam;   a = a+1   } 
    if (s>burnin)Bm_hinge = Bm_hinge + M/(Iters-burnin)
  }
  print(ac[ss] <- a/Iters)
  
  ### LMC
  Bm_lmc = matrix( 0 ,nrow = p)
  h = 1/(p)^3 # 2.4
  M = coefficients(cv.lasso.hqreg)
  for(s in 1:Iters){
    YXm = Y-X%*%M
    M = M + h*tX%*%( YXm > 0 )*mytau + h*tX%*%( YXm <= 0 )*(mytau-1) - h*sum(4*M/(tau^2 + M^2) ) +sqrt(2*h)*rnorm(p)
    if (s>burnin)Bm_lmc = Bm_lmc + M/(Iters-burnin)
  }
  lmc[[ss]] = c(mean( my.quantile.loss(Ytest - Xtest%*%Bm_lmc,tau = mytau)) , 0)
  mala[[ss]] = c(mean( my.quantile.loss(Ytest - Xtest%*%Bm_hinge,tau = mytau)), 0 )
  
  ### horseshoe bayes
  hsquantile = HSBQR(Y,X,quant = mytau, nsave = 1000,nburn = 500,thin = 1,iter = 1000)
  horSH[[ss]] = c(mean( my.quantile.loss(Ytest - Xtest%*%hsquantile,tau = mytau)), 0)
  
  print(ss)
}
save.image("real_Gexpres_scale_tau09_1.rda")


