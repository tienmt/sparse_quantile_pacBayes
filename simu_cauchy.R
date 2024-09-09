### sparse quantile
my.quantile.loss = function(u,tau) (u>0)*tau*u + (u<=0)*(1-tau)*abs(u)
setwd("~/Dropbox/ongoing_works/ON Going/sparse QUANTILE regres/Rcodes")
source('HSBQR.R')
library(rqPen); library(hqreg)
Iters = 30000
burnin = 500
tau = 1  # in the prior
# random data generation
n = 50  # samples
p = 100   # predictors
s0 = 25    # true sparsity

mytau = 0.9

lasso = mala = lmc = horSH = list()
for (ss in 1:100) {
  # generate data 
  X = matrix(rnorm(n*p),nc=p); tX= t(X)
  beta0 = rep(0,p)
  beta0[1:s0] = rnorm(s0,0,sd=1)
  Y = X%*%beta0 + rcauchy(n)
  
  # lasso 
  cv.lasso.hqreg = cv.hqreg(X, Y, FUN = 'hqreg_raw', nfolds = 5,intercept = FALSE, method = 'quantile',tau = mytau,alpha = 1)
  predict.lasso = predict(cv.lasso.hqreg, X, lambda = "lambda.min")
  lasso[[ss]] = c(mean( my.quantile.loss(Y - predict.lasso,tau = mytau)), mean((coefficients(cv.lasso.hqreg)- beta0)^2) )
  
  ### MALA
  Bm_hinge = matrix( 0 ,nrow = p)
  h = 1/(p)^2.4 # 2.4
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
      M = tam;   a = a+1
    } 
    if (s>burnin)Bm_hinge = Bm_hinge + M/(Iters-burnin)
  }
  print(a/Iters)
  
  ### LMC
  Bm_lmc = matrix( 0 ,nrow = p)
  h = 1/(p)^2.5/2 # 2.4
  M = coefficients(cv.lasso.hqreg)
  for(s in 1:Iters){
    YXm = Y-X%*%M
    M = M + h*tX%*%( YXm > 0 )*mytau + h*tX%*%( YXm <= 0 )*(mytau-1) - h*sum(4*M/(tau^2 + M^2) ) +sqrt(2*h)*rnorm(p)
    if (s>burnin)Bm_lmc = Bm_lmc + M/(Iters-burnin)
  }
  lmc[[ss]] = c(mean( my.quantile.loss(Y - X%*%Bm_lmc,tau = mytau)) , mean((Bm_lmc- beta0)^2) )
  mala[[ss]] = c(mean( my.quantile.loss(Y - X%*%Bm_hinge,tau = mytau)), mean((Bm_hinge- beta0)^2) )
  
  ### horseshoe bayes
  hsquantile = HSBQR(Y,X,quant = mytau, nsave = 1000,nburn = 500,thin = 1,iter = 1000)
  horSH[[ss]] = c(mean( my.quantile.loss(Y - X%*%hsquantile,tau = mytau)), mean((hsquantile- beta0)^2))
  
  print(ss)
}
save.image("/home/ahomef/t/thetm/Documents/simCAU_n50p100s25_tau09.rda")

