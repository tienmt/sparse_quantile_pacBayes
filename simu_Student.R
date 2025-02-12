### sparse quantile
my.quantile.loss = function(u,tau) (u>0)*tau*u + (u<=0)*(1-tau)*abs(u)
setwd("~/Dropbox/Apps/Overleaf/(stat n Computing) sparse QUANTILE regres/Rcodes")
source('HSBQR.R')
library(hqreg)
Iters = 30000
burnin = 500
tau = .001  # in the prior
# random data generation
n = 100  # samples
n_test = n*.2
p = 200   # predictors
s0 = 5    # true sparsity

mytau = tau_quantile = 0.9
rho.x = 0.5 ; S = matrix(rho.x, ncol = p, nrow = p);diag(S) = 1;
for (i in 1:(p-1))for (j in (i+1):p)S[i,j] = rho.x^{abs(i-j)}
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)

lasso = mala = lmc = horSH = enet= list()
for (ss in 1:100) {
  # generate data 
  beta0 = rep(0,p)
  beta0[1:s0] = rnorm(s0,0,sd=1)
  xall = matrix(rnorm((n + n_test)*p),nc=p)
  ###  (1 + xall[,1])*
  yall = xall%*%beta0 + rnorm(n +n_test ,sd=3)   # rnorm(n +n_test ,sd=3) # rcauchy(n+n_test) # 2*rt(n+n_test,df=3)
  xtest = xall[(n+1):(n+ n_test),]
  ytest = yall[(n+1):(n+ n_test),]
  X = xall[1:n,] ; tX= t(X)
  Y = yall[1:n,]
  
  # lasso 
  cv.lasso.hqreg = cv.hqreg(X, Y, FUN = 'hqreg_raw', nfolds = 5,intercept = FALSE, method = 'quantile',tau = mytau,alpha = 1)
  predict.lasso = predict(cv.lasso.hqreg, X, lambda = "lambda.min")
  lasso[[ss]] = c(mean( my.quantile.loss(ytest - predict.lasso,tau = tau_quantile)), mean((coefficients(cv.lasso.hqreg)- beta0)^2) )
  # enet
  cv.enet = cv.hqreg(X, Y, FUN = 'hqreg_raw', nfolds = 5,intercept = FALSE, method = 'quantile',tau = mytau,alpha = .5)
  predict.enet = predict(cv.enet, X, lambda = "lambda.min")
  enet[[ss]] = c(mean( my.quantile.loss(ytest - predict.enet,tau = tau_quantile)), mean((coefficients(cv.enet)- beta0)^2) )
  ### horseshoe bayes
  hsquantile = HSBQR(Y,X,quant = mytau, nsave = 1000,nburn = 500,thin = 1,iter = 1000)
  horSH[[ss]] = c(mean( my.quantile.loss(ytest - xtest%*%hsquantile,tau = tau_quantile)), mean((hsquantile- beta0)^2))
  
  ### MALA
  Bm_hinge = matrix( 0 ,nrow = p)
  h = 1/(p)^4 # 2.4
  a = 0  
  M = hsquantile
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
  h = 1/(p)^4/2 # 2.4
  M = hsquantile
  for(s in 1:Iters){
    YXm = Y-X%*%M
    M = M + h*tX%*%( YXm > 0 )*mytau + h*tX%*%( YXm <= 0 )*(mytau-1) - h*sum(4*M/(tau^2 + M^2) ) +sqrt(2*h)*rnorm(p)
    if (s>burnin)Bm_lmc = Bm_lmc + M/(Iters-burnin)
  }
  lmc[[ss]] = c(mean( my.quantile.loss(ytest - xtest%*%Bm_lmc,tau = tau_quantile)) , mean((Bm_lmc- beta0)^2) )
  mala[[ss]] = c(mean( my.quantile.loss(ytest - xtest%*%Bm_hinge,tau = tau_quantile)), mean((Bm_hinge- beta0)^2) )
  
  print(ss)
}
save.image("simuOUT/simST_n100p200s5_tau09_rX.rda")



