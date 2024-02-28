
rm(list = ls())
source("~/Documents/OneDrive - Texas A&M University/Documents/fair-ML/RCodes/mediation_git/semifinal_fun.R")
#source("~/Documents/OneDrive - Texas A&M University/Documents/fair-ML/RCodes/mediation_git/dact2.R")
source("~/Downloads/dact.R")
library(HDMT)
library(locfdr)
library(foreach)
library(doParallel)






comb.fcn = function(X, M, Y, Z, tp, tn, size = 0.05)
{
  m = nrow(Y)
  n = ncol(Y)
  ##Calculating pvalues
  p1 = vector()
  p2 = vector()
  t1 = vector()
  t2 = vector()
  for(i in 1:m)
  {
    obj = lm(M[i,] ~ -1 + X)
    
    obj2 = lm(Y[i,] ~ -1 + M[i,] + X)
    
    t1[i] = summary(obj)$coefficients[1,3]
    t2[i] = summary(obj2)$coefficients[1,3]
    p1[i] = summary(obj)$coefficients[1,4]
    p2[i] = summary(obj2)$coefficients[1,4]
  }
  input_pvalues = cbind(p1, p2)
  pmax = apply(input_pvalues, 1, max)
  #p_dact = DACT(p1, p2, correction = "NULL")
  p_dact = DACT(t1,t2, 100)
  ##null estimation
  
  nullprop = null_estimation(input_pvalues)
  fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                           nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
  threshhold = max(pmax[fdr_hdmt<= size])
  rej1 = pmax <= threshhold
  fdr1 = sum(rej1*tn)/max(1,sum(rej1))
  pow1 = sum(rej1*tp)/sum(tp)
  
  rej2 = p_dact <= size
  fdr2 = sum(rej2*tn)/max(1,sum(rej2))
  pow2 = sum(rej2*tp)/sum(tp)
  
  
  sample2 = est.coeff(X, Y, M, Z)
  alpha = sample2[,1]
  beta = sample2[,2]
  
  ##call a function here to estimate the value of pi_start using storey's method.
  ##pi_start = pi.init(p1, p2, lambda)
  #if(nullprop$alpha10*nullprop$alpha01 == 0)
 # {
  #  print("One of the initial values of pi is zero, resorting to default")
  #  frac = (1-nullprop$alpha00)/3
  #  pi_start = c(nullprop$alpha00, frac, frac, frac)
  #}else{
  #  pi_start = c(nullprop$alpha00, nullprop$alpha10, nullprop$alpha01, 1 - (nullprop$alpha00 + nullprop$alpha10 + nullprop$alpha01))
 #}
  pi_start = c(0.7, 0.1, 0.1, 0.1)
    
  obj = maximization(alpha, beta, X, Y, M, pi_start, maxiter = 1000)
  parm = obj$par
  obj_rej = lfdr(sample2, parm)
  rej3 = obj_rej$rej
  fdr3 = sum(rej3*tn)/max(1,sum(rej3))
  pow3 = sum(rej3*tp)/sum(tp)
  
  return(c(fdr1, fdr2, fdr3, pow1, pow2, pow3))
  
}

##sparse alternative

#tau = 0.5
#tau = 1.5
tau = 2
n = 100
m = 1000
pi = c(0.4, 0.2, 0.2, 0.2)
kap = 1
psi = 2
v = matrix(nrow = 20, ncol = 6)
#cl = makeCluster(4)
#registerDoParallel(cl)
#foreach(l = 1:20, .combine = "rbind", .packages = c("locfdr","HDMT","dplyr", "MASS", "emdbook"))%dopar%{
for(l in 1:20)
{
  X = rnorm(n, 1, sd = 0.75)
  Z = rnorm(n, 0.3, sd = 1)
  generate.obj = generate(m,n,pi, tau,kap, psi, X, Z)
  M = generate.obj$M
  Y = generate.obj$Y
  tp = generate.obj$tp
  tn = generate.obj$tn
 
 v[l,]=comb.fcn(X, M, Y,Z, tp, tn)
 print(l)
}
#stopCluster(cl)
means = colMeans(v, na.rm = T)
std = sqrt(apply(v,2, var, na.rm = T))


means_mat[3,] = means
sd_mat[3,] = std

write.csv(rbind(means_mat, sd_mat), file = "dense_alt_3.csv")
