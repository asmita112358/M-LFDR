rm(list = ls())
source("~/Documents/OneDrive - Texas A&M University/MLFDR/hierarchical-EM.R")
source("~/Documents/OneDrive - Texas A&M University/MLFDR/EM_general_fun.R")
source("~/Downloads/dact.R")

#libraries
library(HDMT)
library(locfdr)
library(qvalue)
library(emdbook)

sim.size = function(tau,size = 0.05)
{
  X = rnorm(n, 1, sd = 0.75)
  Z = rnorm(n, 1.2, sd = 1)
  M = matrix(nrow = m, ncol = n)
  Y = matrix(nrow = m, ncol = n)
  #Y = vector()
  gamma = sample(1:4, m, replace = T, prob = pi)
  del = rnorm(m, 0, 0.5)
  alpha = vector()
  beta = vector()
  tn = vector()
  tp = vector()
  theta = rnorm(m, 2, 0.25)
  
  vec1 = rnorm(m, 0.2*tau, kap)
  vec2 = rnorm(m, 0.3*tau, psi)
  for(i in 1:m)
  {
    if(gamma[i] == 1){  ##h00
      alpha[i] = 0
      beta[i] = 0
      
    }else if(gamma[i] ==3){  ##h01
      alpha[i] = 0
      beta[i] = vec2[i]
      
    }else if(gamma[i] ==2){  ##h10
      alpha[i] = vec1[i]
      beta[i] = 0
      
    }else{    ##h11
      alpha[i] = vec1[i]
      beta[i] = vec2[i]
    }
    
    tn[i] = alpha[i]*beta[i] ==0
    tp[i] = alpha[i]*beta[i] !=0
    
    
    M[i,] = alpha[i]*X + rnorm(n)
    Y[i,] = beta[i]*M[i,]  + rnorm(n)  
    
    
  } 
  
  
  ##Estimate coefficients
  alpha_hat = vector()
  beta_hat = vector()
  var_alpha = c()
  var_beta = c()
  p1 = vector()
  p2 = vector()
  for(i in 1:m)
  {
    obj1 = lm(M[i,] ~ -1 + X )
    obj2 = lm(Y[i,] ~ -1 + M[i,])
    table1 = coef(summary(obj1))
    table2 = coef(summary(obj2))
    
    
    alpha_hat[i] = table1[1,1]
    beta_hat[i] = table2[1,1]
    p1[i] = table1[1,4]
    p2[i] = table2[1,4]
    var_alpha[i] = table1[1,2]^2
    var_beta[i] = table2[1,2]^2
  }
  
  x = cbind(alpha_hat, beta_hat)
  fit = EM_fun(x, k = 4, var_alpha, var_beta)
  pi = fit$lambda
  mu = fit$mu
  k = length(mu)
  sigma = fit$sigma
  lfdr = vector()
  t = matrix(nrow = m, ncol = k)
  for(i in 1:m)
  {
    for(j in 1:k)
    {
      t[i,j] = pi[j]*emdbook::dmvnorm(x[i,], mu[[j]], sigma[j,i,,])
    }
    
    lfdr[i] = (t[i,1] + t[i,2] + t[i,3])/(t[i,1] + t[i,2] + t[i,3] + t[i,4])
    
  }
  st.lfdr<-sort(lfdr)
  k=1
  
  while(k<m && ((1/k)*sum(st.lfdr[1:k])) <= size){
    k=k+1
  }
  k<-k-1
  lfdrk<-st.lfdr[k]
  reject<- lfdr<=lfdrk
  accept<- lfdr>lfdrk
  fdr = sum(reject*tn)/max(1,sum(reject))
  pow = sum(reject*tp)/sum(tp)
  return(c(fdr, pow))
}

m = 1000
n = 100
pi = c(0.6, 0.2, 0.2, 0.2)
tau = 2
kap = 0.7
psi = 1.2
tau = 1.5
sim.size(tau)

