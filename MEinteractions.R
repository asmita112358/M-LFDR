##interaction

rm(list = ls())

source("~/Downloads/MLFDR_Plos/codes/MLFDR/EM_funs.R")
#libraries
library(HDMT)
library(DACT)
library(locfdr)
library(qvalue)
library(emdbook)
library(parallel)
library(Rfast)




sim.size = function(tau, pi, size = 0.05)
{
  X = rbinom(n, 1, 0.2)
  Z = rnorm(n, 0, sd = 1)
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
    
    #w/o confounders
    #M[i,] = alpha[i]*X + rnorm(n)
    #Y[i,] = beta[i]*M[i,]  + rnorm(1,0.5)*X + rnorm(n)  
    
    #with confounder
    # M[i,] = alpha[i]*X + runif(1,0,0.5)*Z + rnorm(n)
    # Y[i,] = beta[i]*M[i,]  + runif(1,0,0.5)*Z + rnorm(1,0.5)*X + rnorm(n)  
    
    #ME interactions
    
    M[i,] = alpha[i]*X + rnorm(n)
    
      
  } 
  indices = sample(1:m, 20)
  theta_ind = theta[indices]/sum(theta[indices])
  M_ind = theta_ind %*% M[indices,]
  for(i in 1:m)
  {
    Y[i,] = beta[i]*M[i,] + rnorm(1,0.5)*X*M_ind + rnorm(n)
  }
  
  
  ##Estimate coefficients
  alpha_hat = vector()
  beta_hat = vector()
  var_alpha = c()
  var_beta = c()
  p1 = vector()
  p2 = vector()
  #pcs
  pc_obj = prcomp(t(M))
  pc_m = pc_obj$x[,1:2]
  for(i in 1:m)
  {
    #w/o confounders
    #obj1 = lm(M[i,] ~ -1 + X )
    #obj2 = lm(Y[i,] ~ -1 + M[i,] + X)
    
    # #with confounders
    # obj1 = lm(M[i,] ~ -1 + X)
    # obj2 = lm(Y[i,] ~ -1 + M[i,] + M[i,]:X)
    
    
    #ME interactions
    obj1 = lm(M[i,] ~  X)
    obj2 = lm(Y[i,] ~  M[i,] + X + pc_m*X)
    
    table1 = coef(summary(obj1))
    table2 = coef(summary(obj2))

    alpha_hat[i] = table1["X",1]
    beta_hat[i] = table2["M[i, ]",1]
    p1[i] = table1["X",4]
    p2[i] = table2["M[i, ]",4]
    var_alpha[i] = table1["X",2]^2
    var_beta[i] = table2["M[i, ]",2]^2
  }
  
  input_pvalues = cbind(p1, p2)
  pmax = apply(input_pvalues, 1, max)
  p_dact = DACT(p1, p2, correction = "JC")
  ##null estimation
  
  nullprop = null_estimation(input_pvalues)
  
  fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                           nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
  threshhold = max(pmax[fdr_hdmt<= size])
  rej1 = pmax <= threshhold
  fdr1 = sum(rej1*tn)/max(1,sum(rej1))
  pow1 = sum(rej1*tp)/sum(tp)
  
  rej2 = qvalue(p_dact, pi0 = 1)$qvalues <= size
  fdr2 = sum(rej2*tn)/max(1,sum(rej2))
  pow2 = sum(rej2*tp)/sum(tp)
  
  
  
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
  return(c(fdr, fdr1, fdr2, pow, pow1, pow2))
}



kap = 1
psi = 2
tau = c(0.5, 1, 1.5)
m = 10000
n = 300
n.sim = 100
pi = c(0.88, 0.05, 0.05, 0.02)
sim.res = matrix(nrow = 3, ncol = 12)
for (k in 1:3) {
  
  # Parallelize the inner loop with mclapply
  temp <- mclapply(
    1:n.sim,
    function(s) {
      sim.size(tau[k], pi)
    },
    mc.cores = detectCores() - 1  # Use available cores minus one
  ) |>
    do.call(what = rbind)
  
  sim.res[k, 1:6]  <- colMeans(temp)
  sim.res[k, 7:12] <- apply(temp, 2, sd)
}


write.csv(sim.res, paste0("m",m,"n",n,"sparse_ME_new.csv"))
sim.res

beepr::beep(4)

##Dense Alt
kap = 1
psi = 2
tau = c(0.5,1,1.5)
m = 10000
n = 300
n.sim = 100
pi = c(0.4, 0.2, 0.2, 0.2)
sim.res = matrix(nrow = 3, ncol = 12)
for (k in 1:3) {
  
  # Parallelize the inner loop with mclapply
  temp <- mclapply(
    1:n.sim,
    function(s) {
      sim.size(tau[k], pi)
      #sim.size(1,pi)
    },
    mc.cores = detectCores() - 1  # Use available cores minus one
  ) |>
    do.call(what = rbind)
  
  sim.res[k, 1:6]  <- colMeans(temp)
  sim.res[k, 7:12] <- apply(temp, 2, sd)
}
write.csv(sim.res, paste0("m",m,"n",n,"dense_ME_new.csv"))
sim.res

beepr::beep(4)

