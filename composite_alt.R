##Simulation for composite alternative

source("~/Downloads/MLFDR_Plos/codes/MLFDR/twostep_EM.R")

library(HDMT)
library(DACT)
library(locfdr)
library(qvalue)
library(emdbook)
library(parallel)

k = 3
m = 1000
n = 100
p = c(0.7, 0.2, 0.1)
q = c(0.8, 0.1, 0.1)
tau = 1
kap1 = 1
kap2 = 2
psi1 = 1.5
psi2 = 2
pi = c()
pi[1] = p[1]*q[1]
pi[2] = (p[2] + p[3])*q[1]
pi[3] = (q[2] + q[3])*p[1]
pi[4] = 1 - pi[1] - pi[2] - pi[3]

sim.size = function(tau, size = 0.05)
{
  X = rbinom(n, 1, 0.2)
  Z = rnorm(n, 0, sd = 1)
  M = matrix(nrow = m, ncol = n)
  Y = matrix(nrow = m, ncol = n)
  gamma = sample(1:3, m, replace = T, prob = p)
  delta = sample(1:3, m, replace = T, prob = q)
  
  alpha <- beta <- tn <- tp <- c()
  
  g1 = rnorm(m, 0.2*tau, kap1)
  g2 = rnorm(m, 1.1*tau, kap2)
  
  h1 = rnorm(m, -1.2*tau, psi1)
  h2 = rnorm(m, 0.3*tau, psi2)
  alpha = 0 + (gamma == 2)*g1 + (gamma == 3)*g2
  beta = 0 + (delta == 2)*h1 + (delta == 3)*h2
  tn = alpha*beta == 0
  tp = alpha*beta != 0
  
  for(i in 1:m)
  {
    M[i,] = alpha[i]*X + rnorm(n)
    Y[i,] = beta[i]*M[i,] + rnorm(n) 
  }
  
  
  #Estimate coefficients
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
  
  
  ##Fit EM algorithm
  x = cbind(alpha_hat, beta_hat)
  
  ##Please check inside the function and edit the starting values as necessary. Pay special care to the starting values for mu.
  #The starting values for mu should be in the order c(0, smaller mean, bigger mean). If this order is messed up you will get incorrect results.
  #Maintaining the above order ensures identifiability of the EM algorithm.
  fit_alpha = EM_comp.h(alpha_hat, var_alpha, k = 3, epsilon = 1e-02, maxit = 10000)
  
  
  fit_beta = EM_comp.h(beta_hat, var_beta, k = 3, epsilon = 1e-02, maxit = 10000)
  
  
  
  
  p = fit_alpha$lambda
  q = fit_beta$lambda
  mu = fit_alpha$mu
  theta = fit_beta$mu
  var_mat.alpha = fit_alpha$var_mat
  var_mat.beta = fit_beta$var_mat
  #var_mat.alpha2 = cbind(var_alpha, var_alpha + kap1, var_alpha + kap2)
  #var_mat.beta2 = cbind(var_beta, var_beta + psi1, var_beta + psi2)
  pi = pi.est(alpha_hat, beta_hat, mu, theta, var_mat.alpha, var_mat.beta)
  lfdr = c()
  z = matrix(nrow = m, ncol = k^2)
  
  for(i in 1:m)
  {
    j = 0
    for(v in 1:k)
    {
      for(u in 1:k)
      {
        j = j + 1
        z[i,j] = 100*pi[j]*emdbook::dmvnorm(c(alpha[i], beta[i]), c(mu[u], theta[v]), 
                                            Sigma = matrix(c(var_mat.alpha[i,u], 0, 0, var_mat.beta[i,v]), nrow = 2)) 
        
      }
    }
    temp = z[i,]*pi
    num = sum(temp[1:4])+ temp[7]
    den = sum(temp)
    lfdr[i] = num/den
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
  fdr3 = sum(reject*tn)/max(1,sum(reject))
  pow3 = sum(reject*tp)/sum(tp)
  ##HDMT and DACT
  input_pvalues = cbind(p1, p2)
  pmax = apply(input_pvalues, 1, max)
  #p_dact = DACT(p1, p2, correction = "NULL")
  print("startdact")
  p_dact = DACT(p1,p2, correction = "JC")
  ##null estimation
  print("starthdmt")
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
  return(c(fdr3, fdr1, fdr2, pow3, pow1, pow2))
}

means_mat = matrix(nrow = 3, ncol = 6)
sd_mat = matrix(nrow = 3, ncol = 6)


m = 10000
n = 300
n.sim = 100
tau = c(0.5, 1, 1.5)
sim.res = matrix(nrow = 3, ncol = 12)
for (k1 in 1:3) {
  # Parallelize the inner loop with mclapply
  temp <- mclapply(
    1:n.sim,
    function(s) {
      sim.size(tau[k1])
    },
    mc.cores = detectCores() - 1  # Use available cores minus one
  ) |>
    do.call(what = rbind)
  
  sim.res[k1, 1:6]  <- colMeans(temp, na.rm = TRUE)
  sim.res[k1, 7:12] <- apply(temp, 2, sd)
}
write.csv(sim.res, paste0("m",m,"n",n,"composite_alt.csv"))
beepr::beep(4)
sim.res






