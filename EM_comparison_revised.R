##Simulation to compare power and time for all methods, including the 2 step EM, for d1 = 1, and d2 = 1

##Data is generated according to the scheme in Simulation 1 - Linear model.


source("~/Downloads/MLFDR_Plos/codes/MLFDR/2step_EM_revised.R")
source("~/Downloads/MLFDR_Plos/codes/MLFDR/EM_funs.R")
library(HDMT)
library(DACT)
library(locfdr)
library(qvalue)
library(emdbook)
library(parallel)


sim.size = function(tau, pi, size = 0.05)
{
  X = rbinom(n, 1, 0.2)
  #Z = rnorm(n, 0, sd = 1)
  M = matrix(nrow = m, ncol = n)
  Y = matrix(nrow = m, ncol = n)
  #Y = vector()
  gamma = sample(1:4, m, replace = T, prob = pi)
  alpha = vector()
  beta = vector()
  tn = vector()
  tp = vector()
  
  
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
    M[i,] = alpha[i]*X + rnorm(n)
    Y[i,] = beta[i]*M[i,]  + rnorm(1,0.5)*X + rnorm(n)
    
    #with confounder
    # M[i,] = alpha[i]*X + runif(1,0,0.5)*Z + rnorm(n)
    # Y[i,] = beta[i]*M[i,]  + runif(1,0,0.5)*Z + rnorm(1,0.5)*X + rnorm(n)  
    
    
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
    #w/o confounders
    obj1 = lm(M[i,] ~ -1 + X )
    obj2 = lm(Y[i,] ~ -1 + M[i,] + X)
    
    #with confounders
    # obj1 = lm(M[i,] ~ -1 + X + Z)
    # obj2 = lm(Y[i,] ~ -1 + M[i,] + Z + X)
    
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
  
  ##Standard MLFDR
  
  x = cbind(alpha_hat, beta_hat)
  fit = EM_fun(x, k = 4, var_alpha, var_beta, epsilon = 0.1)
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
  
  
  
  ##MLFDR-2
  ##Fit EM algorithm
  #x = cbind(alpha_hat, beta_hat)
  
  ##Please check inside the function and edit the starting values as necessary. Pay special care to the starting values for mu.
  #The starting values for mu should be in the order c(0, smaller mean, bigger mean). If this order is messed up you will get incorrect results.
  #Maintaining the above order ensures identifiability of the EM algorithm.
  fit_alpha = EM_comp.h(alpha_hat, var_alpha, k = d1+1, epsilon = 1e-02, maxit = 10000, lambda.init = c(0.7, 0.3))
  
  
  fit_beta = EM_comp.h(beta_hat, var_beta, k = d2+1, epsilon = 1e-02, maxit = 10000, lambda.init = c(0.7, 0.3))
  
  
  
  
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
  z = matrix(nrow = m, ncol = (1+d1)*(1+d2))
  
  for(i in 1:m)
  {
    j = 0
    for(v in 1:(1+d1))
    {
      for(u in 1:(1+d2))
      {
        j = j + 1
        z[i,j] = 100*pi[j]*emdbook::dmvnorm(c(alpha_hat[i], beta_hat[i]), c(mu[u], theta[v]), 
                                            Sigma = matrix(c(var_mat.alpha[i,u], 0, 0, var_mat.beta[i,v]), nrow = 2)) 
        
      }
    }
    temp = z[i,]*pi
    num = sum(temp[1:(1+d1+d2)])
    den = sum(temp)
    lfdr[i] = (num)/den
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
  return(c(fdr, fdr1, fdr2,fdr3, pow, pow1, pow2, pow3))
}

tau = c(0.5, 1, 1.5)
#pi = c(0.4, 0.2, 0.2, 0.2)
pi = c(0.88, 0.05, 0.05, 0.02)
m = 5000
n = 200
kap = 1
psi = 2
d1 = 1
d2 = 1
n.sim = 100
sim.res = matrix(nrow = 3, ncol = 16)
for(k1 in 1:3)
{
  temp <- mclapply(
    1:n.sim,
    function(s) {
      sim.size(tau[k1], pi, 0.05)
    },
    mc.cores = detectCores() - 1  # Use available cores minus one
  ) |>
    do.call(what = rbind)
  
  sim.res[k1, 1:8]  <- colMeans(temp, na.rm = TRUE)
  sim.res[k1, 9:16] <- apply(temp, 2, sd)
}

write.csv(sim.res, file = "EMcomparison_sparse.csv")
beepr::beep(4)


##times comparison
m = 5000
n = 200
pi = c(0.88, 0.05, 0.05, 0.02)
X = rbinom(n, 1, 0.2)
#Z = rnorm(n, 0, sd = 1)
M = matrix(nrow = m, ncol = n)
Y = matrix(nrow = m, ncol = n)
#Y = vector()
gamma = sample(1:4, m, replace = T, prob = pi)
alpha = vector()
beta = vector()
tn = vector()
tp = vector()


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
  M[i,] = alpha[i]*X + rnorm(n)
  Y[i,] = beta[i]*M[i,]  + rnorm(1,0.5)*X + rnorm(n)
  
  #with confounder
  # M[i,] = alpha[i]*X + runif(1,0,0.5)*Z + rnorm(n)
  # Y[i,] = beta[i]*M[i,]  + runif(1,0,0.5)*Z + rnorm(1,0.5)*X + rnorm(n)  
  
  
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
  #w/o confounders
  obj1 = lm(M[i,] ~ -1 + X )
  obj2 = lm(Y[i,] ~ -1 + M[i,] + X)
  
  #with confounders
  # obj1 = lm(M[i,] ~ -1 + X + Z)
  # obj2 = lm(Y[i,] ~ -1 + M[i,] + Z + X)
  
  table1 = coef(summary(obj1))
  table2 = coef(summary(obj2))
  
  
  alpha_hat[i] = table1["X",1]
  beta_hat[i] = table2["M[i, ]",1]
  p1[i] = table1["X",4]
  p2[i] = table2["M[i, ]",4]
  var_alpha[i] = table1["X",2]^2
  var_beta[i] = table2["M[i, ]",2]^2
}

p2[p2==0] = min(p2[p2>0])
start_dact = Sys.time()


p_dact = DACT(p1, p2, correction = "JC")
rej2 = qvalue(p_dact, pi0 = 1)$qvalues <= size

end_dact = Sys.time()
time_dact = end_dact-start_dact
print(time_dact)

start_hdmt = Sys.time()
input_pvalues = cbind(p1, p2)
pmax = apply(input_pvalues, 1, max)

nullprop = null_estimation(input_pvalues)
fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                         nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
threshhold = max(pmax[fdr_hdmt<= size])
rej1 = pmax <= threshhold
end_hdmt = Sys.time()

time_hdmt = end_hdmt - start_hdmt
print(time_hdmt)

start_MLFDR = Sys.time()

x = cbind(alpha_hat, beta_hat)
fit = EM_fun(x, k = 4, var_alpha, var_beta, epsilon = 0.1)
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
end_MLFDR = Sys.time()
time_MLFDR = end_MLFDR-start_MLFDR
print(time_MLFDR)

##MLFDR-2
start_MLFDR2 = Sys.time()
fit_alpha = EM_comp.h(alpha_hat, var_alpha, k = d1+1, epsilon = 1e-02, maxit = 10000, lambda.init = c(0.7, 0.3))


fit_beta = EM_comp.h(beta_hat, var_beta, k = d2+1, epsilon = 1e-02, maxit = 10000, lambda.init = c(0.7, 0.3))





mu = fit_alpha$mu
theta = fit_beta$mu
var_mat.alpha = fit_alpha$var_mat
var_mat.beta = fit_beta$var_mat
#var_mat.alpha2 = cbind(var_alpha, var_alpha + kap1, var_alpha + kap2)
#var_mat.beta2 = cbind(var_beta, var_beta + psi1, var_beta + psi2)
pi = pi.est(alpha_hat, beta_hat, mu, theta, var_mat.alpha, var_mat.beta)
lfdr = c()
z = matrix(nrow = m, ncol = (1+d1)*(1+d2))

for(i in 1:m)
{
  j = 0
  for(v in 1:(1+d1))
  {
    for(u in 1:(1+d2))
    {
      j = j + 1
      z[i,j] = 100*pi[j]*emdbook::dmvnorm(c(alpha_hat[i], beta_hat[i]), c(mu[u], theta[v]), 
                                          Sigma = matrix(c(var_mat.alpha[i,u], 0, 0, var_mat.beta[i,v]), nrow = 2)) 
      
    }
  }
  temp = z[i,]*pi
  num = sum(temp[1:(1+d1+d2)])
  den = sum(temp)
  lfdr[i] = (num)/den
}
st.lfdr<-sort(lfdr)
k=1

while(k<m && ((1/k)*sum(st.lfdr[1:k])) <= size){
  k=k+1
}
k<-k-1
lfdrk<-st.lfdr[k]
reject<- lfdr<=lfdrk

end_MLFDR2 = Sys.time()
time_MLFDR2 = end_MLFDR2-start_MLFDR2

print(time_MLFDR2)
