##Parallel mediator model with survival outcome

library(glmnet)
library(survival)
library(parallel)
library(HDMT)
library(DACT)
library(qvalue)
library(xtable)
library(gtools)
library(lumi)
library(dplyr)
library(ggplot2)
library(tidyr)

##All functions used in the sim (upto line 150)
source("~/Downloads/MLFDR_Plos/codes/MLFDR/EM_funs.R")
cox_inference <- function(x, y, delta, kk){
  n = length(y)
  p = ncol(x)
  #if(typeof(x) == "list") x = matrix(unlist(x), nrow = n, ncol = p)
  if(typeof(x) == "list") x = as.matrix(x, nrow = n, ncol = p)
  pen1 <- cv.glmnet(x,Surv(y, delta),family="cox", nfolds = 10)
  
  
  
  # tuning parameter
  # if divided by 10, coverage probability is lower
  
  s.lambda <- pen1$lambda.min/kk
  
  betahat <- coef(pen1, s=s.lambda)
  expxbeta <- numeric(0)
  for(j in 1:n) expxbeta[j] <- as.numeric(exp(sum(x[j,]*betahat)))
  
  
  as <- rep(0,p)
  C <- diag(rep(1,p))
  #C = Matrix::sparseMatrix(i = 1:p, j = 1:p)
  T2 <- rep(1,p)
  lambda.cv <- (rep(1,p))
  thetahat <- C
  bhat <- rep(0,p)
  vhat <- bhat
  variance <- C
  cl <- matrix(0,p,2)
  ci.length <- rep(0,p)
  
  #############################
  ### function mu0, mu1, mu2
  #############################
  
  mu.all <- function(t, xx, ebeta,
                     nn = nrow(xx), pp = ncol(xx)){
    ind <- which(y >= t)
    
    mu0 <- mean(ebeta[ind])
    
    mu1 <- ebeta[ind]*xx[ind,]
    
    mu2 <- crossprod(mu1, xx[ind,])
    
    if(length(ind)> 1){
      return(list('mu0' = mu0,
                  'mu1' = colSums(mu1)/nn,
                  'mu2' = colSums(mu2)/nn))
    }else{
      mu2 = mu1*sum(xx[ind,])
      return(list('mu0' = mu0,
                  'mu1' = mu1/nn,
                  'mu2' = mu2/nn))
    }
  }
  
  dl <- rep(0,p)
  mu.all.comp = list()
  system.time(mu.all.comp <- mclapply(y, mu.all, xx = x, ebeta = expxbeta, nn = n, pp = p, mc.cores = detectCores()-2))
  #for(i in 1:n)
  #{
  # mu.all.comp[[i]] = mu.all(y[i], x, expxbeta)
  # cat(i)
  #}
  
  mu0 = vector()
  mu1 = matrix(nrow = n, ncol = p)
  mu2 = matrix(nrow = n, ncol = p)
  for(i in 1:n)
  {
    mu0[i] = unlist(mu.all.comp[[i]][1])
    mu1[i,] = unlist(mu.all.comp[[i]][2])
    mu2[i,] = unlist(mu.all.comp[[i]][3])                 
  }
  for(i in 1:n){
    dl <- dl + (x[i,]-unlist(mu.all.comp[[i]][2])/unlist(mu.all.comp[[i]][1]))*delta[i]
    cat(i)
  }
  dl <- -dl/n
  
  ddl <- matrix(0,p,p)
  for(i in 1:n){
    ddl <- ddl + (unlist(mu.all.comp[[i]][3])/unlist(mu.all.comp[[i]][1])-(unlist(mu.all.comp[[i]][2])/unlist(mu.all.comp[[i]][1]))%*%t(unlist(mu.all.comp[[i]][2])/unlist(mu.all.comp[[i]][1])))*delta[i]
    #mu0y = mu0(y[i])
    #ddl <- ddl + (mu2(y[i])/mu0y-tcrossprod(mu1(y[i])/mu0y))*delta[i]
    #rm(mu0y)
    on.exit(gc())
    
    cat(i)
  }
  ddl <- ddl/n
  
  Sigma.hat <- ddl #+ 0.01*diag(p)
  
  #e <- eigen(Sigma.hat)
  #if(n < p) nx <- t(e$vectors[,1:n]%*% diag(sqrt(e$values[1:n])))
  #t(nx) %*% diag(e$values[1:n]) %*% (nx) - Sigma.hat
  
  
  nx <- matrix(0, nrow=n*n, ncol=p)
  
  for(i in 1:n){
    for(j in 1:n){
      nx[((i-1)*n+j),] <- delta[i]*(y[j] >= y[i])*sqrt(expxbeta[j]/unlist(mu.all.comp[[i]][1]))*(x[i,]-unlist(mu.all.comp[[i]][2])/unlist(mu.all.comp[[i]][1]))
    }
  }
  
  
  mclapply(1:p, function(i){
    glmnetfit <- glmnet(nx[,-i], nx[,i], intercept=FALSE)
    #if(k==1)
    
    system.time(lambda.cv[i] <- cv.glmnet(nx[,-i], nx[,i])$lambda.min/kk)
    coeffs <- as.matrix(predict(glmnetfit,nx[,-i], type = "coefficients", s = lambda.cv[i]))[-1,]    
    C[-i,i] <- -coeffs
    T2[i] <- Sigma.hat[i,i] - t(Sigma.hat[i,-i]) %*% coeffs
    
  }, mc.cores = detectCores()-2)
  
  thetahat <- solve(diag(T2))%*%C
  bhat <- as.vector(betahat - thetahat%*%dl)
  variance <- thetahat%*%ddl%*%t(thetahat)
  pval = c()
  for(l in 1:p)
  {
    z = sqrt(n)*bhat[l]/sqrt(abs(variance[l,l]))
    pval[l] = 2*min(pnorm(z), 1-pnorm(z))
    rm(z)
  }
  return(list(bhat = bhat, var = abs(diag(variance))/n, pval = pval))
  
}

sim.size = function(tau, pi, size = 0.05)
{
  X = rnorm(n, 0, sd = 1)
  M = matrix(nrow = m, ncol = n)
  Y = vector()
  gamma = sample(1:4, m, replace = T, prob = pi)
  del = rnorm(m, 1, 0.5)
  alpha = vector()
  beta = vector()
  tn = vector()
  tp = vector()
  
  
   vec1 = rnorm(m, 2*tau, kap)
   vec2 = rnorm(m, 2.5*tau, psi)
  
  for(i in 1:m)
  {
    if(gamma[i] == 1){  ##h00
      alpha[i] = 0
      beta[i] = 0
      
    }else if(gamma[i] ==3){  ##h01
      alpha[i] = 0
      beta[i] = vec2[i]#0.7*tau
      
    }else if(gamma[i] ==2){  ##h10
      alpha[i] = vec1[i]#0.5*tau
      beta[i] = 0
      
    }else{    ##h11
      alpha[i] = vec1[i]#0.5*tau
      beta[i] = vec2[i]#0.7*tau
    }
    
    tn[i] = alpha[i]*beta[i] ==0
    tp[i] = alpha[i]*beta[i] !=0
    
    
    M[i,] = alpha[i]*X + rnorm(n)
    #Y = beta[i]*M[i,] + del[i]*X + 0.3*Z+ rnorm(n) 
    
    
  } 
  #sum(tp)/n
  #Y = colSums(beta*M) + 0.2*rnorm(n)
  u = runif(n)
  ty <- -log(u)*exp(-cbind(t(M), X)%*%c(beta, 0.3))
  
  
  ##Estimating alpha and corresponding p values
  alpha_hat = c()
  del_hat = c()
  p1 = c()
  p2 = c()
  var_alpha = c()
  var_beta = c()
  for(i in 1:m)
  {
    obj1 = lm(M[i,] ~  X)
    alpha_hat[i] = obj1$coefficients["X"]
    p1[i] = coef(summary(obj1))["X",4]
    var_alpha[i] = coef(summary(obj1))["X",2]^2
    on.exit(gc())
  }
  ##Estimating the betahat
  tcens <- rep(sort(ty)[n*censoring],n)
  delta <- (ty <= tcens)*1
  ty[ty>tcens] = tcens[1]
  obj = cox_inference(x = cbind(t(M), X), y = ty, delta, kk = 5)
  beta_hat = obj$bhat[1:m]
  var_beta = obj$var[1:m]
  p2 = obj$pval[1:m]
  ##Some p values are zero, because pnorm gives 1 after a point. Replacing those with a random uniform from 0 to min(non zero p) for computational ease. HDMT and DACT will not work otherwise.
  num.zeroes = sum(p2==0)
  if(num.zeroes>0){
    p2[p2==0] <- runif(num.zeroes,0,min(p2[p2>0]))
  }
  
  
  
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
  counter=1
  
  while(counter<m && ((1/counter)*sum(st.lfdr[1:counter])) <= size){
    counter=counter+1
  }
  counter<-counter-1
  lfdrk<-st.lfdr[counter]
  reject<- lfdr<=lfdrk
  accept<- lfdr>lfdrk
  fdr = sum(reject*tn, na.rm = TRUE)/max(1,sum(reject, na.rm = TRUE))
  pow = sum(reject*tp, na.rm = TRUE)/sum(tp, na.rm = TRUE)
  
  
  ##HDMT and DACT
  
  input_pvalues = cbind(p1, p2)
  pmax = apply(input_pvalues, 1, max)
  #p_dact = DACT(p1, p2, correction = "NULL")
  print("startdact")
  p_dact = DACT(p1,p2,correction = "JC")
  ##null estimation
  print("starthdmt")
  nullprop = null_estimation(input_pvalues)
  fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                           nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
  threshhold = max(pmax[fdr_hdmt<= size])
  rej1 = pmax <= threshhold
  fdr1 = sum(rej1*tn)/max(1,sum(rej1))
  pow1 = sum(rej1*tp)/sum(tp)
  
  rej2 = qvalue(p_dact,pi0=1)$qvalues <= size
  fdr2 = sum(rej2*tn)/max(1,sum(rej2))
  pow2 = sum(rej2*tp)/sum(tp)
  
  
  return(c(fdr, fdr1, fdr2,pow, pow1, pow2))
  
}

tau = c(0.5,1,1.5)
m = 300
n = 100
n.sim = 100
censoring = 0.7
pi = c(0.88, 0.05, 0.05, 0.02)
kap = 1
psi = 2
sim.res = matrix(nrow = 3, ncol = 12)
for (k1 in 1:3) {
  
  # Parallelize the inner loop with mclapply
  temp = matrix(nrow = n.sim,ncol = 6)
  for(iter in 1:n.sim)
  {
    temp[iter,] = sim.size(tau[k1],pi)
  }
  
  sim.res[k1, 1:6]  <- colMeans(temp,na.rm = TRUE)
  sim.res[k1, 7:12] <- apply(temp, 2, sd, na.rm = TRUE)
}


write.csv(sim.res, paste0("m",m,"n",n,"parallel_cox.csv"))
sim.res

beepr::beep(4)
