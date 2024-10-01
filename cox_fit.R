#tcga_server

##defining libraries
source("~/Downloads/dact.R")
source("~/Documents/OneDrive - Texas A&M University/MLFDR/EM_general_fun.R")
library(readr)
#library(tidyverse)
#library(dplyr)
#library(janitor)
library(qvalue)
#library(forecast)
library(compositions)
#library(dplyr)
library(survival)
library(parallel)
library(doParallel)
library(MASS)
library(penalized)
library(Matrix)
library(glmnet)
library(foreach)
library(survminer)
#define functions

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

cox_inference <- function(x, y, delta, kk){
  n = length(y)
  p = ncol(x)
  #if(typeof(x) == "list") x = matrix(unlist(x), nrow = n, ncol = p)
  if(typeof(x) == "list") x = as.matrix(x, nrow = n, ncol = p)
  gc()
  system.time(pen1 <- cv.glmnet(x,Surv(y, delta),family="cox", nfolds = 10))
  
  gc()
  
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
  #mu0 <- function(t){
  #    return(mean(expxbeta*(y>=t)))
  #}
  
  #mu1 <- function(t){
  #   temp <- rep(0,p)
  #   for(j in 1:n){
  #       temp <- temp + (y[j]>=t)*expxbeta[j]*x[j,]
  #   }
  #   return(temp/n)
  #}
  
  #mu2 <- function(t){
  #    temp <- matrix(0,p,p)
  #    for(j in 1:n){
  #       temp <- temp + (y[j]>=t)*expxbeta[j]*(x[j,]%*%t(x[j,]))
  #   }
  #   return(temp/n)
  # }
  ############################
  
  
  dl <- rep(0,p)
  mu.all.comp = list()
  system.time(mu.all.comp <- mclapply(y, mu.all, xx = x, ebeta = expxbeta, nn = n, pp = p, mc.cores = 6))
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
  
  #registerDoParallel(6)
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
    z = sqrt(n)*bhat[l]/sqrt(variance[l,l])
    pval[l] = 2*min(pnorm(z), 1-pnorm(z))
    rm(z)
  }
  return(list(bhat = bhat, var = diag(variance)/n, pval = pval))
  
}
#Read data
data_sub = read_csv("tcga_cox_data_subset2.csv")
n = nrow(data_sub)

meth = data_sub[,-(1:8)]. ####################
surv = data_sub$surv_time
cens = data_sub$censoring
smoking = as.numeric(data_sub$smoking)
confounder = data_sub$age_at_diag ############
#meth = cbind(meth,smoking)
rm(data_sub)
gc()
#mod.lung <- survfit(Surv(surv, cens) ~ 1)
#ggsurvplot(mod.lung, fun = function(y) -log(y))
#meth = scale(meth)
smoking = scale(smoking)
confounder = scale(confounder) #######
p = ncol(meth)
alpha_hat <- beta_hat <- var_alpha <- var_beta <- c()
p1 <- p2 <- c()
for(i in 1:p)
{
  obj1 = lm(scale(meth[,i])~ -1 + smoking + confounder)  ############
  alpha_hat[i] = obj1$coefficients
  var_alpha[i] = coef(summary(obj1))[1,2]^2
  p1[i] = coef(summary(obj1))[1,4]
  #print(i)
  #on.exit(gc())
}
meth = cbind(meth,smoking, confounder)    #########
obj2 = cox_inference(meth, surv, cens, kk = 5)
beta_hat = obj$bhat
var_beta = obj$var
p2 = obj$pval
output = cbind(alpha_hat, beta_hat, var_alpha, var_beta, p1, p2)
write.csv(output, "coef.csv")

