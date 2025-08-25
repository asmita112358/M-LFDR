##cox hdi code
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
source("~/Downloads/MLFDR_Plos/codes/MLFDR/EM_funs.R")
#Run from line 193
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
    z = sqrt(n)*bhat[l]/sqrt(variance[l,l])
    pval[l] = 2*min(pnorm(z), 1-pnorm(z))
    rm(z)
  }
  return(list(bhat = bhat, var = diag(variance)/n, pval = pval))
  
}


##Analysis for TCGA survival data, with different variance cutoffs

meth = readRDS("CpG0.08_cox.RDS")
surv2 = readRDS("clinical_cox.RDS")
cpg_names = colnames(meth)
surv = surv2$surv_time
cens = surv2$censoring
smoking = as.numeric(surv2$smoking)
age= surv2$age_at_diag
gc()
smoking = scale(smoking)
age = scale(age)
p = ncol(meth)

alpha <- beta <- var_alpha <- var_beta <- c()
p1 <- p2 <- c()
for(i in 1:p)
{
  obj1 = lm(meth[,i] ~ smoking + age)
  alpha[i] = obj1$coefficients["smoking"]
  var_alpha[i] = coef(summary(obj1))["smoking",2]^2
  p1[i] = coef(summary(obj1))["smoking",4]
  print(i)
}
meth = cbind(meth,smoking, age)    #########
obj2 = cox_inference(meth, surv, cens, kk = 5)
##saveRDS(obj2, file = "coxfit8.RDS")
beepr::beep(4)
p2 = obj2$pval[1:129]

size = 0.05


###MLFDR
beta = obj2$bhat[1:length(p1)]
var_beta = obj2$var[1:length(p1)]
#Save coxfit
data_cox = data.frame(alpha, beta, var_alpha, var_beta, p1, p2)
saveRDS(data_cox, "cox_fit_2.RDS")
readRDS(data_cox, "cox_fit_2.RDS")
data_cox = readRDS("cox_fit_2.RDS")
obj_mlfdr = EM_fun(cbind(data_cox$alpha, data_cox$beta), k = 4, data_cox$var_alpha, data_cox$var_beta, lambda.init = c(0.4, 0.2, 0.2, 0.2), kappa.init = 1, psi.init = 1, 
                   kappa_int = c(0.03,10), psi_int = c(0.03,10), epsilon = 0.1)


pi = obj_mlfdr$lambda
mu = obj_mlfdr$mu
k = length(mu)
sigma = obj_mlfdr$sigma
lfdr = vector()
m = length(data_cox$alpha)
t = matrix(nrow = m, ncol = k)

x = cbind(data_cox$alpha, data_cox$beta)
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
  print((1/k)*sum(st.lfdr[1:k]))
}

k<-k-1
lfdrk<-st.lfdr[k]

#HDMT
input_pvalues = na.omit(cbind(data_cox$p1, data_cox$p2))
pmax = apply(input_pvalues, 1, max)
nullprop = null_estimation(input_pvalues)
fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                         nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=1)
threshhold = max(pmax[fdr_hdmt<= size])
rej_HDMT = cpg_names[pmax <= threshhold]
rej_HDMT
cpg_data = data.frame(cpg_names, lfdr,pmax)

reject_MLFDR = cpg_data[lfdr<=lfdrk,]%>% arrange(pmax)
reject_HDMT = cpg_data[pmax<=threshhold,]
##DACT
p_dact = DACT(data_cox$p1, data_cox$p2, correction = "JC")
rej_DACT = cpg_data[p_dact <= size,]


venn_list = list(DACT = rej_DACT$cpg_names, HDMT = reject_HDMT$cpg_names,MLFDR = reject_MLFDR$cpg_names)
cox_venn = ggvenn::ggvenn(venn_list, fill_color = c("#F8766D","#00BA38","#619CFF"),
stroke_size = 0.5) + labs(title = "CpG methylation sites identified by all methods")+
  theme(plot.title = element_text(hjust = 0.5))  
cox_venn
other_two = merge(reject_HDMT, rej_DACT, by = "cpg_names", all = TRUE)
unique_to_MLFDR =  anti_join(reject_MLFDR, other_two, by = "cpg_names") %>% arrange(pmax)

ggsave("Cox_venn.png", plot = cox_venn,bg = "white", width = 7, height = 7, dpi = 1200)  # Adjust width and height as needed


print(xtable(unique_to_MLFDR[1:5,1:2], digits = c(0,4,4)), type = "latex")


##Analysis for TCGA survival data, "null scenario"
meth = readRDS("CpG0.08_cox.RDS")
surv2 = readRDS("clinical_cox.RDS")
permute = sample(1:600)
meth = meth[permute,]
cpg_names = colnames(meth)
surv = surv2$surv_time
cens = surv2$censoring
smoking = as.numeric(surv2$smoking)
age= surv2$age_at_diag
gc()
smoking = scale(smoking)
age = scale(age)
p = ncol(meth)

alpha <- beta <- var_alpha <- var_beta <- c()
p1 <- p2 <- c()
for(i in 1:p)
{
  obj1 = lm(meth[,i] ~ smoking + age)
  alpha[i] = obj1$coefficients["smoking"]
  var_alpha[i] = coef(summary(obj1))["smoking",2]^2
  p1[i] = coef(summary(obj1))["smoking",4]
  print(i)
}
meth = cbind(meth,smoking, age)    #########
obj2 = cox_inference(meth, surv, cens, kk = 5)
##saveRDS(obj2, file = "coxfit8.RDS")
beepr::beep(4)
p2 = obj2$pval[1:129]
var_beta = obj2$var[1:129]
beta = obj2$bhat[1:129]

input_pvalues = na.omit(cbind(p1,p2))
pmax = apply(input_pvalues, 1, max)
nullprop = null_estimation(input_pvalues)
fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                         nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=1)
threshhold = max(pmax[fdr_hdmt<= size])
rej_HDMT = cpg_names[pmax <= threshhold]
rej_HDMT
cpg_data = data.frame(cpg_names, lfdr,pmax)

reject_MLFDR = cpg_data[lfdr<=lfdrk,]%>% arrange(pmax)
reject_HDMT = cpg_data[pmax<=threshhold,]



##MLFDR
obj_mlfdr = EM_fun(cbind(alpha, beta), k = 4, var_alpha, var_beta, lambda.init = c(0.4, 0.2, 0.2, 0.2), kappa.init = 1, psi.init = 1, 
                   kappa_int = c(0.03,10), psi_int = c(0.03,10), epsilon = 0.1)


pi = obj_mlfdr$lambda
mu = obj_mlfdr$mu
k = length(mu)
sigma = obj_mlfdr$sigma
lfdr = vector()
m = length(data_cox$alpha)
t = matrix(nrow = m, ncol = k)

x = cbind(alpha, beta)
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
  print((1/k)*sum(st.lfdr[1:k]))
}

k<-k-1
lfdrk<-st.lfdr[k]

##DACT
p_dact = DACT(data_cox$p1, data_cox$p2, correction = "JC")
p_dact_corr = p.adjust(p_dact, "BH")
rej_DACT = cpg_data[p_dact_corr <= size,]



