##EM for composite alternative

##mu is a vector of length 3
##var is a matrix of length m X 3
##lambda is a vector of length 3
##x is the vector of coefficients, length m 
##please double check the dimensions before running the code, some of them are hard coded.
LL.data.h = function(lambda, mu, var, x, k = length(lambda))
{
  
  m = length(x)
  t = matrix(nrow = m, ncol = k)
  
  for(i in 1:m)
  {
    for(j in 1:k)
    {
      t[i,j] = lambda[j]*dnorm(x[i], mu[j], sqrt(var[i,j]))
    }
  }
  return(sum(log(rowSums(t))))
}

##lambda is a vector of length 3
##mu = vector of length 3
#var_coeff = variance vector of length m, obtained from regression object
#x = coefficient vector
#Q matrix defined in notes, m X 3
LL.complete.h = function(kappa1, kappa2, lambda, mu, var_coeff, x, z)
{
  k = length(lambda)
  m = length(var_coeff)
  var_mat <- t <- matrix(nrow = m, ncol = k)
  
  for(i in 1:m)
  {
    for(j in 1:k)
    {
      var_mat[i,j] = var_coeff[i] + (j == 2)*kappa1
      t[i,j] = max(dnorm(x[i], mu[j], sqrt(var_mat[i,j])) , 9e-321)
    }
  }
  
  return(sum(z*log(t))+ sum(t(log(lambda)*t(z))))
  
}

EM_comp.h = function(coeff, var_coeff,k = 3, epsilon = 1e-02, maxit = 10000)
{
  lambda = c(0.8,0.2) ##Initial value of lambda, may need change depending on data
  k = length(lambda)
  m = length(coeff)
  mu = c(0, quantile(coeff, 0.5))  ##Initial value of mean vector
  kappa1 <- 1 #Initial value of variance of priors.
  #kappa2 <- 1
  var_mat <- matrix(nrow = m, ncol = k)
  
  for(i in 1:m)
  {
    for(j in 1:k)
    {
      var_mat[i,j] = var_coeff[i] + (j == 2)*kappa1 #+ (j == 3)*kappa2
    }
  }
  
  diff <- 2
  iter <- 0
  
  ll <- LL.data.h(lambda, mu, var_mat, coeff)
  
  while(diff > epsilon & iter < maxit)
  {
    
    
    #Compute Q
    
    z = matrix(nrow = m, ncol = k)
    for (i in 1:m) {
      for (j in 1:k) {
        
        z[i,j] = lambda[j]*dnorm(coeff[i], mu[j], sqrt(var_mat[i,j]))
       
      }
    }
    z = z/rowSums(z)
    
    ##Update probabilities of each cluster
    lambda.new <- apply(z, 2, mean)
    
    ##Update mu
    w = (z[,2])/(var_coeff + kappa1)
    #v = (z[,3])/(var_coeff + kappa2)
    mu.new = c()
    mu.new[1] = 0
    mu.new[2] = sum(coeff*w)/sum(w)
    #mu.new[3] = sum(coeff*v)/sum(v)
    
    
    ##Update kappa1 and kappa2
    kappa1.new = optimize(LL.complete.h, interval = c(0.0001, 10), kappa2 = 0, lambda = lambda.new, mu = mu.new, var_coeff = var_coeff,x = coeff, z = z, maximum = TRUE )$maximum
   # kappa2.new = optimize(LL.complete, interval = c(0.0001, 10), kappa1 = kappa1.new, lambda = lambda.new, mu = mu.new, var_coeff = var_coeff,x = coeff, z = z, maximum = TRUE )$maximum

    
    ##Update var_mat
    var_mat.new = matrix(nrow = m, ncol = k)
    for(i in 1:m)
    {
      for(j in 1:k)
      {
        var_mat.new[i,j] = var_coeff[i] + (j == 2)*kappa1.new #+ (j == 3)*kappa2.new
      }
    }
    
    ##update all parameters
    lambda <- lambda.new
    mu <- mu.new
    kappa1 <- kappa1.new
    #kappa2 <- kappa2.new
    var_mat <- var_mat.new
    newobsloglik <- LL.data.h(lambda, mu, var_mat, coeff)
    #print(newobsloglik)
    #print(lambda)
    #print(kappa1)
    #print(kappa2)
    
    diff = newobsloglik - ll
    ll <- newobsloglik
    iter <- iter +1
  }
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  cat("number of iterations=", iter, "\n")
  a = list(coeff = coeff, lambda = lambda, mu = mu, var_mat = var_mat, 
           loglik = newobsloglik, posterior = z)
  #class(a) = "mixEM"
  a
}

pi.est = function(alpha, beta, mu, theta, var_mat.alpha, var_mat.beta)
{
  pi.init = rep(1,4)
  pi.new = c(0.6, 0.2, 0.2, 0.2)
  m = length(alpha)
  k = length(mu)
  z = matrix(nrow = m, ncol = k^2)
  
  while(sum((pi.init - pi.new)^2) > 0.01)
  {
    pi.init = pi.new
    for(i in 1:m)
    {
      j = 0
      for(v in 1:k)
      {
        for(u in 1:k)
        {
          j = j + 1
          z[i,j] = pi.init[j]*emdbook::dmvnorm(c(alpha[i], beta[i]), c(mu[u], theta[v]), 
                                               Sigma = matrix(c(var_mat.alpha[i,u], 0, 0, var_mat.beta[i,v]), nrow = 2)) 
          
        }
      }
    }
    #update pi
    z = z/rowSums(z)
    pi.new <- apply(z, 2, mean)
    #print(pi.new)
  }
  return(pi.new)  
}
