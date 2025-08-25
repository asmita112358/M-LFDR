##two step EM, generalized to work for any k

library(NMOF)

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


LL.complete.h = function(kappa, lambda, mu, var_coeff, x, z)
{
  k = length(lambda)
  #if(length(lambda) != length(kappa)) print("Error inside LL.complete.h, length of lambda and kappa does not match")
  #if(kappa[1] != 0)print("Error inside LL.complete.h, first entry of kappa should be 0")
  m = length(var_coeff)
  var_mat <- t <- matrix(nrow = m, ncol = k)
  kappa = c(0, kappa)
  for(i in 1:m)
  {
    for(j in 1:k)
    {
      var_mat[i,j] = var_coeff[i] + sum((j==k)*(kappa[k]))#(j == 2)*kappa1 + (j==3)*kappa2
      t[i,j] = max(dnorm(x[i], mu[j], sqrt(var_mat[i,j])) , 9e-321)
    }
  }
  
  return(sum(z*log(t))+ sum(t(log(lambda)*t(z))))
  
}
# Define a wrapper function for gridsearch
#kappa_values is of length k-1, and the first entry of kappa is forcibly set to 0
objective_function <- function(kappa_values, lambda, mu, var_coeff, x, z) {
  return(-LL.complete.h(kappa_values, lambda, mu, var_coeff, x, z))  # Negate for maximization
}


EM_comp.h = function(coeff, var_coeff,k = 3, epsilon = 1e-02, maxit = 10000, lambda.init = c(0.7,0.2), mu.init = NULL)
{
  
  #Initialize lambda
  if(is.null(lambda.init)){
    stop(print("lambda.init must be provided"))
  }
  lambda = lambda.init ##Initial value of lambda, may need change depending on data
  if(k != length(lambda.init)) message("length of lambda.init is different from k, k is assigned as length(lambda.init)")
  k = length(lambda)
  m = length(coeff)
  
  ##Initialize mean vector
  if(is.null(mu.init)){
    probs = seq(0.1,1,length.out = k-1)
    mu.init = c(0, quantile(coeff, probs))
  }
  mu =  mu.init 
  if(k!= length(mu))message("length of mu.init is different from k")
  
  #kappa = rep(1, k-1) #Initial value of variance of priors.
  
  kappa = c(0, rep(1,k-1))
  
  var_mat <- matrix(nrow = m, ncol = k)
  
  for(i in 1:m)
  {
    for(j in 1:k)
    {
      var_mat[i,j] = var_coeff[i] + sum((j==k)*(kappa[k]))
    }
  }
  
  diff <- 2
  iter <- 0
  
  ll <- LL.data.h(lambda, mu, var_mat, coeff)
  z.update = matrix(nrow = m, ncol = k)
  z.update[,1] = rep(0,m)
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
    mu.new = c()
    mu.new[1] = 0
    for(j in 2:k)
    {
      z.update[,j] = (z[,k])/(var_coeff + kappa[k])
      mu.new[j] = sum(coeff*z.update[,j])/sum(z.update[,j])
    }
    # w = (z[,2])/(var_coeff + kappa1)
    # v = (z[,3])/(var_coeff + kappa2)
    
    
    # mu.new[2] = sum(coeff*w)/sum(w)
    # mu.new[3] = sum(coeff*v)/sum(v)
    
    
    ##Update kappa1 and kappa2
    # kappa1.new = optimize(LL.complete.h, interval = c(0.0001, 50), kappa2 = kappa2, lambda = lambda.new, mu = mu.new, var_coeff = var_coeff,x = coeff, z = z, maximum = TRUE )$maximum
    # kappa2.new = optimize(LL.complete.h, interval = c(0.0001, 50), kappa1 = kappa1.new, lambda = lambda.new, mu = mu.new, var_coeff = var_coeff,x = coeff, z = z, maximum = TRUE )$maximum
    lower_bounds <- c(rep(0.01, k - 1))  # First kappa is fixed at 0
    upper_bounds <- c(rep(10, k - 1))  # Adjust upper limit as needed
    
    grid_results <- gridSearch(fun = objective_function,
                               lower = lower_bounds,
                               upper = upper_bounds,
                               n = 20/(k-1),
                               lambda = lambda.new,
                               mu = mu.new,
                               var_coeff = var_coeff,
                               x = coeff,
                               z = z)
    kappa.new = c(0, grid_results$minlevels)
    ##Update var_mat
    var_mat.new = matrix(nrow = m, ncol = k)
    for(i in 1:m)
    {
      for(j in 1:k)
      {
        var_mat.new[i,j] = var_coeff[i] + sum((j == k)*kappa.new[k])
      }
    }
    
    ##update all parameters
    lambda <- lambda.new
    mu <- mu.new
    kappa <- kappa.new
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
  
  m = length(alpha)
  d1 = length(mu)
  d2 = length(theta)
  z = matrix(nrow = m, ncol = d1*d2)
  pi.init = rep(1,d1*d2)
  pi.new = runif(d1*d2)
  pi.new = pi.new/sum(pi.new)
  while(sum((pi.init - pi.new)^2) > 0.01)
  {
    pi.init = pi.new
    for(i in 1:m)
    {
      j = 0
      for(v in 1:d1)
      {
        for(u in 1:d2)
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

