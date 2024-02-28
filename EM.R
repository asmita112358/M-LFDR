##EM final (hopefully?)

LL.data = function(lambda, mu, x, var)
{
  k = length(lambda)
  m = length(x)
  t = matrix(nrow = m, ncol = k)
  
  for(i in 1:m)
  {
    for(j in 1:k)
    {
      t[i,j] = lambda[j]*dnorm(x[i], mu[j], sqrt(var[i]))
    }
  }
  return(sum(log(max(rowSums(t), 9e-320))))
}

##first entry of mu is always 0
EM = function(x,var, k, lambda, mu, eps = 1e-06, maxit = 1000)
{
  m = length(x)
  diff <- 2
  iter <- 0
  
  ll <- LL.data(lambda, mu, x, var)
  while(diff > epsilon & iter < maxit)
  {
    z = matrix(nrow = m, ncol = k)
    for (i in 1:m) {
      for (j in 1:k) {
        
        z[i,j] = lambda[j]*dnorm(x[i], mu[j], sqrt(var[i]))
        
      }
    }
    z = z/rowSums(z) 
    ##Update probabilities of each cluster
    lambda.new <- apply(z, 2, mean, na.rm = T)
    
    
    ##Update mu
    mu.new = c()
    mu.new = colSums(z*x, na.rm = T)/colSums(z, na.rm = T)
    mu.new[1] = 0
    
    #Update all parameters
    lambda <- lambda.new
    mu <- mu.new
    
    newobsloglik <- LL.data(lambda, mu, x, var)
    print(newobsloglik)
    print(lambda)
    print(mu)
    diff = newobsloglik - ll
    ll <- newobsloglik
    iter <- iter +1
  } 
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  cat("number of iterations=", iter, "\n")
  a = list(coeff = x, lambda = lambda, mu = mu)
  #class(a) = "mixEM"
  a
  
}