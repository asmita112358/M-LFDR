
LL.data = function(coeff_mat, mu, sigma, lambda)
{
  k = length(mu)
  m = nrow(coeff_mat)
  t = matrix(nrow = m, ncol = k)
  
for(i in 1:m)
  {
  for(j in 1:k)
  {
    t[i,j] = lambda[j]*emdbook::dmvnorm(coeff_mat[i,], mu[[j]], sigma[j,i,,])
  }
} 
  return(sum(log(rowSums(t))))
}


LL.complete = function(kappa, psi, var_alpha, var_beta, coeff_mat, mu.new, lambda.new, z)
{
  k = 4
  m = length(var_alpha)
  sigma <- array(0,dim = c(k, m, 2,2))
  t = matrix(nrow = m, ncol = 4)
  for(i in 1:4)
  {
    for(j in 1:m)
    {
      sigma[i,j,1,1] = var_alpha[j] + ifelse(i == 2||i == 4, 1, 0)*kappa
      sigma[i,j,2,2] = var_beta[j] + ifelse(i == 3||i == 4, 1, 0)*psi
      t[j,i] = max(emdbook::dmvnorm(coeff_mat[j,], mu.new[[i]], sigma[i,j,,]), 9e-321)
    }
    
  }
  
  return(sum(z*log(t))+ sum(t(log(lambda.new)*t(z))))
  
}


EM_fun <- function(coeff_mat, k = 4, var_alpha, var_beta ,epsilon = 1e-02, maxit = 10000) 
{
  lambda = c(0.7,0.1,0.1,0.1)
  coeff_mat <- as.matrix(coeff_mat)
  m <- nrow(coeff_mat)
  p <- ncol(coeff_mat)
  #tmp <- mvnormalmix.init(x = x, lambda = lambda, mu = mu, 
  #sigma = sigma, k = k, arbmean = arbmean, arbvar = arbvar)
  #lambda <- c(0.7, 0.1, 0.1, 0.1)
  mu.init = quantile(coeff_mat[,1], 0.99)
  theta.init = quantile(coeff_mat[,2], 0.99)
  psi = 1
  kappa= 1
  
  mu = list(c(0,0), c(mu.init, 0), c(0, theta.init), c(mu.init, theta.init))
  sigma <- array(0,dim = c(k, m, 2,2))
  
  for(j in 1:m)
  {
    for(i in 1:k)
    {
    
      sigma[i,j,1,1] = var_alpha[j] + ifelse(i == 2||i == 4, 1, 0)*kappa
      sigma[i,j,2,2] = var_beta[j] + ifelse(i == 3||i == 4, 1, 0)*psi
    }
    
  }
  
  
  
  diff <- 1
  iter <- 0
  
  
  ll <- LL.data(coeff_mat, mu, sigma, lambda)
  restarts <- 0
  while (diff > epsilon & iter < maxit) {
    
    ##Compute Q
    z = matrix(nrow = m, ncol = k)
    for (i in 1:m) {
      for (j in 1:k) {
        #print(j)
        z[i,j] = lambda[j]*max(emdbook::dmvnorm(coeff_mat[i,], mu[[j]], sigma[j,i,,]), 9e-321)
        #z[i, j] = 1/sum(z.denom)
      }
    }
    z = z/rowSums(z)
    #sing <- sum(is.nan(z))
    lambda.new <- apply(z, 2, mean)
    w = (z[,2] + z[,4])/(var_alpha + kappa)
    v = (z[,3] + z[,4])/(var_beta + psi)
    
    #print(sum(w))
    #print(sum(v))
    m.new = sum(coeff_mat[,1]*w)/sum(w)
    theta.new = sum(coeff_mat[,2]*v)/sum(v)
    #mu.new <- list()
    ##Check this
    #print("ran untill here")
    mu.new <- list(c(0,0), c(m.new, 0), c(0, theta.new), c(m.new, theta.new))
    
    ##update kappa and psi
    kappa.new = optimize(LL.complete, interval = c(0.1,50) , psi = psi, var_alpha = var_alpha, var_beta = var_beta, coeff_mat = coeff_mat, mu.new = mu.new,lambda.new = lambda.new, z = z, maximum = TRUE)$maximum
    psi.new = optimize(LL.complete, interval = c(0.1,50) , kappa = kappa.new, var_alpha = var_alpha, var_beta = var_beta, coeff_mat = coeff_mat, mu.new = mu.new,lambda.new = lambda.new, z = z, maximum = TRUE)$maximum
    ##Update sigma
    
    sigma.new <- array(0,dim = c(k, m, 2,2))
    for(i in 1:k)
    {
      for(j in 1:m)
      {
        sigma.new[i,j,1,1] = var_alpha[j] + ifelse(i == 2||i == 4, 1, 0)*kappa.new
        sigma.new[i,j,2,2] = var_beta[j] + ifelse(i == 3||i == 4, 1, 0)*psi.new
      }
      
    }  
    
    
    #lapply(1:k, function(j) matrix(apply(sapply(1:n, function(i) z[i, j] * (x[i, ] - mu.new[[j]]) %*%t(x[i, ] - mu.new[[j]])), 1, sum), p, p)/sum(z[,j]))
    
    
    lambda <- lambda.new
    mu <- mu.new
    kappa <- kappa.new
    psi <- psi.new
    sigma <- sigma.new
    #print(mu)
    #print(lambda)
    print(kappa)
    print(psi)
    
    newobsloglik <- LL.data(coeff_mat, mu, sigma, lambda)
    print(newobsloglik)
    
    
    diff = newobsloglik - ll
    ll <- newobsloglik
    iter <- iter +1
    
    
  }
  
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  #colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
  cat("number of iterations=", iter, "\n")
  a = list(coeff_mat = coeff_mat, lambda = lambda, mu = mu, sigma = sigma, 
           loglik = newobsloglik, posterior = z, all.loglik = ll, restarts = restarts)
  #class(a) = "mixEM"
  a
}
