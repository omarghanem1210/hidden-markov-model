em.poisson.hmm <- function(delta, transition_matrix, observations, lambda, epsilon = 0.01, n_iter=0){
  n <- length(observations)
  m <- nrow(transition_matrix)

  
  for (i in 1:n_iter) {
    # E step
    #Replace the indicator variables by their conditional expectation
    u <- matrix(0, nrow = m, ncol = n)
    v <- list()
    probabilities <- compute_forward_backward_probabilities(delta, transition_matrix, observations, lambda)
    alpha <- probabilities$alpha
    beta <- probabilities$beta
    likelihood.poisson <- drop(t(alpha[, n]) %*% beta[, n])
    
    for(t in 1:n){
      u[, t] <- (alpha[, t] * beta[, t]) / likelihood.poisson
      t_matrix <- matrix(0, nrow = m, ncol = m)
      
      for (k in 1:m) {
        if (t == 1) {
          break
        }
        p_k_t <- P.poisson.hmm(observations[t], lambda)[k, k]
        t_matrix[, k] <- (alpha[, t-1] %*% diag(transition_matrix[, k]) * p_k_t * drop(beta[k, t])) / likelihood.poisson 
      }
      v[[t]] <- t_matrix
    }
    
    #M step
    delta <- t(u[, 1])
    f <- Reduce('+', v)
    for (j in 1:m) {
      transition_matrix[j, ] <- f[j, ] / sum(f[j, ])
    }
    
    lambda <- c((u[,2:n] %*% observations[2:n]) / rowSums(u[,2:n]))
    
    
  }
  
  
  return(list(delta=delta, transition_matrix=transition_matrix, lambda = lambda))
  
}



compute_forward_backward_probabilities <- function(delta, transition_matrix, observations, lambda){
  #Returns Forward and backward probabilities matrices where 
  #row represents state and column represents time
  
  len <- length(observations)
  m <- nrow(transition_matrix)
    
  alpha <- matrix(0, nrow = m, ncol = len)
  beta <- matrix(0, nrow = m, ncol = len)
  
  alpha[, 1] <- t(delta %*% P.poisson.hmm(observations[1], lambda))
  beta[, len] <- matrix(1, nrow = m, ncol = 1)
  i = 2
  
  for (observation in observations[2:len]) {
    alpha[, i] <- t(alpha[, i-1] %*% transition_matrix %*% P.poisson.hmm(observation, lambda))
    beta[, len - i + 1] <- transition_matrix %*% P.poisson.hmm(observations[len - i + 2], lambda) %*% beta[, len - i + 2]
    i = i + 1
  }
  
  return(list(beta=beta, alpha=alpha))
}

P.poisson.hmm <- function(observation, lambda){
  # Returns a diagonal matrix such that the elements at row i column i is pr(X_t = x given C_t = i)
  diagonal <- dpois(observation, lambda)
  return(diag(diagonal))
}