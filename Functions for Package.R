#Functions to be used in package

#ChatGPT link used by LP: https://chat.openai.com/c/99d38f73-c5da-49d0-82b8-2e8a8cbd5713




#Estimating Coefficients

logistic_regression <- function(X, y, max_iter = 100, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(cbind(1, X)) # Add intercept
  beta <- rep(0, p + 1) # Initial coefficients
  
  for (i in 1:max_iter) {
    p_hat <- 1 / (1 + exp(-X %*% beta))
    W <- diag(as.vector(p_hat * (1 - p_hat)))
    z <- X %*% beta + solve(W) %*% (y - p_hat)
    beta_new <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z
    if (sum(abs(beta_new - beta)) < tol) {
      break
    }
    beta <- beta_new
  }
  return(beta)
}


#Initial Values for Optimization
initial_values <- function(X, y) {
  X <- as.matrix(cbind(1, X)) # Add intercept
  return(solve(t(X) %*% X) %*% t(X) %*% y)
}


#Bootstrap confidence intervals
bootstrap_ci <- function(X, y, alpha = 0.05, n_bootstraps = 20) {
  bootstrapped_betas <- replicate(n_bootstraps, {
    sample_indices <- sample(nrow(X), replace = TRUE)
    logistic_regression(X[sample_indices, ], y[sample_indices])
  })
  lower <- quantile(bootstrapped_betas, alpha / 2)
  upper <- quantile(bootstrapped_betas, 1 - alpha / 2)
  return(list(lower = lower, upper = upper))
}


#Plot of the fitted logistic curve
plot_logistic_curve <- function(X, y, beta) {
  # Assuming X is one-dimensional for plotting
  seq_x <- seq(min(X), max(X), length.out = 100)
  seq_y <- 1 / (1 + exp(-(beta[1] + beta[2] * seq_x)))
  plot(seq_x, seq_y, type = 'l', xlab = 'X', ylab = 'Probability')
  points(X, y)
}

