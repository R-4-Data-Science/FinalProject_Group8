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


# Function to calculate the confusion matrix
calculate_confusion_matrix <- function(y, predictions) {
  confusion_matrix <- table(y, predictions)
  return(confusion_matrix)
}


# Function to calculate prevalence
calculate_prevalence <- function(y) {
  prevalence <- sum(y) / length(y)
  return(prevalence)
}

# Function to calculate accuracy
calculate_accuracy <- function(confusion_matrix) {
  accuracy <- (sum(diag(confusion_matrix))) / length(confusion_matrix)
  return(accuracy)
}

# Function to calculate sensitivity
calculate_sensitivity <- function(confusion_matrix) {
  sensitivity <- confusion_matrix[1, 1] / (sum(confusion_matrix[1, ]))
  return(sensitivity)
}

# Function to calculate specificity
calculate_specificity <- function(confusion_matrix) {
  specificity <- confusion_matrix[2, 2] / (sum(confusion_matrix[2, ]))
  return(specificity)
}

# Function to calculate false discovery rate
calculate_false_discovery_rate <- function(confusion_matrix) {
  false_discovery_rate <- confusion_matrix[0, 1] / (sum(confusion_matrix[0, ]))
  return(false_discovery_rate)
}

# Function to calculate diagnostic odds ratio
calculate_diagnostic_odds_ratio <- function(sensitivity, specificity) {
  diagnostic_odds_ratio <- (sensitivity / (1 - specificity)) / ((1 - sensitivity) / specificity)
  return(diagnostic_odds_ratio)
}

# Function to plot metrics over a grid of cutoff values
plot_metrics_over_cutoff_values <- function(fitted_probabilities, y, cutoff_values) {
  metrics <- matrix(nrow = length(cutoff_values), ncol = 6)
  for (i in 1:length(cutoff_values)) {
    cutoff <- cutoff_values[i]
    predictions <- ifelse(fitted_probabilities > cutoff, 1, 0)
    confusion_matrix <- calculate_confusion_matrix(y, predictions)
    metrics[i, ] <- c(
      calculate_prevalence(y),
      calculate_accuracy(confusion_matrix),
      calculate_sensitivity(confusion_matrix),
      calculate_specificity(confusion_matrix),
      calculate_false_discovery_rate(confusion_matrix),
      calculate_diagnostic_odds_ratio(calculate_sensitivity(confusion_matrix), calculate_specificity(confusion_matrix))
    )
  }
  
  plot(cutoff_values, metrics[, 2], type = "l", col = "blue", xlab = "Cutoff Value", ylab = "Accuracy")
  lines(cutoff_values, metrics[, 3], type = "l", col = "green")
  lines(cutoff_values, metrics[, 4], type = "l", col = "red")
  legend("topright", c("Accuracy", "Sensitivity", "Specificity"), col = c("blue", "green", "red"))
}


# Reference: https://bard.google.com/chat/139b518c839ed4ce