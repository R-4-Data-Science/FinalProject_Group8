#Functions to be used in package

#ChatGPT link used by LP: https://chat.openai.com/c/99d38f73-c5da-49d0-82b8-2e8a8cbd5713
#Reference Used by Felix: https://bard.google.com/chat/139b518c839ed4ce

#Estimating Coefficients

logistic_regression <- function(X, y, max_iter = 100, tol = 1e-6) {
  if (!is.matrix(X) || length(y) != nrow(X)) {
    stop("Invalid inputs for X or y")
  }
  
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(cbind(1, X))  # Add intercept
  beta <- rep(0, p + 1)        # Initial coefficients
  log_likelihood_prev <- -Inf
  
  for (i in 1:max_iter) {
    p_hat <- 1 / (1 + exp(-X %*% beta))
    W <- diag(as.vector(p_hat * (1 - p_hat)))
    z <- X %*% beta + (y - p_hat) / (p_hat * (1 - p_hat))
    
    # Handling potential singularity in the Hessian matrix
    tryCatch({
      Hessian <- t(X) %*% W %*% X
      if (det(Hessian) == 0) {
        stop("Matrix inversion not possible due to singularity.")
      }
      beta_new <- solve(Hessian) %*% t(X) %*% W %*% z
    }, error = function(e) {
      message("Error in matrix inversion: ", e$message)
      break
    })
    
    # Log-Likelihood
    log_likelihood <- sum(y * log(p_hat) + (1 - y) * log(1 - p_hat))
    
    # Convergence check
    if (abs(log_likelihood - log_likelihood_prev) < tol && sum(abs(beta_new - beta)) < tol) {
      message("Convergence achieved after ", i, " iterations.")
      break
    }
    
    if (i == max_iter) {
      warning("Maximum number of iterations reached without convergence.")
    }
    
    log_likelihood_prev <- log_likelihood
    beta <- beta_new
  }
  
  return(list(coefficients = beta, log_likelihood = log_likelihood, iterations = i))
}




#Initial Values for Optimization
initial_values <- function(X, y) {
  X <- as.matrix(cbind(1, X)) # Add intercept
  
  tryCatch({
    XtX <- t(X) %*% X
    if (det(XtX) == 0) {
      stop("Matrix inversion not possible due to singularity.")
    }
    return(solve(XtX) %*% t(X) %*% y)
  }, error = function(e) {
    message("Error in matrix inversion: ", e$message)
  })
}


#Bootstrap confidence intervals
bootstrap_ci <- function(X, y, alpha = 0.05, n_bootstraps = 20) {
  # Input validation
  if (!is.matrix(X) || length(y) != nrow(X)) {
    stop("Invalid inputs for X or y")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("Invalid alpha level")
  }
  if (n_bootstraps <= 0) {
    stop("Number of bootstraps must be positive")
  }
  
  # Initialize matrix to store bootstrap coefficients
  p <- ncol(X)
  bootstrapped_betas <- matrix(nrow = n_bootstraps, ncol = p + 1)
  
  # Bootstrap procedure
  for (i in 1:n_bootstraps) {
    sample_indices <- sample(nrow(X), replace = TRUE)
    # Extracting only the coefficients from the logistic regression output
    regression_result <- logistic_regression(X[sample_indices, ], y[sample_indices])
    bootstrapped_betas[i, ] <- regression_result$coefficients
  }
  
  # Calculate confidence intervals
  lower <- apply(bootstrapped_betas, 2, quantile, probs = alpha / 2)
  upper <- apply(bootstrapped_betas, 2, quantile, probs = 1 - alpha / 2)
  
  return(list(lower = lower, upper = upper))
}




#Plot of the fitted logistic curve
plot_logistic_curve <- function(X, y, beta, predictor_index = 1) {
  if (!is.matrix(X) || length(y) != nrow(X)) {
    stop("Invalid inputs for X or y")
  }
  if (predictor_index > ncol(X)) {
    stop("Predictor index out of bounds")
  }
  
  x_range <- range(X[, predictor_index])
  seq_x <- seq(from = x_range[1], to = x_range[2], length.out = 100)
  
  # Using column means for other predictors
  X_means <- apply(X, 2, mean)
  X_temp <- matrix(X_means, nrow = 100, ncol = ncol(X), byrow = TRUE)
  X_temp[, predictor_index] <- seq_x
  
  z <- cbind(1, X_temp) %*% beta
  p_hat <- 1 / (1 + exp(-z))
  
  # Adjusting plot's y-axis to range from 0 to 1
  plot(seq_x, p_hat, type = 'l', ylim = c(0, 1), col = 'blue', 
       xlab = names(X)[predictor_index], ylab = 'Predicted Probability')
  points(X[, predictor_index], y, col = 'red', pch = 20)
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
  # Ensure matrix has the correct dimensions
  if (!is.matrix(confusion_matrix) || nrow(confusion_matrix) != 2 || ncol(confusion_matrix) != 2) {
    stop("Invalid confusion matrix")
  }
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  return(accuracy)
}



# Function to calculate sensitivity
calculate_sensitivity <- function(confusion_matrix) {
  sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
  return(sensitivity)
}



# Function to calculate specificity
calculate_specificity <- function(confusion_matrix) {
  # Checking the structure of the confusion matrix
  if (!is.matrix(confusion_matrix) || nrow(confusion_matrix) != 2 || ncol(confusion_matrix) != 2) {
    stop("Invalid confusion matrix")
  }
  
  TN <- confusion_matrix[1, 1]  # True Negatives
  FP <- confusion_matrix[1, 2]  # False Positives
  
  if ((TN + FP) == 0) {
    return(NA)  # Avoid division by zero
  }
  
  specificity <- TN / (TN + FP)
  return(specificity)
}



# Function to calculate false discovery rate
calculate_false_discovery_rate <- function(confusion_matrix) {
  if (sum(confusion_matrix[, 2]) == 0) {
    return(NA) # Avoid division by zero
  }
  false_discovery_rate <- confusion_matrix[1, 2] / sum(confusion_matrix[, 2])
  return(false_discovery_rate)
}



# Function to calculate diagnostic odds ratio
calculate_diagnostic_odds_ratio <- function(confusion_matrix) {
  TP <- confusion_matrix[2, 2]
  TN <- confusion_matrix[1, 1]
  FP <- confusion_matrix[1, 2]
  FN <- confusion_matrix[2, 1]
  
  if (TP == 0 || TN == 0 || FP == 0 || FN == 0) {
    return(NA) # Avoid division by zero
  }
  
  diagnostic_odds_ratio <- (TP / FN) / (FP / TN)
  return(diagnostic_odds_ratio)
}


# Function to plot metrics over a grid of cutoff values
plot_selected_metrics_over_cutoffs <- function(fitted_probabilities, y, selected_metrics) {
  cutoff_values <- seq(0.1, 0.9, by = 0.1)  # Generate cutoff values
  metrics <- matrix(nrow = length(cutoff_values), ncol = 6)
  colnames(metrics) <- c("Prevalence", "Accuracy", "Sensitivity", "Specificity", "False Discovery Rate", "Diagnostic Odds Ratio")
  
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
      calculate_diagnostic_odds_ratio(confusion_matrix)
    )
  }
  
  # Plotting selected metrics
  for (metric in selected_metrics) {
    if (metric %in% colnames(metrics)) {
      plot(cutoff_values, metrics[, metric], type = "l", col = which(colnames(metrics) == metric),
           xlab = "Cutoff Value", ylab = metric)
      lines(cutoff_values, metrics[, metric], type = "l", col = which(colnames(metrics) == metric))
    }
  }
  legend("topright", legend = selected_metrics, col = 1:length(selected_metrics), lty = 1)
}

################################################################################
#Test functions defined above
################################################################################
# Generate a synthetic dataset
set.seed(123)  # For reproducibility
n <- 100       # Number of observations
p <- 2         # Number of predictors

# Generating predictors
X <- matrix(rnorm(n * p), ncol = p)
colnames(X) <- c("Predictor1", "Predictor2")

# Generating response
beta_true <- c(-0.5, 0.75, 1.25) # True coefficients including intercept
logits <- cbind(1, X) %*% beta_true
prob <- 1 / (1 + exp(-logits))
y <- rbinom(n, 1, prob)

# Testing logistic_regression
logistic_model <- logistic_regression(X, y)
print("Logistic Regression Model:")
print(logistic_model)

# Testing initial_values
initial_beta <- initial_values(X, y)
print("Initial Coefficients:")
print(initial_beta)

# Testing bootstrap_ci
ci <- bootstrap_ci(X, y)
print("Bootstrap Confidence Intervals:")
print(ci)

# Testing plot_logistic_curve
print("Plotting Logistic Curve:")
plot_logistic_curve(X, y, logistic_model$coefficients, predictor_index = 1)

# Generating fitted probabilities for testing other functions
fitted_probabilities <- 1 / (1 + exp(-cbind(1, X) %*% logistic_model$coefficients))

# Testing calculate_confusion_matrix and other metrics
conf_matrix <- calculate_confusion_matrix(y, ifelse(fitted_probabilities > 0.5, 1, 0))

print("Confusion Matrix:")
print(conf_matrix)

# Calculating and printing various metrics
prevalence <- calculate_prevalence(y)
accuracy <- calculate_accuracy(conf_matrix)
sensitivity <- calculate_sensitivity(conf_matrix)
specificity <- calculate_specificity(conf_matrix)
false_discovery_rate <- calculate_false_discovery_rate(conf_matrix)
diagnostic_odds_ratio <- calculate_diagnostic_odds_ratio(conf_matrix)

metrics_list <- list(
  prevalence = prevalence, 
  accuracy = accuracy, 
  sensitivity = sensitivity, 
  specificity = specificity, 
  false_discovery_rate = false_discovery_rate, 
  diagnostic_odds_ratio = diagnostic_odds_ratio
)

print("Metrics:")
print(metrics_list)

# Testing plot_selected_metrics_over_cutoffs for each metric
metrics_names <- c("Prevalence", "Accuracy", "Sensitivity", "Specificity", "False Discovery Rate", "Diagnostic Odds Ratio")
for (metric in metrics_names) {
  print(paste("Plotting", metric))
  plot_selected_metrics_over_cutoffs(fitted_probabilities, y, metric)
}
