#' Logistic Regression Using Numerical Optimization
#'
#' Performs logistic regression using numerical optimization without relying on existing logistic regression functions in R.
#'
#' @param X A matrix of predictor variables.
#' @param y A vector of binary response variables.
#' @param max_iter Maximum number of iterations for optimization (default: 100).
#' @param tol Tolerance for convergence check (default: 1e-6).
#'
#' @return A list containing coefficients, log-likelihood, and number of iterations.
#' @examples
#' X <- matrix(rnorm(200), ncol = 2)
#' y <- rbinom(100, 1, 0.5)
#' result <- logistic_regression(X, y)
#' @export
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




#' Initial Values for Logistic Regression Optimization
#'
#' Computes the initial values for the optimization of logistic regression coefficients using least squares.
#'
#' @param X A matrix of predictor variables.
#' @param y A vector of response variables.
#'
#' @return A vector of initial values for optimization.
#' @examples
#' X <- matrix(rnorm(200), ncol = 2)
#' y <- rbinom(100, 1, 0.5)
#' init_vals <- initial_values(X, y)
#' @export
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


#' Bootstrap Confidence Intervals for Logistic Regression
#'
#' Calculates bootstrap confidence intervals for logistic regression coefficients.
#'
#' @param X A matrix of predictor variables.
#' @param y A vector of response variables.
#' @param alpha Significance level for confidence intervals (default: 0.05).
#' @param n_bootstraps Number of bootstrap samples (default: 20).
#'
#' @return A list containing lower and upper confidence intervals.
#' @examples
#' X <- matrix(rnorm(200), ncol = 2)
#' y <- rbinom(100, 1, 0.5)
#' ci <- bootstrap_ci(X, y)
#' @export
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




#' Plot the Fitted Logistic Regression Curve
#'
#' Plots the logistic regression curve based on fitted model coefficients.
#'
#' @param X A matrix of predictor variables.
#' @param y A vector of response variables (binary).
#' @param beta A vector of estimated coefficients from logistic regression.
#' @param predictor_index Index of the predictor to be used for plotting (default: 1).
#'
#' @examples
#' X <- matrix(rnorm(200), ncol = 2)
#' y <- rbinom(100, 1, 0.5)
#' beta <- logistic_regression(X, y)$coefficients
#' plot_logistic_curve(X, y, beta)
#' @export
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



#' Calculate the Confusion Matrix
#'
#' Computes the confusion matrix for binary classification predictions.
#'
# Function to calculate the confusion matrix
#' Title
#'
#'  @param y True binary responses.
#' @param predictions Binary predictions from the model.
#' @return A confusion matrix.
#' @examples
#' true_values <- rbinom(100, 1, 0.5)
#' predicted_values <- rbinom(100, 1, 0.5)
#' conf_matrix <- calculate_confusion_matrix(true_values, predicted_values)
#' @export
calculate_confusion_matrix <- function(y, predictions) {
  confusion_matrix <- table(y, predictions)
  return(confusion_matrix)
}




#' Calculate Prevalence
#'
#' Computes the prevalence of the positive class in binary response data.
#'
#' @param y Binary response vector.
#'
#' @return Prevalence value.
#' @examples
#' response <- rbinom(100, 1, 0.7)
#' prevalence <- calculate_prevalence(response)
#' @export
calculate_prevalence <- function(y) {
  prevalence <- sum(y) / length(y)
  return(prevalence)
}



#' Calculate Accuracy of Predictions
#'
#' Computes the accuracy from the confusion matrix for binary classification.
#'
#' @param confusion_matrix A confusion matrix from binary classification predictions.
#'
#' @return Accuracy value.
#' @examples
#' predictions <- rbinom(100, 1, 0.5)
#' true_labels <- rbinom(100, 1, 0.5)
#' conf_matrix <- calculate_confusion_matrix(true_labels, predictions)
#' accuracy <- calculate_accuracy(conf_matrix)
#' @export
calculate_accuracy <- function(confusion_matrix) {
  # Ensure matrix has the correct dimensions
  if (!is.matrix(confusion_matrix) || nrow(confusion_matrix) != 2 || ncol(confusion_matrix) != 2) {
    stop("Invalid confusion matrix")
  }
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  return(accuracy)
}


#' Calculate Sensitivity
#'
#' Computes sensitivity from the confusion matrix for binary classification.
#'
# Function to calculate sensitivity
#' Title
#'
#' @param confusion_matrix A confusion matrix from binary classification predictions.
#' @return Sensitivity value.
#' @examples
#' predictions <- rbinom(100, 1, 0.5)
#' true_labels <- rbinom(100, 1, 0.5)
#' conf_matrix <- calculate_confusion_matrix(true_labels, predictions)
#' sensitivity <- calculate_sensitivity(conf_matrix)
calculate_sensitivity <- function(confusion_matrix) {
  sensitivity <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
  return(sensitivity)
}



#' Calculate Specificity
#'
#' Computes specificity from the confusion matrix for binary classification.
#'
#' @param confusion_matrix A confusion matrix from binary classification predictions.
#'
#' @return Specificity value.
#' @examples
#' predictions <- rbinom(100, 1, 0.5)
#' true_labels <- rbinom(100, 1, 0.5)
#' conf_matrix <- calculate_confusion_matrix(true_labels, predictions)
#' specificity <- calculate_specificity(conf_matrix)
#' @export
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



#' Calculate False Discovery Rate
#'
#' Computes the false discovery rate from the confusion matrix for binary classification.
#'
#' @param confusion_matrix A confusion matrix from binary classification predictions.
#'
#' @return False Discovery Rate value.
#' @examples
#' predictions <- rbinom(100, 1, 0.5)
#' true_labels <- rbinom(100, 1, 0.5)
#' conf_matrix <- calculate_confusion_matrix(true_labels, predictions)
#' fdr <- calculate_false_discovery_rate(conf_matrix)
#' @export
calculate_false_discovery_rate <- function(confusion_matrix) {
  if (sum(confusion_matrix[, 2]) == 0) {
    return(NA) # Avoid division by zero
  }
  false_discovery_rate <- confusion_matrix[1, 2] / sum(confusion_matrix[, 2])
  return(false_discovery_rate)
}



#' Calculate Diagnostic Odds Ratio
#'
#' Computes the diagnostic odds ratio from the confusion matrix for binary classification.
#'
#' @param confusion_matrix A confusion matrix from binary classification predictions.
#'
#' @return Diagnostic Odds Ratio value.
#' @examples
#' predictions <- rbinom(100, 1, 0.5)
#' true_labels <- rbinom(100, 1, 0.5)
#' conf_matrix <- calculate_confusion_matrix(true_labels, predictions)
#' dor <- calculate_diagnostic_odds_ratio(conf_matrix)
#' @export
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


#' Plot Selected Metrics Over Various Cutoff Values
#'
#' Plots selected metrics (e.g., accuracy, sensitivity) evaluated over a range of cutoff values.
#'
#' @param fitted_probabilities Predicted probabilities from the logistic regression model.
#' @param y True binary responses.
#' @param selected_metrics A vector of selected metric names to plot.
#'
#' @examples
#' fitted_probs <- logistic_regression(X, y)$fitted_values
#' plot_selected_metrics_over_cutoffs(fitted_probs, y, c("Accuracy", "Sensitivity"))
#' @export
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







