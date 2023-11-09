#' @title Evaluate Linear Loss Function
#'
#' @description This function delivers the value of the loss function for a given balue of \code{beta}
#' @param resp A \code{vector} of dimension n.
#' @param pred A \code{matrix} containing predictors.
#' @param beta A \code{vector} containing coefficients.
#' @param norm A \code{character} defining loss to use (default is `L2`).
#' @return A \code{numeric} giving value of loss at \code{beta}
#' @author Rob Molinari
#' @importFrom stats
#' @export
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(600), nrow = 200) # true matrix of predictors
#' epsilon <- rnorm(200, 0, sd = 0.25) # observation error
#' X <- cbind(rep(1, 200), X) # add a column of ones to generate the response with an intercept value
#' beta <- c(-1, 4, -5, 2) # true coefficient values we want to estimate (first element is the intercept)
#' y <- X%*%beta + epsilon # generate response y based on this model
#' beta_ls(resp = y, pred = X, beta = beta, norm = "L2")
beta_ls <- function(resp, pred, beta, norm = "L2") {
  
  if(norm == "L2") {
    
    out <- t(resp - pred%*%beta)%*%(resp - pred%*%beta) # L2-norm (i.e. least-squares)
    
  } else {
    
    out <- sum(abs(resp - pred%*%beta)) # L1-norm (robust but inefficient)
    
  }
  
  return(out)
  
}