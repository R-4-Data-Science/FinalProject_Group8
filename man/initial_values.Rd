% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions for Package.R
\name{initial_values}
\alias{initial_values}
\title{Initial Values for Logistic Regression Optimization}
\usage{
initial_values(X, y)
}
\arguments{
\item{X}{A matrix of predictor variables.}

\item{y}{A vector of response variables.}
}
\value{
A vector of initial values for optimization.
}
\description{
Computes the initial values for the optimization of logistic regression coefficients using least squares.
}
\examples{
X <- matrix(rnorm(200), ncol = 2)
y <- rbinom(100, 1, 0.5)
init_vals <- initial_values(X, y)
}
