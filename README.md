Logistic Regression Package - Group 8
================
Lucas Parvin, Priyadarshni Jitendrabhai Patel, & Felix Satognon
2023-11-30

# Project 8 Final Project Package

## How Our Package Works:

### Logistic Regression Using Numerical Optimization

This function performs logistic regression using numerical optimization.
It does not rely on existing logistic regression functions in R.

#### Example Usage

``` r
#NOTES
# Don't include any code

#Having saved this file within the main folder of your package, you can compile this document using the pkgdown package (this can be installed by executing the following command in your console: devtools::install_github("hadley/pkgdown")). When the latter package is installed and loaded, all you need to do is execute the command pkgdown::build_site() to obtain the following web-page:

# Include link to website after it has been created

# Follow all steps in section 8.1

# Only include descriptions of funIctions that user will use.

#Overview

#This R package provides a set of functions for performing logistic regression, evaluating model performance, and visualizing results. The package is designed to be user-friendly, allowing users to perform logistic regression without relying on existing functions in R. Additionally, it includes tools for computing bootstrap confidence intervals, plotting the logistic regression curve, calculating confusion matrices, and more.

#Functions

#logistic_regression


#This function performs logistic regression using numerical optimization. Users can input a matrix of predictor variables (X) and a vector of binary response variables (y). Additional parameters include max_iter for the maximum number of iterations and tol for the convergence tolerance.



#initial_values

#Computes the initial values for logistic regression optimization using least squares. Users provide a matrix of predictor variables (X) and a vector of response variables (y).

#bootstrap_ci

#Calculates bootstrap confidence intervals for logistic regression coefficients. Users input a matrix of predictor variables (X), a vector of response variables (y), and optional parameters alpha for the significance level and n_bootstraps for the number of bootstrap samples.


#plot_logistic_curve

#Plots the logistic regression curve based on fitted model coefficients. Users provide a matrix of predictor variables (X), a vector of response variables (y), and the estimated coefficients from logistic regression (beta).

#calculate_confusion_matrix

#Computes the confusion matrix for binary classification predictions. Users input true binary responses (y) and predicted values (predictions).
```
