% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions for Package.R
\name{calculate_sensitivity}
\alias{calculate_sensitivity}
\title{Calculate Sensitivity}
\usage{
calculate_sensitivity(confusion_matrix)
}
\arguments{
\item{confusion_matrix}{A confusion matrix from binary classification predictions.}
}
\value{
Sensitivity value.
}
\description{
Computes sensitivity from the confusion matrix for binary classification.
}
\details{
Title
}
\examples{
predictions <- rbinom(100, 1, 0.5)
true_labels <- rbinom(100, 1, 0.5)
conf_matrix <- calculate_confusion_matrix(true_labels, predictions)
sensitivity <- calculate_sensitivity(conf_matrix)
}
