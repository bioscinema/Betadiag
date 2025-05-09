% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_distance.R
\name{check_distance}
\alias{check_distance}
\title{Check distance matrix properties}
\usage{
check_distance(D)
}
\arguments{
\item{D}{A symmetric distance matrix (numeric matrix).}
}
\value{
A list with the following components:
\describe{
  \item{is.metric}{Whether the distance matrix satisfies the triangle inequality (0/1).}
  \item{is.Euclidean}{Whether the distance matrix corresponds to a Euclidean space (0/1).}
  \item{collinearity.score}{A numeric vector of collinearity scores for each sample.}
  \item{nonlinearity.score}{A numeric vector of nonlinearity scores for each sample.}
  \item{FNI}{Fraction of negative eigenvalues (from Gram matrix), same as `rate` in `Euclidean.Check()`.}
}
}
\description{
This function checks whether a given distance matrix is metric and Euclidean, 
and computes collinearity and nonlinearity scores, fraction of negative inertia (FNI) for each sample.
}
