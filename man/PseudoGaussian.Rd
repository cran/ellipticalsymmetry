% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PseudoGaussian.R
\name{PseudoGaussian}
\alias{PseudoGaussian}
\title{Pseudo-Gaussian test for elliptical symmetry}
\usage{
PseudoGaussian(X, location = NA)
}
\arguments{
\item{X}{A numeric matrix.}

\item{location}{A vector of location parameters.}
}
\value{
An object of class \code{"htest"} containing the following components:
\item{\code{statistic}}{The value of the test statistic.}
\item{\code{pvalue}}{The p-value of the test.}
\item{\code{alternative}}{A character string describing the alternative hypothesis.}
\item{\code{method}}{A character string indicating what type of test was performed.}
}
\description{
Tests for elliptical symmetry: specified and unspecified location.
}
\details{
Note that \code{location} allows the user to specify the known location.
The default is set to \code{NA} which means that the unspecified location test will be performed unless the user specifies location.
}
\section{Background}{

Pseudo-Gaussian tests for elliptical symmetry are based on Le Cam’s theory of statistical experiments.
They are most efficient against a multivariate form of Fechner-type asymmetry.
These tests require finite moments of order 4 and they have a simple asymptotic chi-squared distribution
under the null hypothesis of ellipticity.
}

\examples{

## sepal width and length of the versicolor subset of the Iris data
X = datasets::iris[51:100,1:2]

PseudoGaussian(X)

}
\references{
Cassart, D., Hallin, M. & Paindaveine, D., (2008). Optimal detection of Fechner-asymmetry. \emph{Journal of Statistical Planning and Inference}, \bold{138}, 2499-2525.

Cassart, D., (2007). Optimal tests for symmetry. Ph.D. thesis, Univ. libre de Bruxelles, Brussels.
}
