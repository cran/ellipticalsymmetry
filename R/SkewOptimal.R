f = function(v, k, x) {
  val = (1 + x * x / v) ^ (-(v + k) / 2)
  return(val)
}

fprime = function(v, k, x) {
  # derivative of f
  val = (-(v + k) / 2) * (1 + x * x / v) ^ (-(v + k) / 2 - 1) * 2 * x /
    v
  return(val)
}

phi_t = function(v, k, x)
  # phi = -f'/f
  return(-fprime(v, k, x) / f(v, k, x))

phiprime_t = function(v, k, x) {
  val = (v + k) * (v - x * x) / (v + x * x) ^ 2
  return(val)
}

getK_t = function(v, d, norms, n) {
  #v=dof; k=dimension; d = norm (definisano je dolje u EllipticalSymLU funkciji); n = sample size
  summ = 0
  for (i in seq(1, n)) {
    summ = summ + phiprime_t(v, d, norms[i]) + (d - 1) * phi_t(v, d, norms[i]) / norms[i]
  }
  summ = summ / n
  return(summ)
}


phi_logistic = function(x){
  val = 2*x*(1-exp(-x^2))/(1+exp(-x^2))
  return(val)
}


phiprime_logistic = function(x) {
  val = (2*(1-exp(-x^2))*(1+exp(-x^2)) + 8*x^2*exp(-x^2))/(1+exp(-x^2))^2
  return(val)
}


getK_logistic = function(d, norms, n) {
  # k=dimension; d = norm (definisano je dolje u EllipticalSymLU funkciji); n = sample size
  summ = 0
  for (i in seq(1, n)) {
    summ = summ + phiprime_logistic(norms[i]) + (d - 1) * phi_logistic(norms[i]) / norms[i]
  }
  summ = summ / n
  return(summ)
}


phi_powerExp = function(lambda, x){
  val = lambda*x^(2*lambda -1)
  return(val)
}


phiprime_powerExp = function(lambda, x) {
  val = lambda*(2*lambda - 1)*x^(2*lambda - 2)
  return(val)
}


getK_powerExp = function(lambda, d, norms, n) {
  # k=dimension; d = norm (definisano je dolje u EllipticalSymLU funkciji); n = sample size
  summ = 0
  for (i in seq(1, n)) {
    summ = summ + phiprime_powerExp(lambda, norms[i]) + (d - 1) * phi_powerExp(lambda, norms[i]) / norms[i]
  }
  summ = summ / n
  return(summ)
}


SkewOptimalLK = function(X, location){
  n = dim(X)[1]
  d = dim(X)[2]

  sigma = tyler_cov(X, location)
  X = sweep(X, 2, location)
  sigma_inv = solve(sigma)
  theta = colMeans(X)
  statistic = n*t(theta)%*%sigma_inv%*%theta

  p_val = 1 - stats::pchisq(statistic, df = d)
  output = list(statistic = statistic[[1]], p.value = p_val[[1]])
  return(output)
}

SkewOptimalLU = function(X, f = 't', param = NaN) {
  # do this for 2.1, 4, 8 and 10  #simulacije su pokazale da su najbolji rez za 2.1 ili 4, to mozemo staviti kao input
  data_size = dim(X)
  n = data_size[1]
  d = data_size[2]

  theta = colMeans(X)
  #sigmahat=tyler.shape(X);
  sigma = tyler_cov(X, theta)

  sigma_root = spd_matrix_pow(sigma, -1/2)
  Z = sweep(X, 2, theta)%*%sigma_root
  #d = apply(Z, 1, norm, type = "2")
  norms = apply(Z, 1, function(x) norm(x, type='2'))
  U = Z / norms

  idmatr = diag(d)  # identity matrix
  gama = matrix(0, nrow = d, ncol = d)
  delta = matrix(0, nrow = d, ncol = d)
  if(f == 't'){
    if (any(is.nan(param))){
      param = 4
    }
    J = getK_t(param, d, norms, n)
    for (i in 1:n) {
      phi_val = phi_t(param, d, norms[i])
      gama = gama + (norms[i] - d * phi_val / J) ^ 2 / d * idmatr
      delta = delta + (norms[i] - d * phi_val / J) * U[i,]
    }
  }

  else if(f == 'logistic'){
    if (!any(is.nan(param))){
      warning('param is specified but it will not be used')
    }
    J = getK_logistic(d, norms, n)
    for (i in 1:n) {
      phi_val = phi_logistic(norms[i])
      gama = gama + (norms[i] - d * phi_val / J) ^ 2 / d * idmatr
      delta = delta + (norms[i] - d * phi_val / J) * U[i, ]
    }
  }

  else if(f == 'powerExp'){
    if (any(is.nan(param))){
      param = 0.5
    }
    J = getK_powerExp(param, d, norms, n)
    for (i in 1:n) {
      phi_val = phi_powerExp(param, norms[i])
      gama = gama + (norms[i] - d * phi_val / J) ^ 2 / d * idmatr
      delta = delta + (norms[i] - d * phi_val / J) * U[i, ]
    }
  }
  else{
    stop('f has to take one of the following three values: t, logistic, powerExp')
  }

  statistic = t(delta) %*% solve(gama) %*% delta
  p_val = 1 - stats::pchisq(statistic, df = d)
  output = list(statistic = statistic[[1]], p.value = p_val[[1]])
  return(output)
}


#' Tests for elliptical symmetry by Babic et al.
#'
#' @description Tests for elliptical symmetry: specified and unspecified location.
#'
#' @param X A numeric matrix.
#' @param location A vector of location parameters.
#' @param f A string that specifies the type of the radial density based on which the test is built. Currently supported options are \code{"t"}, \code{"logistic"} and \code{"powerExp"}.
#' The default is set to \code{"t"}.
#' @param param A parameter that is used when \code{f = "t"} and \code{f = "powerExp"}.
#' The default value of \code{param} represents the degrees of freedom of the multivariate t distribution and it is set to 4.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'  \item{\code{statistic}}{The value of the test statistic.}
#'  \item{\code{pvalue}}{The p-value of the test.}
#'  \item{\code{alternative}}{A character string describing the alternative hypothesis.}
#'  \item{\code{method}}{A character string indicating what type of test was performed.}
#'
#' @details
#' \code{X} and \code{location} are the only input arguments for the specified location test.
#' The default value for \code{location} is set to \code{NaN} which implies that the unspecified location test will be performed
#' unless the user specifies location.
#'
#' For the unspecified location test, besides the data matrix \code{X}, the input arguments are \code{f} and \code{param}.
#' The \code{f} argument is a string that specifies the type of the radial density based on which the test is built.
#' Currently supported options are: \code{"t"} for the radial density of the multivariate t distribution,
#' \code{"logistic"} for the multivariate logistic and \code{"powerExp"} for the radial density of the multivariate power-exponential distribution.
#' Note that the default is set to \code{"t"}.
#' The role of the \code{param} argument is as follows.
#' If \code{f = "t"} then \code{param} denotes the degrees of freedom of the multivariate t distribution.
#' Given that the default radial density is \code{"t"}, it follows that the default value of \code{param}
#' represents the degrees of freedom of the multivariate t distribution and it is set to 4.
#' Note also that the degrees of freedom have to be greater than 2.
#' If \code{f = "powerExp"} then \code{param} denotes the kurtosis parameter. In that case the value of \code{param}
#' has to be different from 1, because for the multivariate power exponential distribution, kurtosis parameter equal to 1 corresponds
#' to the multivariate Gaussian distribution (the Gaussian \code{f} is excluded due to a singular information matrix).
#' The default value is set to 0.5.
#'
#' @section Background:
#' Tests for elliptical symmetry both for specified and unspecified location. These tests are based on
#' Le Cam’s theory of statistical experiments and they are optimal against generalized skew-elliptical alternatives,
#' but they remain quite powerful under a much broader class of non-elliptical distributions.
#' The test statistic for the specified location scenario has a very simple form and an asymptotic chi-squared distribution.
#' When location is not specified, these tests have a simple asymptotic chi-squared
#' distribution under the null hypothesis of ellipticity, they are affine-invariant, computationally fast,
#' have a simple and intuitive form, only require finite moments of order 2.
#'
#'
#'
#' @references
#' Babic, S., Gelbgras, L., Hallin, M., & Ley, C. (2019). Optimal tests for elliptical symmetry: specified and unspecified location.
#'
#'
#' @examples
#'
#' ## sepal width and length of the versicolor subset of the Iris data
#' X = datasets::iris[51:100, 1:2]
#'
#'
#' ## location unspecifed test based on the radial density of the multivariate t distribution
#' SkewOptimal(X)
#'
#' ## location unspecifed test based on the radial density of the logistic distribution
#' SkewOptimal(X, f='logistic')
#'
#' ## location unspecifed test based the radial density of the power exponential distribution
#' SkewOptimal(X, f='powerExp')
#'
#' @export
SkewOptimal = function(X, location = NaN, f = 't', param = NaN) {

  dname = deparse(substitute(X))

  if(!is.matrix(X)) {
    warning_message = paste("coercing '", dname, "' to a matrix.", sep = "")
    warning(warning_message)
    X = as.matrix(X)
    if (!(is.matrix(X) && length(X) > 1)){
      stop("X is not in the valid matrix form.")
    }
  }

  else if(!is.numeric(X)){
    stop('X has to take numeric values')
  }


  if (any(is.nan(location))){
    output = SkewOptimalLU(X, f, param)
  }
  else{
    output = SkewOptimalLK(X, location)
  }
  names(output$statistic) = 'statistic'

  res <- list(method = 'SkewOptimal test for elliptical symmetry',
              data.name = dname,
              statistic = output$statistic,
              p.value = output$p.value,
              alternative = 'the distribution is not ellipticaly symmetric')
  class(res) <- "htest"
  return(res)
}



