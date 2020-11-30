
getAngle = function(x, y){
  val = atan2(y, x)
  if(val < 0){
    val = 2*pi + val
  }
  return(val)
}

unispehre = function(n, d){
  points = matrix(0, nrow = n, ncol = d)
  for (i in seq(d)){
    samp = stats::rnorm(n)
    for (j in seq(n)){
      points[j,i] = samp[j]
    }
  }
  norms = apply(points, 1, function(x) norm(x, type='2'))
  points = points/norms
  return (points)
}

getPermutationIndex = function(permutation){
  perm_string = paste(permutation, collapse = '')
  string_size=nchar(perm_string)
  if(string_size == 1){
    return(1)
  }
  ind=0
  for (i in seq(1, string_size)){
    if (substr(perm_string, i, i) == paste(string_size)){
      ind=i
      break
    }
  }
  new_perm = paste(substr(perm_string, 1, ind-1), substr(perm_string, ind + 1, string_size), sep = '')
  return((ind - 1)*factorial(string_size - 1) + getPermutationIndex(new_perm))
}

getExactPValue = function(d, c, statistic){
  cutoffs = stats::qchisq(seq(0, 1, length.out = c + 1), df = d)
  a = replicate(c, 0)
  b = replicate(c, 0)
  for (i in seq(1, c)){
    a[i] = stats::pchisq(cutoffs[i + 1], df = d + 1) - stats::pchisq(cutoffs[i] , df = d + 1)
    b[i] = stats::pchisq(cutoffs[i + 1], df = d + 2) - stats::pchisq(cutoffs[i] , df = d + 2)
  }
  astar = 2*c*sum(a^2)/pi
  bstar = 4*c*sum(b^2)/(pi*pi)
  df1 = c*(2^d - 1) - d*(d + 1)/2
  df2 = d
  df3 = d*(d - 1)/2
  sample_size = 1000000
  sample1 = stats::rchisq(sample_size, df1)
  sample2 = stats::rchisq(sample_size, df2)
  sample3 = stats::rchisq(sample_size, df3)
  sample = sample1 + (1 - astar)*sample2 + (1 - bstar)*sample3
  emp_cdf = stats::ecdf(sample)
  p_val = 1 - emp_cdf(statistic)
  return(p_val)
}

bitsToInt = function(x) {
  packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
}


getstatisticHP = function(X, g, c, sector){
  data_size = dim(X)
  n = data_size[1]
  d = data_size[2]

  theta = colMeans(X)
  sigma = (n - 1)*stats::cov(X)/n
  L = solve(chol(sigma))
  Z = sweep(X, 2, theta)%*%L
  norms = apply(Z, 1, function(x) norm(x, type='2'))
  emp_cdf = stats::ecdf(norms)
  if (sector == 'orthants'){
    counts = matrix(0, nrow=2 ^ d, ncol = c)
    for (i in seq(n)){
      sample = Z[i,]
      bin = (sign(sample) + 1)/2
      ind1 = bitsToInt(bin) + 1
      sample_norm = norm(sample, type='2')
      ind2 = min(floor(emp_cdf(sample_norm)*c) + 1, c)
      counts[ind1, ind2] = counts[ind1, ind2] + 1

    }
    expected = n/(2^d*c)
  }
  if (sector == 'bivariate_angles'){
    counts = matrix(0, nrow=g, ncol = c)
    for (i in seq(n)){
      sample = Z[i,]
      angle = getAngle(sample[1], sample[2])
      ind1 = min(floor(angle/(2*pi)*g) + 1, g)
      sample_norm = norm(sample, type='2')
      ind2 = min(floor(emp_cdf(sample_norm)*c) + 1, c)
      counts[ind1, ind2] = counts[ind1, ind2] + 1
    }
    expected = n/(g*c)
  }
  if (sector == 'permutation'){
    fact = factorial(d)
    counts = matrix(0, nrow=fact, ncol = c)
    for (i in seq(n)){
      sample = Z[i,]
      ind1 = getPermutationIndex(order(sample))
      sample_norm = norm(sample, type='2')
      ind2 = min(floor(emp_cdf(sample_norm)*c) + 1, c)
      counts[ind1, ind2] = counts[ind1, ind2] + 1
    }
    expected = n/(fact*c)
  }
  statistic = 0
  for (count in counts){
    statistic = statistic + (expected - count)^2/expected
  }
  return(statistic)
}

bootstrapHP = function(X, g, c, sector, R, statistic){
  data_size = dim(X)
  n = data_size[1]
  d = data_size[2]

  theta = colMeans(X)
  sigma = (n - 1)*stats::cov(X)/n
  L = solve(chol(sigma))
  norms = apply(sweep(X, 2, theta)%*%L, 1, function(x) norm(x, type='2'))
  bootstrap_statistic = replicate(R, 0)
  for (i in seq(R)){
    U = unispehre(n, d=d)
    boot_norms = sample(norms, replace=T)
    Xsim = U*boot_norms
      bootstrap_statistic[i] = getstatisticHP(Xsim, g, c, sector)
  }
  emp_cdf = stats::ecdf(bootstrap_statistic)
  p_val = 1 - emp_cdf(statistic)
  return(p_val)
}

#' Huffer and Park's test for elliptical symmetry
#'
#' @description Pearson chi-square type test for elliptical symmetry.
#'
#' @param X A numeric matrix.
#' @param c The number of spherical shells that are used to divide the space.
#' @param R The number of bootstrap replications.
#' @param sector A string that specifies the type of sectors used to divide the space. Currently supported options are \code{"orthants"}, \code{"permutation"} and \code{"bivariate_angles"}.
#' @param g A parameter that is used if \code{sector = "bivariate_angles"}. It denotes the number of regions used to divide the plane.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'  \item{\code{statistic}}{The value of the test statistic.}
#'  \item{\code{pvalue}}{The p-value of the test.}
#'  \item{\code{alternative}}{A character string describing the alternative hypothesis.}
#'  \item{\code{method}}{A character string indicating what type of test was performed.}
#'
#' @details
#' Huffer and Park (2007) propose a Pearson chi-square type test with multi-dimensional cells.
#' After dividing the space into \code{c} spherical shells  and \code{g} sectors (in total \code{gc} cells),
#' and after determining the observed cell counts, the test statistic is easily computed.
#' \code{sectors} is an option that allows the user to specify the type of sectors used to divide the space.
#' Currently supported options are \code{"orthants"}, \code{"permutation"} and \code{"bivariate_angles"},
#' the last one being available only in dimension 2. The \code{g} argument indicates the number of sectors.
#' The user has to choose \code{g} only if \code{sectors = "bivariate_angles"} and it denotes the number of regions used to divide the plane.
#' In this case, regions consist of points whose angle in polar coordinates is between \eqn{2(m-1)\pi/g} and \eqn{2m\pi/g} for \eqn{m} in \eqn{(1,..., g)}.
#' If \code{sectors} is set to \code{"orthants"}, then \code{g} is fixed and equal to \eqn{2^d}, while for \code{sectors = "permutation"} \code{g} is \eqn{d}!.
#' No matter what type of sectors is chosen, the user has to specify the number of spherical shells that are used to divide the space, which is \code{c}.
#' The value of \code{c} should be such that the average cell counts \eqn{n/(gc)} are not too small.
#'
#' The asymptotic distribution is available only under \code{sectors = "orthants"} when the underlying distribution is close to normal.
#' Otherwise, bootstrap procedures are required and the user can freely choose the number of bootstrap replicates, denoted as \code{R}.
#' Note that by default \code{sectors} is set to \code{"orthants"} and \code{R = NaN}, which means that the non-bootstrap version of the test
#' will be performed unless the user specifies \code{R}.
#'
#' @references
#' Huffer, Fred W., & Park, C., (2007). A test for elliptical symmetry. \emph{Journal of Multivariate Analysis}, \bold{98}(2), 256-281.
#'
#' @examples
#'
#'
#' ## sepal width and length of the versicolor subset of the Iris data
#' X = datasets::iris[51:100, 1:2]
#'
#' ## the non-bootstrap test
#' HufferPark(X, c = 2)
#'
#' ## the bootstrap tests
#' HufferPark(X, c = 2, R = 10, sector='orthants')
#'
#' HufferPark(X[,1:2], c = 3, R = 10, sector='bivariate_angles', g = 3)
#'
#' HufferPark(X, c = 2, R = 10, sector='permutation')
#' @export
HufferPark = function(X, c, R = NaN, sector='orthants', g=NaN){

  if (is.nan(R)){
    sector = 'orthants'
  }
  else if (!sector %in% c('orthants', 'bivariate_angles', 'permutation')){
    stop('when the bootstrap is chosen (R is not NaN), sector has to take one of the following three values: orthants, bivariate_angles, permutation')
  }

  if (sector == 'bivariate_angles' && is.nan(g)){
    stop('g has to be specified.')
  }

  else if (sector != 'bivariate_angles' && !is.nan(g)){
    warning('g is specified but it will not be used.')
  }

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

  statistic = getstatisticHP(X, g, c, sector)[[1]]

  k = dim(X)[2]

  if (is.nan(R)){
    p_val = getExactPValue(k, c, statistic)[[1]]

  }
  else{
    p_val = bootstrapHP(X, g, c, sector, R, statistic)[[1]]
  }
  names(statistic) = 'statistic'
  res <- list(method = 'Test for elliptical symmetry by Huffer and Park',
              data.name = dname,
              statistic = statistic,
              p.value = p_val,
              alternative = 'the distribution is not ellipticaly symmetric')
  class(res) <- "htest"
  return(res)
}


