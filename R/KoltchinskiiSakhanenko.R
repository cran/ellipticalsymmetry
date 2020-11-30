
getstatisticKS = function(X){


  data_size = dim(X)
  n = data_size[1]
  d = data_size[2]

  theta = colMeans(X)
  sigma = (n - 1)*stats::cov(X)/n
  sigma_root = spd_matrix_pow(sigma, -1/2)
  Z = sweep(X, 2, theta)%*%sigma_root
  norms = apply(Z, 1, function(x) norm(x, type='2'))
  norms_order = order(norms)

  deg1 = d
  deg2 = choose(d + 1, 2) - 1
  deg3 = choose(d + 2, 3) - d
  deg4 = choose(d + 3, 4) - choose(d + 1, 2)

  deg = deg1 + deg2 + deg3 + deg4
  func_eval = matrix(0, deg, n)

  indexi = 1
  for (i in seq2(1, d)){
    indexj = 1
    for (w in norms_order){
      sample = Z[w,]
      sample_norm = norm(sample, type='2')
      func_eval[indexi,indexj] = func1deg1(sample/sample_norm, i)
      indexj = indexj + 1
    }
    indexi = indexi + 1
  }
  for (i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')
        func_eval[indexi,indexj] = func1deg2(sample/sample_norm, i, j)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }

  for (k in seq2(2, d)){
    indexj = 1
    for (w in norms_order){
      sample = Z[w,]
      sample_norm = norm(sample, type='2')
      func_eval[indexi,indexj] = func2deg2(sample/sample_norm, k)
      indexj = indexj + 1
    }
    indexi = indexi + 1
  }


  for(i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      for (k  in seq2(j + 1, d)){
        indexj = 1
        for (w in norms_order){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')

          func_eval[indexi,indexj] = func1deg3(sample/sample_norm, i, j, k)
          indexj = indexj + 1
        }
        indexi = indexi + 1
      }
    }
  }
  for(k in seq2(2, d)){
    for (r  in seq2(k + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')

        func_eval[indexi,indexj] = func2deg3(sample/sample_norm, k, r)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }

  for(j in seq2(1, d)){
    for (k  in seq2(j + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')

        func_eval[indexi,indexj] = func3deg3(sample/sample_norm, j, k)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }
  for (r  in seq2(2, d)){
    indexj = 1
    for (w in norms_order){
      sample = Z[w,]
      sample_norm = norm(sample, type='2')

      func_eval[indexi,indexj] = func4deg3(sample/sample_norm, r)
      indexj = indexj + 1
    }
    indexi = indexi + 1
  }

  for (r  in seq2(2, d)){
    indexj = 1
    for (w in norms_order){
      sample = Z[w,]
      sample_norm = norm(sample, type='2')

      func_eval[indexi,indexj] = func1deg4(sample/sample_norm, r)
      indexj = indexj + 1
    }
    indexi = indexi + 1
  }

  for(i in seq2(1, d)){
    for (r  in seq2(i + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')

        func_eval[indexi,indexj] = func2deg4(sample/sample_norm, i, r)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }
  for(i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      for (r  in seq2(j + 1, d)){
        indexj = 1
        for (w in norms_order){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')

          func_eval[indexi,indexj] = func3deg4(sample/sample_norm, i, j, r)
          indexj = indexj + 1
        }
        indexi = indexi + 1
      }
    }
  }
  for(r in seq2(2, d)){
    for (s  in seq2(r + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')

        func_eval[indexi,indexj] = func4deg4(sample/sample_norm, r, s)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }
  for(i in seq2(1, d)){
    for (j  in seq2(i + 1, d)){
      for (r  in seq2(j + 1, d)){
        for (s  in seq2(r + 1, d)){
          indexj = 1
          for (w in norms_order){
            sample = Z[w,]
            sample_norm = norm(sample, type='2')

            func_eval[indexi,indexj] = func5deg4(sample/sample_norm, i, j, r, s)
            indexj = indexj + 1
          }
          indexi = indexi + 1
        }
      }
    }
  }
  for (k  in seq2(2, d)){
    for (r  in seq2(k + 1, d)){
      for (s  in seq2(r + 1, d)){
        indexj = 1
        for (w in norms_order){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')

          func_eval[indexi,indexj] = func6deg4(sample/sample_norm, k, r, s)
          indexj = indexj + 1
        }
        indexi = indexi + 1
      }
    }
  }
  for (j  in seq2(1, d)){
    for (r  in seq2(j + 1, d)){
      for (s  in seq2(r + 1, d)){
        indexj = 1
        for (w in norms_order){
          sample = Z[w,]
          sample_norm = norm(sample, type='2')

          func_eval[indexi,indexj] = func7deg4(sample/sample_norm, j, r, s)
          indexj = indexj + 1
        }
        indexi = indexi + 1
      }
    }
  }
  for(r in seq2(2, d)){
    for (s  in seq2(r + 1, d)){
      indexj = 1
      for (w in norms_order){
        sample = Z[w,]
        sample_norm = norm(sample, type='2')

        func_eval[indexi,indexj] = func8deg4(sample/sample_norm, r, s)
        indexj = indexj + 1
      }
      indexi = indexi + 1
    }
  }
  statistic = 0

  func_eval_sums = matrix(0, deg, n)

  for (i in seq(deg)){
    val = 0
    for (j in seq(n)){
      val = val + func_eval[i,j]
      func_eval_sums[i,j] = val^2
    }
  }

  statistic = sqrt(max(colSums(func_eval_sums))/n)
  return(statistic)

}



bootstrapKS = function(X, R, statistic){
  data_size = dim(X)
  n = data_size[1]
  d = data_size[2]

  theta = colMeans(X)
  sigma = (n - 1)*stats::cov(X)/n
  sigma_root = spd_matrix_pow(sigma, -1/2)
  Z = sweep(X, 2, theta)%*%sigma_root
  norms = apply(Z, 1, function(x) norm(x, type='2'))
  bootstrap_statistic = replicate(R, 0)
  for (i in seq(R)){
    U = unispehre(n, d=d)
    boot_norms = sample(norms, replace=T)
    Xsim = U*boot_norms
    bootstrap_statistic[i] = getstatisticKS(Xsim)
  }
  emp_cdf = stats::ecdf(bootstrap_statistic)
  p_val = 1 - emp_cdf(statistic)
  return(p_val)
}

#' Koltchinskii and Sakhanenko's test for elliptical symmetry
#'
#' @description Test for elliptical symmetry.
#'
#' @param X A numeric matrix.
#' @param R The number of bootstrap replications.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'  \item{\code{statistic}}{The value of the test statistic.}
#'  \item{\code{pvalue}}{The p-value of the test.}
#'  \item{\code{alternative}}{A character string describing the alternative hypothesis.}
#'  \item{\code{method}}{A character string indicating what type of test was performed.}
#'
#' @section Background:
#' Koltchinskii and Sakhanenko (2000) proposed a class of omnibus bootstrap tests for elliptical symmetry
#' that are affine invariant and consistent against any fixed alternative. This test is based on spherical harmonics.
#'
#' @references
#' Koltchinskii, V., & Sakhanenko, L., (2000). Testing for ellipsoidal symmetry of a multivariate distribution. \emph{High Dimensional Probability II}, 493-510, Springer.
#'
#' Sakhanenko, L., (2008). Testing for ellipsoidal symmetry: A comparison study. \emph{Computational Statistics & Data Analysis}, \bold{53}(2), 565-581.
#'
#' @examples
#'
#' ## sepal width and length of the versicolor subset of the Iris data
#' X = datasets::iris[51:100, 1:2]
#'
#' KoltchinskiiSakhanenko(X, R = 10)
#' @export
KoltchinskiiSakhanenko = function(X, R=1000){

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

  statistic = getstatisticKS(X)[[1]]
  p_val = bootstrapKS(X, R, statistic)[[1]]
  names(statistic) = 'statistic'
  res <- list(method = 'Test for elliptical symmetry by Koltchinskii and Sakhanenko',
              data.name = dname,
              statistic = statistic,
              p.value = p_val,
              alternative = 'the distribution is not ellipticaly symmetric')
  class(res) <- "htest"
  return(res)
}
