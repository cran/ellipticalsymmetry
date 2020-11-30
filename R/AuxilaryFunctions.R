spd_matrix_pow = function(A, pow){
  if(mean(abs(A - t(A))) > 1e-6){
    stop("The matrix is not symmetric.")
  }

  eig = eigen(A)
  for (val in eig$values){
    if (val < 0){
      stop("The matrix is not possitive semi-definite.")
    }
  }
  new_eig_vals = eig$values ^ pow
  Q = eig$vectors
  D = diag(new_eig_vals)
  res = Q %*% D %*% t(Q)
  return(res)

}



seq2 = function(a,b){
  if (b >= a)  {
    return (seq(a, b))
  }
  else return()
}

unispehre = function(n, d){
  points = matrix(0, nrow = n, ncol = d)
  for (i in seq(d)){
    samp = rnorm(n)
    for (j in seq(n)){
      points[j,i] = samp[j]
    }
  }
  norms = apply(points, 1, function(x) norm(x, type='2'))
  points = points/norms
  return (points)
}

tyler_cov = function(X, location){
  n = dim(X)[1]
  k = dim(X)[2]
  sigmahat = ICSNP::tyler.shape(X, location=location);
  b=array(0,dim=c(1,n));
  invsigma=solve(sigmahat);
  for (i in 1:n){b[i]=(X[i,]-location)%*%invsigma%*%(X[i,]-location)}
  sigmahat=sigmahat*mean(b)/k
  return(sigmahat)
}
