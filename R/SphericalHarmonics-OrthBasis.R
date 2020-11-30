
func1deg1 = function(vec, i) {
  d = length(vec)
  c1 = d
  const = sqrt(c1)
  res = vec[i]
  return (const * res)
}

func1deg2 = function(vec, i, j) {
  d = length(vec)
  c1 = d * (d + 2)
  const = sqrt(c1)
  res = vec[i] * vec[j]
  return (const * res)
}

func2deg2 = function(vec, k) {
  d = length(vec)
  c1 = d * (d + 2)
  S = sum(vec[1:k - 1] ^ 2)
  const = sqrt(c1 / (2 * k * (k - 1)))
  res = S - (k - 1) * vec[k] ^ 2
  return (const * res)
}


func1deg3 = function(vec, i, j, k) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4)
  const = sqrt(c1)
  res = vec[i] * vec[j] * vec[k]
  return (const * res)
}

func2deg3 = function(vec, k, r) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4)
  S = sum(vec[1:k - 1] ^ 2)
  const = sqrt(c1 / (2 * k * (k - 1)))
  res = vec[r] * (S - (k - 1) * vec[k] ^ 2)
  return(const * res)
}

func3deg3 = function(vec, j, k) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4)
  S = sum(vec[1:k - 1] ^ 2)
  const = sqrt(c1 / (2 * (k + 1) * (k + 2)))
  res = vec[j] * (S - (k + 1) * vec[k] ^ 2)
  return(const * res)
}

func4deg3 = function(vec, r) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4)
  S = sum(vec[1:r - 1] ^ 2)
  const = sqrt(c1 / (6 * (r - 1) * (r + 2)))
  res = (r - 1) * vec[r] ^ 3 - 3 * vec[r] * S
  return(const * res)
}

func1deg4 = function(vec, r) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(vec[1:r - 1] ^ 2)
  const = sqrt(c1 / (24 * (r - 1) * (r + 1) * (r + 2) * (r + 4)))
  res = 6 * (r + 1) * vec[r] ^ 2 * S - (r ^ 2 - 1) * vec[r] ^ 4 - 3 * S ^
    2
  return(const * res)
}

func2deg4 = function(vec, i, r) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(vec[1:r - 1] ^ 2)
  const = sqrt(c1 / (6 * (r + 1) * (r + 4)))
  res = vec[r] * vec[i] * (3 * S - (r + 1) * vec[r] ^ 2)
  return(const * res)
}

func3deg4 = function(vec, i, j, r) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(vec[1:r - 1] ^ 2)
  const = sqrt(c1 / (2 * (r + 3) * (r + 4)))
  res = vec[i] * vec[j] * (S - (r + 3) * vec[r] ^ 2)
  return(const * res)
}
func4deg4 = function(vec, r, s) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  Sr = sum(vec[1:r - 1] ^ 2)
  Ss = sum(vec[1:s - 1] ^ 2)
  const = sqrt(c1 / (4 * r * (r - 1) * (s + 3) * (s + 4)))
  res = (Ss - (s + 3) * vec[s] ^ 2) * (Sr - (r - 1) * vec[r] ^ 2)
  return(const * res)
}

func5deg4 = function(vec, i, j, r, s) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  const = sqrt(c1)
  res = vec[i] * vec[j] * vec[r] * vec[s]
  return (const * res)
}

func6deg4 = function(vec, k, r, s) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(vec[1:k - 1] ^ 2)
  const = sqrt(c1 / (2 * k * (k - 1)))
  res = vec[r] * vec[s] * (S - (k - 1) * vec[k] ^ 2)
  return(const * res)
}

func7deg4 = function(vec, j, r, s) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(vec[1:r - 1] ^ 2)
  const = sqrt(c1 / (2 * (r + 1) * (r + 2)))
  res = vec[j] * vec[s] * (S - (r + 1) * vec[r] ^ 2)
  return(const * res)
}

func8deg4 = function(vec, r, s) {
  d = length(vec)
  c1 = d * (d + 2) * (d + 4) * (d + 6)
  S = sum(vec[1:r - 1] ^ 2)
  const = sqrt(c1 / (6 * (r - 1) * (r + 2)))
  res = vec[r] * vec[s] * (3 * S - (r - 1) * vec[r] ^ 2)
  return(const * res)
}
