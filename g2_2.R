## fix lambda0 = 3
g2 <- function(X, Y){
  ## step 1: data preparation
  idx = order(X)
  x = X[idx]
  y = Y[idx]
  # normalize
  y = y - mean(y)
  y = sqrt(n)*y/sqrt(sum(y^2))
  ## step 2: main algorithm
  n = length(X)
  m = ceiling(sqrt(n))
  lambda = -3*log(n)/2
  alpha = exp(lambda)
  # initialize three sequences
  Mi = numeric(n)
  Bi = numeric(n)
  Ti = numeric(n)
  Mi[1] = 0
  Bi[1] = Ti[1] = 1
  for (i in m:n){
    bi = 0
    ti = 0
    if (i < 2*m)
    { # do not ignore
      Mi[i] = Mi[1]
      Bi[i] = Bi[1]
      Ti[i] = Ti[1]
      next
    }
    #mi = numeric(i-2*m+1)
    mi = rep(-Inf, i-2*m+1)
    kk = 0
    for (k in c(1, seq(m+1, i-m+1))){
      kk = kk + 1
      # regression y on x for k:i
      xx = x[k:i]
      yy = y[k:i]
      xx2 = xx-mean(xx)
      yy.hat = sum(xx2*yy)/sum(xx2^2)*xx2 + mean(yy)
      sigma2.hat = var(yy - yy.hat)
      l.ki = -(i-k)*log(sigma2.hat)/2
      mi[kk] = lambda + Mi[k] + l.ki
      L.ki = exp(l.ki)
      bi = bi + Bi[k] # no need to multiple alpha
      ti = ti + Ti[k]*L.ki # no need to multiple alpha
    }
    Mi[i] = max(mi)
    Bi[i] = bi
    Ti[i] = ti
  }
  ## step 3: final result
  res = list(g2m = 1-exp(-2/n*(Mi[n]-lambda)),
             g2t = 1-(Ti[n]/Bi[n])^{-2/n})
  return(res)
}