library(Rcpp)
N = 225
Rcpp::sourceCpp("g20.cpp")
x = 1:N
gap = array(NA, c(2, 20))
for (i in 1:20)
{
  y = sin(1:N) + rnorm(N, sd = i)
  df = data.frame(x=x, y=y)
  res1 = g2cpp(df)
  res2 = g2(x, y)
  gap[1, i] = res1$g2m - res2$g2m
  gap[2, i] = res1$g2t - res2$g2t
}
gap