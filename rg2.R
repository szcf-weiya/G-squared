library(Rcpp)
Rcpp::sourceCpp("g20.cpp")


## statpower
N = 100
num.noise = 6 # 30
num.type = 2 # 8
n1 = 200
n2 = 200
system.time(
  {
    res = StatPower(N, num.noise, num.type, n1, n2)
  }
)

## test 
x = 1:N
y = 1:N + rnorm(N, sd = 50)
df = data.frame(x=x, y=y)
test(df)
g2(x,y)

## compare results
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

## compare time 
system.time(
  {
    for (i in 1:20)
      g2(x,y)
  }
)
system.time(
  {
    for (i in 1:20)
      g2cpp(df)
  }
)
system.time(
  {
    for (i in 1:20)
      test(df)
  }
)