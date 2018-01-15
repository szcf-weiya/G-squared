library(Rcpp)
Rcpp::sourceCpp("g2.cpp")
x = 1:100
y = 1:100+rnorm(100)
df = data.frame(x=x, y=y)
