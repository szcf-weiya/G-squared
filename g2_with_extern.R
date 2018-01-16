library(Rcpp)
dyn.load("g2_with_extern.so")
StatPower <- function(N, num.noise, num.type, n1, n2){
  .Call("StatPower", N, num.noise, num.type, n1, n2)
}
StatPower(225, 20, 8, 1000, 1000)
