library(Rcpp)
dyn.load("g2_with_extern.so")
StatPower <- function(N, num.noise, num.type, n1, n2){
  .Call("StatPower", N, num.noise, num.type, n1, n2)
}
system.time(
  {
    res = StatPower(225, 30, 8, 1000, 100)
  }
)
#StatPower(100, 4, 8, 100, 100)
#StatPower(225, 1, 8, 500, 100)


## save 
write.csv(res$cor, "res_cor.csv")
write.csv(res$g2m, "res_g2m.csv")
write.csv(res$g2t, "res_g2t.csv")