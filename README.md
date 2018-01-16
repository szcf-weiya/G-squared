# G-squared
Implementation of Generalized R-squared

## Compare the results by R and Rcpp

The code is as follows:
```r
library(Rcpp)
N = 1000
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
```

We would get the following results

![](compare_res.png)

Our original goal is to speed up the program by Rcpp. The following results demonstrate the Rcpp program works!


![](compare_time.png)

## Reference
1. Wang, X., Jiang, B., & Liu, J. S. (2017). Generalized R-squared for detecting dependence. Biometrika,​ ​ 104(1),​ ​ 129-139.
2. Simon N. and Tibshirani R. :http://statweb.stanford.edu/~tibs/reshef/script.R
