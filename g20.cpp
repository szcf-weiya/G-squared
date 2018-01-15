#include <iostream>
using namespace std;
#include <math.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>

// [[Rcpp::export]]
Rcpp::List genXY(const int n, const double epsi, const int type, const int resimulate)
{
  gsl_rng *r;
  RcppGSL::vector<double> x(n), y(n);
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);
  double noise = gsl_ran_gaussian(r, epsi);
  double tmp, xtmp;
  for (size_t i = 0; i < n; i++)
  {
    xtmp = gsl_rng_uniform(r);
    gsl_vector_set(x, i, xtmp);
    if (type == 1){ //linear
      gsl_vector_set(y, i, xtmp+noise);
    } else if (type == 2){ // quadradic
      gsl_vector_set(y, i, pow(xtmp, 2)+noise);
    } else if (type == 3){ // cubic
      gsl_vector_set(y, i, pow(xtmp, 3)+noise);
    } else if (type == 4){ // radical
      gsl_vector_set(y, i, sqrt(xtmp) + noise);
    } else if (type == 5){ // low freq sine
      gsl_vector_set(y, i, sin(2*PI*xtmp));
    } else if (type == 6){ // triangle
       tmp = 1-abs(xtmp);
       if (tmp < 0){
         gsl_vector_set(y, i, 0);
       } else {
         gsl_vector_set(y, i, tmp);
       }
    } else if (type == 7){
      // high freq sine
      gsl_vector_set(y, i, sin(20*PI*xtmp));
    } else if (type == 8){
      // step function
      gsl_vector_set(y, i, floor(xtmp/0.2));
    } else {
      return R_NilValue;
    }
  }
  if (resimulate == 1){
    for (size_t i = 0; i < n; i++){
      gsl_vector_set(x, i, gsl_rng_uniform(r));
    }
  }
  Rcpp::DataFrame res = Rcpp::DataFrame::create(Rcpp::Named("x") = x, Rcpp::Named("y") = y);
  x.free();
  y.free();
  gsl_rng_free(r);
  return(res);
}

// [[Rcpp::export]]
Rcpp::List g2(Rcpp::DataFrame ds)
{
  Rcpp::DataFrame D(ds);
  RcppGSL::vector<double> y = D["y"];
  RcppGSL::vector<double> x = D["x"];
  size_t n = x->size;
  size_t m = ceil(sqrt(n));
  double lambda = -3*log(n)/2;
  cout << lambda << endl;
  size_t *order = (size_t*)malloc(sizeof(size_t)*n);
  gsl_sort_index(order, x->data, 1, n);
  gsl_vector *tmp = gsl_vector_alloc(n);
  gsl_vector_memcpy(tmp, x);
  // sort (xi, yi) by the ascending order of x
  for (size_t i = 0; i < n; i++){
    gsl_vector_set(x, i, gsl_vector_get(tmp, order[i]));
    cout << gsl_vector_get(x, i) << " " << endl;
  }
  gsl_vector_memcpy(tmp, y);
  for (size_t i = 0; i < n; i++){
    gsl_vector_set(y, i, gsl_vector_get(tmp, order[i]));
  }
  // normalize
  gsl_vector_add_constant(y, -1.0*gsl_stats_mean(y->data, 1, n)); // pay attention to the sign
  gsl_vector_scale(y, 1./gsl_stats_sd(y->data, 1, n));
  //double sum_yy;
  //gsl_blas_ddot(y, y, &sum_yy);
  //sum_yy = gsl_blas_dnrm2(y);
  //gsl_vector_scale(y, sqrt(n)/sum_yy);
  // initialize three sequences
  //RcppGSL::vector<double> Mi(n);
  gsl_vector *Mi = gsl_vector_calloc(n);
  RcppGSL::vector<double> Bi(n);
  RcppGSL::vector<double> Ti(n);
  Mi[0] = 0;
  Bi[0] = 1;
  Ti[0] = 1;
  double bi, ti;
  for(size_t i = m - 1; i < n; i++){
    bi = 0;
    ti = 0;
    if (i < 2*m)
    {
      Mi[i] = Mi[0];
      Bi[i] = Bi[0];
      Ti[i] = Ti[0];
      continue;
    }
    gsl_vector *mi = gsl_vector_calloc(i-2*m+2);
    gsl_vector_add_constant(mi, -10000000);
    // construct k
    RcppGSL::vector<int> k(i-2*m+2);
    k[0] = 0;
    for (size_t ii = 1; ii < i-2*m+2; ii++){
      k[ii] = m -1 + ii;
    }

    for (size_t kk = 0; kk < i-2*m+2; kk++){
      // regression y on x for k:i
      RcppGSL::vector<double> xx(i-k[kk]+1);
      RcppGSL::vector<double> yy(i-k[kk]+1);
      RcppGSL::vector<double> yyhat(i-k[kk]+1);
      for (size_t ki = k[kk]-1; ki < i; ki++){
        xx[ki-k[kk]+1] = gsl_vector_get(x, ki);
        yy[ki-k[kk]+1] = gsl_vector_get(y, ki);
      }
      gsl_vector_add_constant(xx, -1.*gsl_stats_mean(xx->data, 1, i-k[kk]+1));
      double sum_xy, sum_xx;
      gsl_blas_ddot(xx, yy, &sum_xy);
      gsl_blas_ddot(xx, xx, &sum_xx);
      //sum_xx = pow(gsl_blas_dnrm2(xx), 2);
      gsl_vector_memcpy(yyhat, xx);
      gsl_vector_scale(yyhat, sum_xy/sum_xx);
      gsl_vector_add_constant(yyhat, gsl_stats_mean(yy->data, 1, i-k[kk]+1));
      gsl_vector_sub(yyhat, yy);
      double sigma2hat;
      sigma2hat = gsl_stats_variance(yyhat->data, 1, i-k[kk]+1);
      double lki = -(i-k[kk])*log(sigma2hat)/2;
      gsl_vector_set(mi, kk, lambda + Mi[k[kk]] + lki);
      //cout << "mi = " << gsl_vector_get(mi, kk) << endl;
      cout << lambda << " " << Mi[k[kk]] << endl;
      double Lki = exp(lki);
      bi = bi + Bi[k[kk]];
      ti = ti + Ti[k[kk]]*Lki;
      xx.free();
      yy.free();
      yyhat.free();
    }
    Mi[i] = gsl_vector_max(mi);
    Bi[i] = bi;
    Ti[i] = ti;
    gsl_vector_free(mi);
    k.free();
  }
  // step final result
  return(Rcpp::List::create(
    Rcpp::Named("g2m") = 1-exp(-2.0/n*(Mi[n-1]-lambda)),
    Rcpp::Named("g2t") = 1-pow((Ti[n-1]/Bi[n-1]), -2./n)
  ));

}
