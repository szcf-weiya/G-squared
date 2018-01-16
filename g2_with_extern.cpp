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
#include <omp.h>

int genXY(const int n, const double epsi, const int type, const int resimulate, double *x, double *y, size_t seed)
{
  gsl_rng *r;
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, seed);
  double noise;
  double tmp;
  for (size_t i = 0; i < n; i++)
  {
    x[i] = gsl_rng_uniform(r);
    noise = gsl_ran_gaussian(r, epsi);
    if (type == 0){ //linear
      y[i] = x[i] + noise;
    } else if (type == 1){ // quadradic
      y[i] = pow(x[i], 2) + noise;
    } else if (type == 2){ // cubic
      y[i] = pow(x[i], 3) + noise;
    } else if (type == 3){ // radical
      y[i] = sqrt(x[i]) + noise;
    } else if (type == 4){ // low freq sine
      y[i] = sin(2*PI*x[i]) + noise;
    } else if (type == 5){ // triangle
       tmp = 1 - abs(x[i]);
       if (tmp < 0){
         y[i] = 0 + noise;
       } else {
         y[i] = tmp + noise;
       }
    } else if (type == 6){
      // high freq sine
      y[i] = sin(20*PI*x[i]) + noise;
    } else if (type == 7){
      // step function
      y[i] = floor(x[i]/0.2) + noise;
    } else {
      return 0;
    }
  }
  if (resimulate == 1){
    for (size_t i = 0; i < n; i++){
      x[i] = gsl_rng_uniform(r);
    }
  }
  gsl_rng_free(r);
  return 1;
}

double myexp(const double x)
{
  double tmp = exp(x);
  if (tmp > 1e300) /// max double: 1.79769e+308
    tmp = 1e300;
  return tmp;
}

void add_constant(double *x, const int n, const double constant)
{
  for(size_t i = 0; i < n; i++)
  {
    x[i] += constant;
  }
}

void scale(double *x, const int n, const double f)
{
  for (size_t i = 0; i < n; i++)
  {
    x[i] *= f;
  }
}

double ddot(const double *x, const double *y, const int n)
{
  double res = 0.0;
  for (size_t i = 0; i < n; i++)
  {
    res += x[i]*y[i];
  }
  return res;
}

void sub(double *x, const double *y, const int n)
{
  for (size_t i = 0; i < n; i++)
    x[i] -= y[i];
}

void mycopy(double *des, const double *src, const int n)
{
  for (size_t i = 0; i < n; i++)
    des[i] = src[i];
}

void sortxy(double *x, double *y, const int n)
{
  // sort (xi, yi) by the ascending order of x
  size_t *order = new size_t[n]; // see https://github.com/szcf-weiya/G-squared/issues/3
  gsl_sort_index(order, x, 1, n);
  double *xtmp = new double[n];
  double *ytmp = new double[n];
  mycopy(xtmp, x, n);
  mycopy(ytmp, y, n);
  for (size_t i = 0; i < n; i++){
    x[i] = xtmp[order[i]];
    y[i] = ytmp[order[i]];
  }
  delete [] order;
  delete [] xtmp;
  delete [] ytmp;
}

void normalize(double *y, const int n)
{
  double mu = 0, sd = 0;
  for (size_t i = 0; i < n; i++)
  {
    mu += y[i];
    sd += y[i] * y[i];
  }
  mu = mu/n;
  sd = sqrt((sd - n*mu*mu)/(n-1));
  for (size_t i = 0; i < n; i++)
  {
    y[i] = (y[i] - mu)/sd;
  }
}

void g2(double *x, double *y, const int n, double *g2m, double *g2t)
{
  size_t m = ceil(sqrt(n));
  double lambda = -3*log(n)/2;
  // sort y by the ascending order of x
  sortxy(x, y, n);
  // normalize
  normalize(y, n);
  // cout << gsl_stats_variance(y, 1, n) << endl;// an efficient debug point
  // initialize three sequences
  double *Mi = new double[n];
  double *Bi = new double[n];
  double *Ti = new double[n];
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
    // construct k
    int *k = new int[i-2*m+2];
    k[0] = 0;
    for (size_t ii = 1; ii < i-2*m+2; ii++){
      k[ii] = m -1 + ii;
    }

    for (size_t kk = 0; kk < i-2*m+2; kk++){
      // regression y on x for k:i
      double *xx = new double[i-k[kk]+1];
      double *yy = new double[i-k[kk]+1];
      for (size_t ki = k[kk]; ki <= i; ki++){ // pay attention to the idx
        xx[ki-k[kk]] = x[ki];
        yy[ki-k[kk]] = y[ki];
      }

      add_constant(xx, i-k[kk]+1, -1.*gsl_stats_mean(xx, 1, i-k[kk]+1));
      double sum_xy = ddot(xx, yy, i-k[kk]+1);
      double sum_xx = ddot(xx, xx, i-k[kk]+1);

      scale(xx, i-k[kk]+1, sum_xy/sum_xx);

      add_constant(xx, i-k[kk]+1, gsl_stats_mean(yy, 1, i-k[kk]+1));
      // yyhat - yy

      sub(xx, yy, i-k[kk]+1);

      double sigma2hat = gsl_stats_variance(xx, 1, i-k[kk]+1);
      double lki = -1.0 * (i-k[kk])*log(sigma2hat)/2;
      gsl_vector_set(mi, kk, lambda + Mi[k[kk]] + lki);
      //double Lki = exp(lki);
      double Lki = myexp(lki);
      bi = bi + Bi[k[kk]];
      ti = ti + Ti[k[kk]]*Lki;

      delete [] xx;
      delete [] yy;
    }
    Mi[i] = gsl_vector_max(mi);
    Bi[i] = bi;
    Ti[i] = ti;
    gsl_vector_free(mi);
    delete [] k;
  }
  // final result
  *g2m = 1 - exp(-2.0/n*(Mi[n-1]-lambda));
  *g2t = 1 - pow((Ti[n-1]/Bi[n-1]), -2./n);
  delete [] Mi;
  delete [] Bi;
  delete [] Ti;
}

double calculate_p(const double* v, const int n, const double cut)
{
  int count = 0;
  for(size_t i = 0; i < n; i++)
  {
    if (v[i] > cut)
      count++;
  }
  return count/(n*1.0);
}


extern "C" SEXP StatPower(SEXP R_n, SEXP R_num_noise, SEXP R_num_type, SEXP R_n1, SEXP R_n2)
{
  int n = Rcpp::as<int> (R_n);
  int num_noise = Rcpp::as<int> (R_num_noise);
  int num_type = Rcpp::as<int> (R_num_type);
  int n1 = Rcpp::as<int> (R_n1);
  int n2 = Rcpp::as<int> (R_n2);
  // construct noise level
  double noise;
  Rcpp::NumericMatrix power_cor(num_noise, num_type);
  Rcpp::NumericMatrix power_g2m(num_noise, num_type);
  Rcpp::NumericMatrix power_g2t(num_noise, num_type);
  for (size_t i = 0; i < num_noise; i++)
  {
    noise = sqrt(1.0/(0.2/num_noise*(i+1)) - 1);
    for (size_t j = 0; j < num_type; j++)
    {
      double *val_cor = new double[n1];
      double *val_g2m = new double[n1];
      double *val_g2t = new double[n1];
      # pragma omp parallel for schedule(dynamic)
      for (size_t k = 0; k < n1; k++)
      {
        //Rcpp::DataFrame ds = genXY(n, noise, j, 1);
        //RcppGSL::vector<double> y = ds["y"];
        //RcppGSL::vector<double> x = ds["x"];
        double *x = new double[n];
        double *y = new double[n];
        genXY(n, noise, j, 1, x, y,  k + 10*(n1+n2)*j + 100*num_noise*i);
        // correlation
        val_cor[k] = pow(gsl_stats_correlation(x, 1, y, 1, n), 2);
        // g2
        g2(x, y, n, val_g2m+k, val_g2t+k);
        delete [] x;
        delete [] y;
      }
      // caculate the rejection cutoffs
      gsl_sort(val_cor, 1, n1);
      gsl_sort(val_g2m, 1, n1);
      gsl_sort(val_g2t, 1, n1);
      double cut_cor = gsl_stats_quantile_from_sorted_data(val_cor, 1, n1, 0.95);
      double cut_g2m = gsl_stats_quantile_from_sorted_data(val_g2m, 1, n1, 0.95);
      double cut_g2t = gsl_stats_quantile_from_sorted_data(val_g2t, 1, n1, 0.95);

      delete [] val_cor;
      delete [] val_g2m;
      delete [] val_g2t;

      double *val_cor2 = new double[n2];
      double *val_g2m2 = new double[n2];
      double *val_g2t2 = new double[n2];
      // for alternative hypothesis
      # pragma omp parallel for schedule(dynamic)
      for (size_t k = 0; k < n2; k++)
      {
        double *x2 = new double[n];
        double *y2 = new double[n];
        genXY(n, noise, j, 0, x2, y2, k+n1 + 10*(n1+n2)*j + 100*num_noise*i);
        // correlation
        val_cor2[k] = pow(gsl_stats_correlation(x2, 1, y2, 1, n), 2);
        // g2
        g2(x2, y2, n, val_g2m2+k, val_g2t2+k);
        //cout << val_cor2[k] << " " << val_g2m2[k] << " " << val_g2t2[k] << endl;
        delete [] x2;
        delete [] y2;
      }

      //cout << cut_cor << " " << cut_g2m << " " << cut_g2t << endl;
      power_cor(i, j) = calculate_p(val_cor2, n2, cut_cor);
      power_g2m(i, j) = calculate_p(val_g2m2, n2, cut_g2m);
      power_g2t(i, j) = calculate_p(val_g2t2, n2, cut_g2t);

      delete [] val_cor2;
      delete [] val_g2m2;
      delete [] val_g2t2;
    }
  }
  return (
    Rcpp::List::create(
      Rcpp::Named("cor") = power_cor,
      Rcpp::Named("g2m") = power_g2m,
      Rcpp::Named("g2t") = power_g2t
    )
  );
}
