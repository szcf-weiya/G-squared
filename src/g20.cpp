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
    gsl_rng_set(r, seed*100000);
    for (size_t i = 0; i < n; i++){
      x[i] = gsl_rng_uniform(r);
    }
  }
  gsl_rng_free(r);
  return 1;
}

/*
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
*/

double myexp(const double x)
{
  double tmp = exp(x);
  if (tmp > 1e300) /// max double: 1.79769e+308
    tmp = 1e300;
  return tmp;
}


// [[Rcpp::export]]
Rcpp::List g2cpp(Rcpp::DataFrame ds)
{
  Rcpp::DataFrame D(ds);
  RcppGSL::vector<double> y = D["y"];
  RcppGSL::vector<double> x = D["x"];
  size_t n = x->size;
  size_t m = ceil(sqrt(n));
  double lambda = -3*log(n)/2;
  size_t *order = (size_t*)malloc(sizeof(size_t)*n);
  gsl_sort_index(order, x->data, 1, n);
  gsl_vector *tmp = gsl_vector_alloc(n);
  gsl_vector_memcpy(tmp, x);
  // sort (xi, yi) by the ascending order of x
  for (size_t i = 0; i < n; i++){
    gsl_vector_set(x, i, gsl_vector_get(tmp, order[i]));
  }
  gsl_vector_memcpy(tmp, y);
  for (size_t i = 0; i < n; i++){
    gsl_vector_set(y, i, gsl_vector_get(tmp, order[i]));
  }
  free(order);
  // normalize
  gsl_vector_add_constant(y, -1.0*gsl_stats_mean(y->data, 1, n)); // pay attention to the sign
  gsl_vector_scale(y, 1./gsl_stats_sd(y->data, 1, n));
  //double sum_yy;
  //gsl_blas_ddot(y, y, &sum_yy);
  //sum_yy = gsl_blas_dnrm2(y);
  //gsl_vector_scale(y, sqrt(n)/sum_yy);
  // initialize three sequences
  RcppGSL::vector<double> Mi(n);
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
      for (size_t ki = k[kk]; ki <= i; ki++){ // pay attention to the idx
        xx[ki-k[kk]] = gsl_vector_get(x, ki);
        yy[ki-k[kk]] = gsl_vector_get(y, ki);
      }
      gsl_vector_add_constant(xx, -1.*gsl_stats_mean(xx->data, 1, i-k[kk]+1));
      double sum_xy, sum_xx;
      gsl_blas_ddot(xx, yy, &sum_xy);
      gsl_blas_ddot(xx, xx, &sum_xx);
      gsl_vector_memcpy(yyhat, xx);
      gsl_vector_scale(yyhat, sum_xy/sum_xx);
      gsl_vector_add_constant(yyhat, gsl_stats_mean(yy->data, 1, i-k[kk]+1));
      //cout << gsl_stats_mean(yy->data, 1, i-k[kk]+1) << endl;

      //for (size_t j = 0; j < i-k[kk]+1; j++){
      //  cout << xx[j] << " " << yy[j] << " " << yyhat[j] << endl;
      //}
      gsl_vector_sub(yyhat, yy);
      double sigma2hat;
      sigma2hat = gsl_stats_variance(yyhat->data, 1, i-k[kk]+1);
      //cout << "sigma2hat" << sigma2hat << endl;
      double lki = -1.0* (i-k[kk])*log(sigma2hat)/2;
      gsl_vector_set(mi, kk, lambda + Mi[k[kk]] + lki);
      //double Lki = exp(lki);
      double Lki = myexp(lki);
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
  // step 3:  final result
  return(Rcpp::List::create(
    Rcpp::Named("g2m") = 1-exp(-2.0/n*(Mi[n-1]-lambda)),
    Rcpp::Named("g2t") = 1-pow((Ti[n-1]/Bi[n-1]), -2./n)
  ));
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
  cout << count << endl;
  return count*1.0/n;
}

// [[Rcpp::export]]
Rcpp::List test(Rcpp::DataFrame ds)
{
  Rcpp::DataFrame D(ds);
  RcppGSL::vector<double> y = D["y"];
  RcppGSL::vector<double> x = D["x"];
  int n = 20;
  RcppGSL::vector<double> yy(n);
  RcppGSL::vector<double> xx(n);
  double g2m, g2t;
  g2(x->data, y->data, x->size, &g2m, &g2t);
  genXY(n, 1, 1, 0, xx->data, yy->data, 1234);
  return (
    Rcpp::List::create(
      Rcpp::Named("g2m") = g2m,
      Rcpp::Named("g2t") = g2t,
      Rcpp::Named("x") = xx,
      Rcpp::Named("y") = yy
    )
  );
}


// [[Rcpp::export]]
Rcpp::List StatPower(const int n, const int num_noise, const int num_type, const int n1, const int n2)
{
  // construct noise level
  double noise;
  Rcpp::NumericMatrix power_cor(num_noise, num_type);
  Rcpp::NumericMatrix power_g2m(num_noise, num_type);
  Rcpp::NumericMatrix power_g2t(num_noise, num_type);
  for (size_t i = 0; i < num_noise; i++)
  {
    //noise = sqrt(1.0/(0.001/num_noise*(i+1)) - 1);
    noise = exp(1+num_noise);
    for (size_t j = 0; j < num_type; j++)
    {
      double *val_cor = new double[n1];
      double *val_g2m = new double[n1];
      double *val_g2t = new double[n1];

      for (size_t k = 0; k < n1; k++)
      {
        //Rcpp::DataFrame ds = genXY(n, noise, j, 1);
        //RcppGSL::vector<double> y = ds["y"];
        //RcppGSL::vector<double> x = ds["x"];
        double *x = new double[n];
        double *y = new double[n];
        genXY(n, noise, j, 1, x, y, k + 10*(n1+n2)*j + 100*num_noise*i);
        // correlation
        val_cor[k] = pow(gsl_stats_correlation(x, 1, y, 1, n), 2);
        // g2
        g2(x, y, n, val_g2m+k, val_g2t+k);
        if (k == n1/2)
        {
          cout << x[n/3] << " " << x[n/2] << " ;" << y[n/3] << " " << y[n/2] <<endl;
          cout << val_cor[k] << " "<< val_g2m[k] << " "<< val_g2t[k] << endl;
        }

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
