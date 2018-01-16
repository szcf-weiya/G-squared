#include <iostream>
#include <gsl/gsl_sort.h>

void test1(const double *x, const int n)
{
  size_t *order = (size_t*)malloc(sizeof(size_t)*n);
  gsl_sort_index(order, x, 1, n);
  free(order);
}

void test2(const double *x, const int n)
{
  size_t *order = new size_t[n];
  gsl_sort_index(order, x, 1, n);
  delete [] order;
}

int main()
{
  int  n = 5;
  double x[] = {5, 3, 2, 1, 4};
//  test1(x, n);
  test2(x, n);
  return 0;
}
