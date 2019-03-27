#include "cpptest.h"
#include <cstddef>

// [[Rcpp::export]]
double add(double a, double b) {
  return a + b;
}

// [[Rcpp::export]]
double do_double(double x) {
  doubler y = doubler(x);
  return y.run();
}
