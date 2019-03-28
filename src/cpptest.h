#ifndef _CPPTEST_CPPTEST_H_
#define _CPPTEST_CPPTEST_H_

#include "read_array.h"
#include <Rcpp.h>
#include <boost/multi_array.hpp>

double add(double a, double b);

class doubler {
private:
  double number;
public:
  // doubler(double x) {
  //   number = x;
  // }
  doubler(double x) : number(x) {}
  double run() {
    number = number * 2;
    return number;
  }
};

std::vector<double> push_array(std::vector<double> arr);

std::list< std::vector<double> > push_list_arrays(std::list< std::vector<double> > lst);

#endif
