#ifndef _CPPTEST_READ_ARRAY_H_
#define _CPPTEST_READ_ARRAY_H_

#include <boost/multi_array.hpp>
#include <Rcpp.h>

typedef boost::multi_array<double, 3> array_3d;
array_3d read_array(std::vector<double> arr);
std::list<array_3d> read_arrays_list(std::list< std::vector<double> > arr_list);

#endif
