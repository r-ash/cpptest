#ifndef _CPPTEST_READ_ARRAY_H_
#define _CPPTEST_READ_ARRAY_H_

#include <boost/multi_array.hpp>
#include <Rcpp.h>

typedef boost::multi_array<double, 1> array_1d;
typedef boost::multi_array<double, 2> array_2d;
typedef boost::multi_array<double, 3> array_3d;
typedef boost::multi_array<double, 4> array_4d;


array_1d read_1d_array(std::vector<double> arr, int dimensions);
array_2d read_2d_array(std::vector<double> arr, std::vector<int> dimensions);
array_3d read_3d_array(std::vector<double> arr, std::vector<int> dimensions);
array_4d read_4d_array(std::vector<double> arr, std::vector<int> dimensions);

std::list<array_3d> read_arrays_list(std::list< std::vector<double> > arr_list);

#endif
