#ifndef _CPPTEST_READ_ARRAY_H_
#define _CPPTEST_READ_ARRAY_H_

#include <boost/multi_array.hpp>
#include <boost/optional.hpp>
//#include <list>
#include <Rcpp.h>
#include "consts.h"

typedef boost::multi_array<double, 1> array_1d;
typedef boost::multi_array<double, 2> array_2d;
typedef boost::multi_array<double, 3> array_3d;
typedef boost::multi_array<double, 4> array_4d;
typedef boost::multi_array<double, 5> array_5d;

array_2d readEntrantPrev(SEXP entrantPrev);
array_2d readEntrantArtCoverage(SEXP entrantArtCoverage);

array_1d read_1d_array(std::vector<double> arr, int dim1);
array_2d read_2d_array(std::vector<double> arr, int dim1, int dim2);
array_3d read_3d_array(std::vector<double> arr, int dim1, int dim2, int dim3);
array_4d read_4d_array(std::vector<double> arr, int dim1, int dim2, int dim3, int dim4);

std::list<array_3d> read_arrays_list(std::list< std::vector<double> > arr_list);

#endif
