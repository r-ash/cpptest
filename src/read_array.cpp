#include "read_array.h"

array_3d read_array(std::vector<double> arr) {
  array_3d A(boost::extents[3][4][2]);
  std::copy(arr.begin(), arr.end(), A.data());
  return A;
}

std::list<array_3d> read_arrays_list(std::list< std::vector<double> > arr_list) {
  std::list<array_3d> arrays;
  for (std::list< std::vector<double> >::iterator it = arr_list.begin(); it != arr_list.end(); it++) {
    arrays.push_back(read_array(*it));
  }
  return arrays;
}
