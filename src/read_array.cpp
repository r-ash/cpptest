#include "read_array.h"

array_1d read_1d_array(std::vector<double> arr, int dimensions) {
  array_1d A(boost::extents[dimensions]);
  std::copy(arr.begin(), arr.end(), A.data());
  return A;
}

array_2d read_2d_array(std::vector<double> arr, std::vector<int> dimensions) {
  array_2d A(boost::extents[dimensions[0]][dimensions[1]]);
  std::copy(arr.begin(), arr.end(), A.data());
  return A;
}

array_3d read_3d_array(std::vector<double> arr, std::vector<int> dimensions) {
  array_3d A(boost::extents[dimensions[0]][dimensions[1]][dimensions[2]]);
  std::copy(arr.begin(), arr.end(), A.data());
  return A;
}

array_4d read_4d_array(std::vector<double> arr, std::vector<int> dimensions) {
  array_4d A(boost::extents[dimensions[0]][dimensions[1]][dimensions[2]][dimensions[3]]);
  std::copy(arr.begin(), arr.end(), A.data());
  return A;
}

std::list<array_3d> read_arrays_list(std::list< std::vector<double> > arr_list) {
  std::list<array_3d> arrays;
  std::vector<int> dims = {3,4,2};
  for (std::list< std::vector<double> >::iterator it = arr_list.begin(); it != arr_list.end(); it++) {
    arrays.push_back(read_3d_array(*it, dims));
  }
  return arrays;
}
