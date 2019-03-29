#include "read_array.h"

array_2d readBasePopulation(std::vector<double> basePop) {
	array_2d basePopulation(boost::extents[SEXES][MODEL_AGES]);
	std::copy(basePop.begin(), basePop.end(), basePopulation.data());
	return basePopulation;
}

array_1d read_1d_array(std::vector<double> arr, int dimensions) {
	array_1d A(boost::extents[dimensions]);
	std::copy(arr.begin(), arr.end(), A.data());
	return A;
}

array_2d read_2d_array(std::vector<double> arr, int dim1, int dim2) {
	array_2d A(boost::extents[dim1][dim2]);
	std::copy(arr.begin(), arr.end(), A.data());
	return A;
}

array_3d read_3d_array(std::vector<double> arr, int dim1, int dim2, int dim3) {
	array_3d A(boost::extents[dim1][dim2][dim3]);
	std::copy(arr.begin(), arr.end(), A.data());
	return A;
}

array_4d read_4d_array(std::vector<double> arr, int dim1, int dim2, int dim3, int dim4) {
	array_4d A(boost::extents[dim1][dim2][dim3][dim4]);
	std::copy(arr.begin(), arr.end(), A.data());
	return A;
}

std::list<array_3d> read_arrays_list(std::list< std::vector<double> > arr_list) {
	std::list<array_3d> arrays;
	for (std::list< std::vector<double> >::iterator it = arr_list.begin(); it != arr_list.end(); it++) {
		arrays.push_back(read_3d_array(*it, 3, 4, 2));
	}
	return arrays;
}
