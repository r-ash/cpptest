#include "cpptest.h"
#include "read_array.h"

typedef boost::multi_array<double, 3> array_type;

// [[Rcpp::export]]
double add(double a, double b) {
  return a + b;
}

// [[Rcpp::export]]
double do_double(double x) {
  doubler y = doubler(x);
  return y.run();
}

// [[Rcpp::export]]
std::vector<double> modify_std_vector(std::vector<double> vector) {
  vector[2] = 25;
  return vector;
}
// // [[Rcpp::export]]
// int computeGCD(int a, int b) {
//   return boost::integer::gcd(a, b);
// }

// boost::multi_array<double, 3> modify_3d_array(boost::multi_array<double, 3> array) {
//   std::list< std::vector<double> >::iterator it = list.begin();
//   std::advance(it, 2);
//   it = vec;
//   return list;
// }

// [[Rcpp::export]]
void get_array() {
  array_type A(boost::extents[3][4][2]);
  A[0][0][0] = 1;
  A[0][0][1] = 2;
  A[0][1][0] = 3;
  A[0][1][1] = 4;
  A[0][2][0] = 5;
  A[0][2][1] = 6;
  A[0][3][0] = 7;
  A[0][3][1] = 8;
  A[1][0][0] = 9;
  Rcpp::Rcout << "hi" << A.origin() << "\n";
}


// [[Rcpp::export]]
std::vector<double> push_array(std::vector<double> arr) {
  array_3d A = read_3d_array(arr, 3, 4, 2);
  Rcpp::Rcout << "hi" << A.origin() << "\n";

  std::vector<double> ret(3 * 4 * 2);
  // std::copy would be nicer but the iterators don't want to work for us.
  double * data = A.data();
  for (size_t i = 0; i < A.num_elements(); i++) {
    ret[i] = data[i];
  }
  ret[12] = 55;
  return ret;
}

// [[Rcpp::export]]
std::list< std::vector<double> > push_list_arrays(std::list< std::vector<double> > lst) {
  std::list<array_3d> arrays = read_arrays_list(lst);

  // return back
  std::list< std::vector<double> > ret;
  std::vector<double> ret_array(3 * 4 * 2);
  for (std::list<array_3d>::iterator it = arrays.begin(); it != arrays.end(); it++) {
    double * data = (*it).data();
    for (size_t j = 0; j < (*it).num_elements(); j++) {
      ret_array[j]  = data[j];
      if (j == 12) {
        ret_array[j] = 55;
      }
    }
    ret.push_back(ret_array);
  }
  return ret;
}
