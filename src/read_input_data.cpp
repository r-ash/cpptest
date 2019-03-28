#include "read_input_data.h"

#define AGE_START 15

#define NG 2
#define pAG 66
#define pDS 2

#define pIDX_FERT 0
#define pAG_FERT 35
#define pIDX_15TO49 0
#define pAG_15TO49  35
#define pIDX_15PLUS 0
#define pAG_15PLUS  66

#define hAG 9
#define hDS 7
#define hTS 3

#define hIDX_FERT 0
#define hAG_FERT 8
#define hIDX_15TO49 0
#define hAG_15TO49  8
#define hIDX_15PLUS 0
#define hAG_15PLUS  9

#define hIDX_CD4_350 2

#define MALE 0
#define FEMALE 1

#define HIVN 0
#define HIVP 1

#define ART0MOS 0
#define ART6MOS 1
#define ART1YR 2

#define ART_STAGE_PROG_RATE 2.0 // HARD CODED: ART stage progression rate

#define EPP_RSPLINE 0
#define EPP_RTREND 1
#define EPP_DIRECTINCID 2  // annual direct incidence inputs (as Spectrum)

#define INCIDMOD_EPPSPEC 0
#define INCIDMOD_TRANSM 1

#define INCIDPOP_15TO49 0 // age range corresponding to incidence input
#define INCIDPOP_15PLUS 1

std::list<array_1d> read_1d_arrays(std::list< std::vector<double> > arr_list) {
  std::list<array_1d> arrays;
  for (std::list< std::vector<double> >::iterator it = arr_list.begin(); it != arr_list.end(); it++) {
    arrays.push_back(read_1d_array(*it, arr_list.size()));
  }
  return arrays;
}

// [[Rcpp::export]]
void test_1d_arrays(std::list< std::vector<double> > arr_list) {
	std::list<array_1d> arrays = read_1d_arrays(arr_list);
	Rcpp::Rcout << "list of size" << arrays.size() << "\n";
	return;
}

std::list<array_2d> read_2d_arrays(std::list< std::vector<double> > arr_list, std::list< std::vector<int> > dimensions) {
  std::list<array_2d> arrays;
  std::list< std::vector<double> >::iterator it1 = arr_list.begin();
  std::list< std::vector<int> >::iterator it2 = dimensions.begin();
  while(it1 != arr_list.end() && it2 != dimensions.end()) {
    arrays.push_back(read_2d_array(*it1, *it2));
    it1++;
    it2++;
  }
  return arrays;
}

// [[Rcpp::export]]
void test_2d_arrays(std::list< std::vector<double> > arr_list, std::list< std::vector<int> > dimensions) {
	std::list<array_2d> arrays = read_2d_arrays(arr_list, dimensions);
	Rcpp::Rcout << "list of size" << arrays.size() << "\n";
	return;
}

std::list<array_3d> read_3d_arrays(std::list< std::vector<double> > arr_list, std::list< std::vector<int> > dimensions) {
  std::list<array_3d> arrays;
  std::list< std::vector<double> >::iterator it1 = arr_list.begin();
  std::list< std::vector<int> >::iterator it2 = dimensions.begin();
  while(it1 != arr_list.end() && it2 != dimensions.end()) {
    arrays.push_back(read_3d_array(*it1, *it2));
    it1++;
    it2++;
  }
  return arrays;
}

//void test_3d_arrays(std::list< std::vector<double> > arr_list, std::list< std::vector<int> > dimensions) {
//	std::list<array_3d> array = read_3d_arrays(arr_list);
//}

std::list<array_4d> read_4d_arrays(std::list< std::vector<double> > arr_list, std::list< std::vector<int> > dimensions) {
  std::list<array_4d> arrays;
  std::list< std::vector<double> >::iterator it1 = arr_list.begin();
  std::list< std::vector<int> >::iterator it2 = dimensions.begin();
  while(it1 != arr_list.end() && it2 != dimensions.end()) {
    arrays.push_back(read_4d_array(*it1, *it2));
    it1++;
    it2++;
  }
  return arrays;
}

//void test_4d_arrays(std::list< std::vector<double> > arr_list, std::list< std::vector<int> > dimensions) {
//	std::list<array_4d> array = read_4d_arrays(arr_list);
//}