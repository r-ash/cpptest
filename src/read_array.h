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

array_2d readBasePopulation(std::vector<double> arr);
array_1d readAgeGroupsSpan(std::vector<double> ageGroupsSp);
array_2d readEntrantPrev(SEXP entrantPrev);
array_1d readVertTransLag(std::vector<double> vertTransLag);
array_1d readPaedSurveyLag(std::vector<double> paedSurveyLag);
array_2d readEntrantPopulation(std::vector<double> entrantPopulation);
array_2d readBirthsLag(std::vector<double> birthsLag);
array_2d readCumulativeSurvey(std::vector<double> cumulativeSurvey);
array_2d readCumulativeNetMigr(std::vector<double> cumulativeNetMigr);
array_3d readPaedSurvCd4Dist(std::vector<double> paedSurvCd4Dist);
array_2d readEntrantArtCoverage(SEXP entrantArtCoverage);
array_4d readPaedSurvArtCd4Dist(std::vector<double> paedSurvArtCd4Dist);
array_3d readSurvivalRate(std::vector<double> survivalRate);
array_3d readNetMigration(std::vector<double> netMigration);
array_2d readAsfr(std::vector<double> asfRate);
array_2d readSexRatioAtBirth(std::vector<double> readSexRatioAtBirth);

array_1d read_1d_array(std::vector<double> arr, int dim1);
array_2d read_2d_array(std::vector<double> arr, int dim1, int dim2);
array_3d read_3d_array(std::vector<double> arr, int dim1, int dim2, int dim3);
array_4d read_4d_array(std::vector<double> arr, int dim1, int dim2, int dim3, int dim4);

std::list<array_3d> read_arrays_list(std::list< std::vector<double> > arr_list);

#endif
