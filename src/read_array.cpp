#include "read_array.h"

array_2d readBasePopulation(std::vector<double> basePop) {
	array_2d basePopulation(boost::extents[SEXES][MODEL_AGES]);
	std::copy(basePop.begin(), basePop.end(), basePopulation.data());
	return basePopulation;
}

array_1d readAgeGroupsSpan(std::vector<double> ageGroupsSp) {
	array_1d ageGroupsSpan(boost::extents[AGE_GROUPS]);
	std::copy(ageGroupsSp.begin(), ageGroupsSp.end(), ageGroupsSpan.data());
	return ageGroupsSpan;
}

array_2d readEntrantPrev(SEXP entrantPrev) {
	std::vector<double> eP = Rcpp::as< std::vector<double> >(entrantPrev);
	array_2d entPrev(boost::extents[PROJECTION_YEARS][SEXES]);
	std::copy(eP.begin(), eP.end(), entPrev.data());
	return entPrev;
}

array_1d readVertTransLag(std::vector<double> vertTransLag) {
	array_1d vertTrans(boost::extents[PROJECTION_YEARS]);
	std::copy(vertTransLag.begin(), vertTransLag.end(), vertTrans.data());
	return vertTrans;
}

array_1d readPaedSurveyLag(std::vector<double> paedSurveyLag) {
	array_1d paedSurvLag(boost::extents[PROJECTION_YEARS]);
	std::copy(paedSurveyLag.begin(), paedSurveyLag.end(), paedSurvLag.data());
	return paedSurvLag;
}

array_2d readEntrantPopulation(std::vector<double> entrantPopulation) {
	array_2d entrantPop(boost::extents[PROJECTION_YEARS][SEXES]);
	std::copy(entrantPopulation.begin(), entrantPopulation.end(), entrantPop.data());
	return entrantPop;
}

array_2d readBirthsLag(std::vector<double> birthsLag) {
	array_2d birthLag(boost::extents[PROJECTION_YEARS][SEXES]);
	std::copy(birthsLag.begin(), birthsLag.end(), birthLag.data());
	return birthLag;
}

array_2d readCumulativeSurvey(std::vector<double> cumulativeSurvey) {
	array_2d cumSurv(boost::extents[PROJECTION_YEARS][SEXES]);
	std::copy(cumulativeSurvey.begin(), cumulativeSurvey.end(), cumSurv.data());
	return cumSurv;
}

array_2d readCumulativeNetMigr(std::vector<double> cumulativeNetMigr) {
	array_2d cumNetMigr(boost::extents[PROJECTION_YEARS][SEXES]);
	std::copy(cumulativeNetMigr.begin(), cumulativeNetMigr.end(), cumNetMigr.data());
	return cumNetMigr;
}

array_3d readPaedSurvCd4Dist(std::vector<double> paedSurvCd4Dist) {
	array_3d paedSurvCd4Distrib(boost::extents[PROJECTION_YEARS][SEXES][CD4_STAGES]);
	std::copy(paedSurvCd4Dist.begin(), paedSurvCd4Dist.end(), paedSurvCd4Distrib.data());
	return paedSurvCd4Distrib;
}

array_2d readEntrantArtCoverage(SEXP entrantArtCoverage) {
	array_2d entPrev(boost::extents[PROJECTION_YEARS][SEXES]);
	if (entrantArtCoverage != R_NilValue) {
		std::vector<double> artCov = Rcpp::as< std::vector<double> >(entrantArtCoverage);
		std::copy(artCov.begin(), artCov.end(), entPrev.data());
	} else {
		std::fill_n(entPrev.data(), entPrev.num_elements(), 0);
	}
	return entPrev;
}

array_4d readPaedSurvArtCd4Dist(std::vector<double> paedSurvArtCd4Dist) {
	array_4d paedSurvArtCd4Distrib(boost::extents[PROJECTION_YEARS][SEXES][CD4_STAGES][TREATMENT_STAGES]);
	std::copy(paedSurvArtCd4Dist.begin(), paedSurvArtCd4Dist.end(), paedSurvArtCd4Distrib.data());
	return paedSurvArtCd4Distrib;
}

array_3d readSurvivalRate(std::vector<double> survivalRate) {
	array_3d survRate(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]);
	std::copy(survivalRate.begin(), survivalRate.end(), survRate.data());
	return survRate;
}

array_3d readNetMigration(std::vector<double> netMigration) {
	array_3d netMigr(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]);
	std::copy(netMigration.begin(), netMigration.end(), netMigr.data());
	return netMigr;
}

array_2d readAsfr(std::vector<double> asfRate) {
	array_2d asfr(boost::extents[PROJECTION_YEARS][FERT_AGES]);
	std::copy(asfRate.begin(), asfRate.end(), asfr.data());
	return asfr;
}

array_2d readSexRatioAtBirth(std::vector<double> readSexRatioAtBirth) {
	array_2d sexRatioBirth(boost::extents[PROJECTION_YEARS][SEXES]);
	std::copy(readSexRatioAtBirth.begin(), readSexRatioAtBirth.end(), sexRatioBirth.data());
	return sexRatioBirth;
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
