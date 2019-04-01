#include "model.h"

// timeSteps here should be SIM_YEARS constant, used var for testing
// [[Rcpp::export]]
std::vector<double> runModel(std::vector<double> basePop, std::vector<double> ageGroupsSpan, int timeArtStart,
                             std::vector<double> entrantPrev, int timeSteps) {
	array_2d basePopulation = readBasePopulation(basePop);
	array_1d ageGroupsSp = readAgeGroupsSpan(ageGroupsSpan);

	Model model = Model(basePopulation, ageGroupsSp, timeArtStart);
	if (entrantPrev.size() > 0) {
		array_2d entPrev = readEntrantPrev(entrantPrev);
		model.setEntrantPrev(entPrev);
	}

	for (int t = 1; t < timeSteps; t++) {
		Rcpp::Rcout << "Looping over time step " << t << "\n";
		model.agePopulation(t);
	}

	// Convert back to view for R for testing
	std::vector<double> ret_array(PROJECTION_YEARS * DISEASE_STATUS * SEXES * MODEL_AGES);
	double * data = model.getPopulation().data();
	for (size_t i = 0; i < (model.getPopulation()).num_elements(); i++) {
		ret_array[i] = data[i];
	}
	return ret_array;
}

// [[Rcpp::export]]
std::vector<double> runModel(std::vector<double> basePop, std::vector<double> ageGroupsSpan, int timeArtStart,
                             int timeSteps) {
	std::vector<double> entrantPrev = {};
	return runModel(basePop, ageGroupsSpan, timeArtStart, entrantPrev, timeSteps);
}

// Notes before putting down there are some issues here with compilation. The problems are arising because entrantPrev may be null in the R code
// overloading methods here doesn't seem to work but also trying to use boost::optional<> types doesn't work as it isn't an SEXP type so Rcpp
// cannot convert it for being alled from R. Using Rcpp::Nullable type also doesn't work because this can't wrap an std::vector.