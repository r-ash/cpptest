#include "model.h"

// timeSteps here should be SIM_YEARS constant, used var for testing
// [[Rcpp::export]]
std::vector<double> runModel(std::vector<double> basePop, int timeSteps) {
	array_2d basePopulation = readBasePopulation(basePop);
	Model model = Model(basePopulation);
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
