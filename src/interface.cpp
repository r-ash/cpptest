#include "interface.h"

// timeSteps here should be SIM_YEARS constant, used var for testing
// [[Rcpp::export]]
std::vector<double> runModel(std::vector<double> basePop, std::vector<double> ageGroupsSpan, int timeArtStart,
                             SEXP entrantPrev, std::vector<double> vertTransLag, std::vector<double> paedSurveyLag,
                             bool populationAdjust, std::vector<double> entrantPop, std::vector<double> birthLag,
                             std::vector<double> cumSurv, std::vector<double> cumNetMigr, double netMigrHivProb,
                             std::vector<double> paedSurvCd4Distrib, SEXP entrantArtCoverage,
                             std::vector<double> paedSurvArtCd4Distrib, std::vector<double> survRate,
                             std::vector<double> netMigr, std::vector<double> asfRate,
                             std::vector<double> sexRatioBirth, int timeSteps) {
	array_2d basePopulation = read_2d_array(basePop, SEXES, MODEL_AGES);
	array_2d entrantPopulation = read_2d_array(entrantPop, PROJECTION_YEARS, SEXES);
	array_2d birthsLag = read_2d_array(birthLag, PROJECTION_YEARS, SEXES);
	array_2d cumulativeSurvey = read_2d_array(cumSurv, PROJECTION_YEARS, SEXES);
	array_2d cumulativeNetMigr = read_2d_array(cumNetMigr, PROJECTION_YEARS, SEXES);
	array_3d paedSurvCd4Dist = read_3d_array(paedSurvCd4Distrib, PROJECTION_YEARS, SEXES, CD4_STAGES);
	array_2d entrantArtCov = readEntrantArtCoverage(entrantArtCoverage);
	array_4d paedSurvArtCd4Dist = read_4d_array(paedSurvArtCd4Distrib, PROJECTION_YEARS, SEXES, CD4_STAGES, TREATMENT_STAGES);
	array_3d survivalRate = read_3d_array(survRate, PROJECTION_YEARS, SEXES, MODEL_AGES);
	array_3d netMigration = read_3d_array(netMigr, PROJECTION_YEARS, SEXES, MODEL_AGES);
	array_2d asfr = read_2d_array(asfRate, PROJECTION_YEARS, FERT_AGES);
	array_2d sexRatioAtBirth = read_2d_array(sexRatioBirth, PROJECTION_YEARS, SEXES);

	Model model = Model(basePopulation, ageGroupsSpan, vertTransLag, paedSurveyLag, populationAdjust,
	                    entrantPopulation, birthsLag, cumulativeSurvey, cumulativeNetMigr, netMigrHivProb,
	                    paedSurvCd4Dist, entrantArtCov, paedSurvArtCd4Dist, survivalRate, netMigration,
	                    asfr, sexRatioAtBirth, timeArtStart);
	if (entrantPrev != R_NilValue) {
		array_2d entPrev = readEntrantPrev(entrantPrev);
		model.setEntrantPrev(entPrev);
	}

	for (int t = 1; t <= timeSteps; t++) {
		//Rcpp::Rcout << "Looping over time step " << t << "\n";
		model.agePopulation(t);
		model.deathsAndMigration(t);
		model.fertility(t);
	}

	// Convert back to view for R for testing
	std::vector<double> ret_array(PROJECTION_YEARS * DISEASE_STATUS * SEXES * MODEL_AGES);
	double * data = model.getPopulation().data();
	for (size_t i = 0; i < (model.getPopulation()).num_elements(); i++) {
		ret_array[i] = data[i];
	}
	return ret_array;
}
