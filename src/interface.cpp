#include "interface.h"

// timeSteps here should be SIM_YEARS constant, used var for testing
// [[Rcpp::export]]
std::vector<double> runModel(std::vector<double> basePop, std::vector<double> ageGroupsSpan, int timeArtStart,
                             SEXP entrantPrev, std::vector<double> vertTransLag, std::vector<double> paedSurveyLag,
                             bool populationAdjust, std::vector<double> entrantPop, std::vector<double> birthLag,
                             std::vector<double> cumSurv, std::vector<double> cumNetMigr, double netMigrHivProb,
                             std::vector<double> paedSurvCd4Distrib, SEXP entrantArtCoverage,
                             std::vector<double> paedSurvArtCd4Distrib, int timeSteps) {
	array_2d basePopulation = readBasePopulation(basePop);
	array_1d ageGroupsSp = readAgeGroupsSpan(ageGroupsSpan);
	array_1d vertTrans = readVertTransLag(vertTransLag);
	array_1d paedSurvLag = readPaedSurveyLag(paedSurveyLag);
	array_2d entrantPopulation = readEntrantPopulation(entrantPop);
	array_2d birthsLag = readBirthsLag(birthLag);
	array_2d cumulativeSurvey = readCumulativeSurvey(cumSurv);
	array_2d cumulativeNetMigr = readCumulativeNetMigr(cumNetMigr);
	array_3d paedSurvCd4Dist = readPaedSurvCd4Dist(paedSurvCd4Distrib);
	array_2d entrantArtCov = readEntrantArtCoverage(entrantArtCoverage);
	array_4d paedSurvArtCd4Dist = readPaedSurvArtCd4Dist(paedSurvArtCd4Distrib);

	Model model = Model(basePopulation, ageGroupsSp, vertTrans, paedSurvLag, populationAdjust,
	                    entrantPopulation, birthsLag, cumulativeSurvey, cumulativeNetMigr, netMigrHivProb,
	                    paedSurvCd4Dist, entrantArtCov, paedSurvArtCd4Dist, timeArtStart);
	if (entrantPrev != R_NilValue) {
		array_2d entPrev = readEntrantPrev(entrantPrev);
		model.setEntrantPrev(entPrev);
	}

	for (int t = 1; t < timeSteps; t++) {
		//Rcpp::Rcout << "Looping over time step " << t << "\n";
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
