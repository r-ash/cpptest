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
                             std::vector<double> sexRatioBirth, int hivStepsPerYear, std::vector<double> cd4Prog,
                             std::vector<double> cd4InitDist, std::vector<double> cd4Mort, std::vector<double> incrrAges,
                             double relinfectArt, SEXP iota, std::vector<double> incrrSex, SEXP incidMod, int eppMod,
                             int scaleCd4Mort, std::vector<double> projSteps, SEXP tsEpidemicStart, std::vector<int> artCd4EligId,
                             std::vector<double> specPopPercentElig, std::vector<double> pregnantWomenArtElig,
                             double who34PercentElig, SEXP rSplineRVec, SEXP rTrendBeta, SEXP rTrendTStab, SEXP rTrendR0,
                             int timeSteps) {
	Rcpp::Rcout << "Got args \n";
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
	array_3d cd4Progression = read_3d_array(cd4Prog, SEXES, AGE_GROUPS, CD4_STAGES - 1);
	array_3d cd4InitialDist = read_3d_array(cd4InitDist, SEXES, AGE_GROUPS, CD4_STAGES);
	array_3d cd4Mortality = read_3d_array(cd4Mort, SEXES, AGE_GROUPS, CD4_STAGES);
	array_3d incrrAge = read_3d_array(incrrAges, PROJECTION_YEARS, SEXES, MODEL_AGES);
	Rcpp::Rcout << "initialising model \n";

	ArtData art = ArtData(entrantArtCov, paedSurvArtCd4Dist, artCd4EligId,
	                      specPopPercentElig, pregnantWomenArtElig, who34PercentElig);
	State state = State(basePopulation, art.population);
	Model model = Model(state, art, ageGroupsSpan, vertTransLag, paedSurveyLag, populationAdjust,
	                    entrantPopulation, birthsLag, cumulativeSurvey, cumulativeNetMigr, netMigrHivProb,
	                    paedSurvCd4Dist, survivalRate, netMigration,
	                    asfr, sexRatioAtBirth, hivStepsPerYear, cd4Progression, cd4InitialDist,
	                    cd4Mortality, incrrAge, timeArtStart, relinfectArt, eppMod, scaleCd4Mort,
	                    projSteps);
	if (entrantPrev != R_NilValue) {
		array_2d entPrev = readEntrantPrev(entrantPrev);
		model.setEntrantPrev(entPrev);
	}
	if (incidMod != R_NilValue) {
		model.setIncrrSex(incrrSex);
	}

	if (eppMod != EPP_DIRECTINCID) {
		model.intialiseNonDirectIncid(Rcpp::as<double>(iota), Rcpp::as<double>(tsEpidemicStart));
	}
	if (eppMod == EPP_RSPLINE) {
		model.initialiseRSpline(Rcpp::as < std::vector<double> >(rSplineRVec));
	} else if (eppMod == EPP_RTREND) {
		model.initialiseRTrend(Rcpp::as< std::vector<double> >(rTrendBeta), Rcpp::as<double>(rTrendTStab),
		                       Rcpp::as<double>(rTrendR0));
	}

	for (int t = 1; t <= timeSteps; t++) {
		Rcpp::Rcout << "Looping over time step " << t << "\n";
		model.agePopulation(t);
		model.deathsAndMigration(t);
		model.fertility(t);
		model.updateModelState(t);
	}

	// Convert back to view for R for testing
	std::vector<double> ret_array(PROJECTION_YEARS * DISEASE_STATUS * SEXES * MODEL_AGES);
	double * data = model.state.outputPopulation.data();
	for (size_t i = 0; i < model.state.outputPopulation.num_elements(); i++) {
		ret_array[i] = data[i];
	}
	return ret_array;
}
