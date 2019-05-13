#ifndef _CPPTEST_MODEL_H_
#define _CPPTEST_MODEL_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "read_array.h"
#include "consts.h"
#include "state.h"
#include "ArtData.h"

class Model {
private:
	std::vector<double> ageGroupsSpan;
	int timeArtStart;
	bool useEntrantPrev;
	array_2d entrantPrev;
	std::vector<double> previousPregnancyLag;
	std::vector<double> vertTransLag;
	std::vector<double> paedSurveyLag;
	bool populationAdjust;
	array_2d entrantPopulation;
	array_2d birthsLag;
	array_2d cumulativeSurvey;
	array_2d cumulativeNetMigr;
	double netMigrHivProb;
	array_3d paedSurvCd4Dist;
	array_3d survivalRate;
	array_3d netMigration;
	array_2d asfr;
	array_2d sexRatioAtBirth;
	const int HIVSTEPS_PER_YEAR;
	const int DT;
	array_3d cd4Progression;
	array_3d cd4InitialDist;
	array_3d cd4Mortality;
	array_3d incrrAge;
	double relinfectArt;
	double iota;
	std::vector<double> incrrSex;
	bool incidMod;
	int eppMod;
	int scaleCd4Mortality;
	std::vector<double> projectionSteps;
	double timeEpidemicStart;
	double prevalenceCurrent;
	// Store prevalence at last time step for r-trend model
	double prevalenceLast;
	std::vector<int> ageGroupsStart;
	std::vector<double> rSplineRVec;
	std::vector<double> rTrendBeta;
	double rTrendTStab;
	double rTrendR0;
	// Outputs
	array_1d entrantPrevOut;
	array_1d incidence15to49;
	array_1d prevalence15to49;
	array_1d incrate15To49;
	array_1d rVec;
	array_2d infections;

public:
	State state;
	ArtData art;
	Model(State modelState, ArtData artData, std::vector<double> ageGroupsSp, std::vector<double> vertTLag, std::vector<double> paedSurvLag,
	      bool popAdjust, array_2d entrantPop, array_2d birthLag, array_2d cumSurvey,
	      array_2d cumNetMigr, double netMigrationHivProb, array_3d paedSurvCd4Distrib,
	      array_3d survRate,
	      array_3d netMigr, array_2d asfRate, array_2d sexRatioBirth, int hivStepsPerYear,
	      array_3d cd4Prog, array_3d cd4InitDist, array_3d cd4Mort, array_3d incrrAges, int tArtStart,
	      double artRelinfect, int eppModel, int scaleCd4Mort, std::vector<double> projSteps)
		: state(modelState),
		  art(artData),
		  entrantPrev(boost::extents[PROJECTION_YEARS][SEXES]),
		  previousPregnancyLag(PROJECTION_YEARS, 0.0),
		  entrantPopulation(boost::extents[PROJECTION_YEARS][SEXES]),
		  birthsLag(boost::extents[PROJECTION_YEARS][SEXES]),
		  cumulativeSurvey(boost::extents[PROJECTION_YEARS][SEXES]),
		  cumulativeNetMigr(boost::extents[PROJECTION_YEARS][SEXES]),
		  paedSurvCd4Dist(boost::extents[PROJECTION_YEARS][SEXES][CD4_STAGES]),
		  survivalRate(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  netMigration(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  asfr(boost::extents[PROJECTION_YEARS][FERT_AGES]),
		  sexRatioAtBirth(boost::extents[PROJECTION_YEARS][SEXES]),
		  HIVSTEPS_PER_YEAR(hivStepsPerYear),
		  DT(1.0 / HIVSTEPS_PER_YEAR),
		  cd4Progression(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES - 1]),
		  cd4InitialDist(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		  cd4Mortality(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		  incrrAge(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  ageGroupsStart(AGE_GROUPS, 0),
		  entrantPrevOut(boost::extents[PROJECTION_YEARS]),
		  incidence15to49(boost::extents[PROJECTION_YEARS]),
		  prevalence15to49(boost::extents[(PROJECTION_YEARS - 1) * hivStepsPerYear]),
		  incrate15To49(boost::extents[(PROJECTION_YEARS - 1) * hivStepsPerYear]),
		  rVec(boost::extents[(PROJECTION_YEARS - 1) * hivStepsPerYear]),
		  infections(boost::extents[SEXES][MODEL_AGES]) {

		timeArtStart = tArtStart;
		populationAdjust = popAdjust;
		netMigrHivProb = netMigrationHivProb;
		survivalRate = survRate;
		vertTransLag = vertTLag;
		paedSurveyLag = paedSurvLag;
		entrantPopulation = entrantPop;
		birthsLag = birthLag;
		cumulativeSurvey = cumSurvey;
		cumulativeNetMigr = cumNetMigr;
		paedSurvCd4Dist = paedSurvCd4Distrib;
		netMigration = netMigr;
		asfr = asfRate;
		sexRatioAtBirth = sexRatioBirth;
		cd4Progression = cd4Prog;
		cd4InitialDist = cd4InitDist;
		cd4Mortality = cd4Mort;
		incrrAge = incrrAges;
		relinfectArt = artRelinfect;
		eppMod = eppModel;
		scaleCd4Mortality = scaleCd4Mort;
		projectionSteps = projSteps;
		ageGroupsSpan = ageGroupsSp;
		incidMod = FALSE;
		prevalenceCurrent = 0.0;
		useEntrantPrev = FALSE;

		for (int ageGroup = 1; ageGroup < AGE_GROUPS; ageGroup++) {
			ageGroupsStart[ageGroup] = ageGroupsStart[ageGroup - 1] + ageGroupsSpan[ageGroup - 1];
		}

		// Prepare outputs
		prevalence15to49[0] = 0.0;

	};

	void updateModelState(int t) {
		state.updatePopulation(t);
		state.updateArtPopulation(art.population);
		state.updateHivPopulation();
		state.updateNaturalDeaths();
		state.updateInfections();
	};

	void setEntrantPrev(array_2d entPrev) {
		useEntrantPrev = TRUE;
		entrantPrev = entPrev;
	};

	void setIncrrSex(std::vector<double> incrrSexRatio) {
		incidMod = TRUE;
		incrrSex = incrrSexRatio;
	}

	void initialiseRSpline(std::vector<double> rSplineVec) {
		rSplineRVec = rSplineVec;
	};

	void initialiseRTrend(std::vector<double> beta, double tStab, double r0) {
		rTrendBeta = beta;
		rTrendTStab = tStab;
		rTrendR0 = r0;
	};

	void intialiseNonDirectIncid(double iotaModel, double tsEpidemicStart) {
		iota = iotaModel;
		timeEpidemicStart = tsEpidemicStart;
	};

	void agePopulation(int t);

	void deathsAndMigration(int t);

	void fertility(int t);

	void diseaseProgression(int t);

	double calcRtrendRt(int hivStep, int ts);

	void calcInfectionsEppSpectrum(double iota, int t, int hts, int ts, double infectionsBySexAge[SEXES][MODEL_AGES]);

};

#endif