#ifndef _CPPTEST_MODEL_H_
#define _CPPTEST_MODEL_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "read_array.h"
#include "consts.h"
#include "state.h"
#include "ArtData.h"
#include "Cd4Data.h"
#include "InfectionData.h"
#include "PopulationData.h"

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
	double netMigrHivProb;
	const int HIVSTEPS_PER_YEAR;
	const int DT;
	double relinfectArt;
	double iota;
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
	Cd4Data cd4;
	InfectionData infection;
	PopulationData pop;
	Model(State modelState, ArtData artData, Cd4Data cd4Data, InfectionData infectionData, PopulationData populationData,
	      std::vector<double> ageGroupsSp, std::vector<double> vertTLag, std::vector<double> paedSurvLag,
	      bool popAdjust, double netMigrationHivProb,
	      int hivStepsPerYear,
	      int tArtStart,
	      double artRelinfect, int eppModel, int scaleCd4Mort, std::vector<double> projSteps)
		: state(modelState),
		  art(artData),
		  cd4(cd4Data),
		  infection(infectionData),
		  pop(populationData),
		  entrantPrev(boost::extents[PROJECTION_YEARS][SEXES]),
		  previousPregnancyLag(PROJECTION_YEARS, 0.0),
		  HIVSTEPS_PER_YEAR(hivStepsPerYear),
		  DT(1.0 / HIVSTEPS_PER_YEAR),
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
		vertTransLag = vertTLag;
		paedSurveyLag = paedSurvLag;
		relinfectArt = artRelinfect;
		eppMod = eppModel;
		scaleCd4Mortality = scaleCd4Mort;
		projectionSteps = projSteps;
		ageGroupsSpan = ageGroupsSp;
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
		infection.incidMod = TRUE;
		infection.incrrSex = incrrSexRatio;
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