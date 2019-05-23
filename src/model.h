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
#include "HivData.h"
#include "Metadata.h"

class Model {
private:
	bool useEntrantPrev;
	array_2d entrantPrev;
	std::vector<double> previousPregnancyLag;
	const int DT;
	double iota;
	double timeEpidemicStart;
	std::vector<double> rSplineRVec;
	std::vector<double> rTrendBeta;
	double rTrendTStab;
	double rTrendR0;
	array_2d birthsLag;
	// Outputs
	array_1d entrantPrevOut;
	array_1d incidence15to49;
	array_1d prevalence15to49;
	array_1d incrate15To49;
	array_1d rVec;
	array_2d infections;

public:
	State state;
	const ArtData art;
	const Cd4Data cd4;
	const InfectionData infection;
	const PopulationData pop;
	const HivData hiv;
	const Metadata meta;
	Model(State modelState, const ArtData& artData, const Cd4Data& cd4Data, const InfectionData& infectionData,
	      const PopulationData& populationData, const HivData& hivData, const Metadata& metadata, array_2d birthLag)
		: state(modelState),
		  art(artData),
		  cd4(cd4Data),
		  infection(infectionData),
		  pop(populationData),
		  hiv(hivData),
		  meta(metadata),
		  birthsLag(boost::extents[PROJECTION_YEARS][SEXES]),
		  entrantPrev(boost::extents[PROJECTION_YEARS][SEXES]),
		  previousPregnancyLag(PROJECTION_YEARS, 0.0),
		  DT(1.0 / metadata.hivStepsPerYear),
		  entrantPrevOut(boost::extents[PROJECTION_YEARS]),
		  incidence15to49(boost::extents[PROJECTION_YEARS]),
		  prevalence15to49(boost::extents[(PROJECTION_YEARS - 1) * metadata.hivStepsPerYear]),
		  incrate15To49(boost::extents[(PROJECTION_YEARS - 1) * metadata.hivStepsPerYear]),
		  rVec(boost::extents[(PROJECTION_YEARS - 1) * metadata.hivStepsPerYear]),
		  infections(boost::extents[SEXES][MODEL_AGES]) {

		birthsLag = birthLag;
		useEntrantPrev = false;


		// Prepare outputs
		prevalence15to49[0] = 0.0;

	};

	void updateModelState(int t) {
		state.updatePopulation(t);
		state.updateArtPopulation();
		state.updateHivPopulation();
		state.updateNaturalDeaths();
		state.updateInfections();
	};

	void setEntrantPrev(array_2d entPrev) {
		useEntrantPrev = true;
		entrantPrev = entPrev;
	};

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