#ifndef _CPPTEST_STRATEGYMODEL_H_
#define _CPPTEST_STRATEGYMODEL_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "consts.h"
#include "AgeMethod.h"
#include "state.h"
#include "ArtData.h"
#include "Cd4Data.h"
#include "InfectionData.h"
#include "PopulationData.h"
#include "HivData.h"
#include "Metadata.h"

class StrategyModel {
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
	std::unique_ptr<AgeMethod> thisAgeMethod;
public:
	State state;
	const ArtData art;
	const Cd4Data cd4;
	const InfectionData infection;
	const PopulationData pop;
	const HivData hiv;
	const Metadata meta;
	StrategyModel(State modelState, const ArtData& artData, const Cd4Data& cd4Data, const InfectionData& infectionData,
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

	void setAgeMethod(AgeMethod *ageMethod) {
		thisAgeMethod.reset(ageMethod);
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

	void runModel(int timeSteps) {
		for (int t = 1; t <= timeSteps; t++) {
			// thisAgeMethod.agePopulation(t, state, hiv, pop, meta, cd4, art, useEntrantPrev,
			//                             entrantPrev, previousPregnancyLag, birthsLag,
			//                             entrantPrevOut);
			thisAgeMethod->agePopulation();
			updateModelState(t);
		}
	};
};

#endif