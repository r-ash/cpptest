#ifndef _CPPTEST_POPULATIONDATA_H_
#define _CPPTEST_POPULATIONDATA_H_

#include "read_array.h"
#include "consts.h"

struct PopulationData {
	array_2d entrantPopulation;
	array_2d cumulativeSurvey;
	array_2d cumulativeNetMigr;
	array_3d survivalRate;
	array_3d netMigration;
	array_2d asfr;
	array_2d sexRatioAtBirth;

	PopulationData(array_2d entrantPop, array_2d cumSurvey,
	               array_2d cumNetMigr, array_3d survRate, array_3d netMigr,
	               array_2d asfRate, array_2d sexRatioBirth)
		: entrantPopulation(boost::extents[PROJECTION_YEARS][SEXES]),
		  cumulativeSurvey(boost::extents[PROJECTION_YEARS][SEXES]),
		  cumulativeNetMigr(boost::extents[PROJECTION_YEARS][SEXES]),
		  survivalRate(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  netMigration(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  asfr(boost::extents[PROJECTION_YEARS][FERT_AGES]),
		  sexRatioAtBirth(boost::extents[PROJECTION_YEARS][SEXES])  {

		entrantPopulation = entrantPop;
		cumulativeSurvey = cumSurvey;
		cumulativeNetMigr = cumNetMigr;
		survivalRate = survRate;
		netMigration = netMigr;
		asfr = asfRate;
		sexRatioAtBirth = sexRatioBirth;
	}

};

#endif