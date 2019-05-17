#ifndef _CPPTEST_METADATA_H_
#define _CPPTEST_METADATA_H_

#include "read_array.h"
#include "consts.h"

struct Metadata {
	int timeArtStart;
	double relinfectArt;
	int eppMod;
	std::vector<double> projectionSteps;
	bool populationAdjust;
	std::vector<double> ageGroupsSpan;
	int scaleCd4Mortality;
	const int hivStepsPerYear;
	Metadata(int tArtStart, double artRelinfect, int eppModel, std::vector<double> projSteps,
	         bool popAdjust, std::vector<double> ageGroupsSp, int scaleCd4Mort, int hivStepsPerYr)
		: hivStepsPerYear(hivStepsPerYr) {
		timeArtStart = tArtStart;
		relinfectArt = artRelinfect;
		eppMod = eppModel;
		projectionSteps = projSteps;
		populationAdjust = popAdjust;
		ageGroupsSpan = ageGroupsSp;
		scaleCd4Mortality = scaleCd4Mort;
	}

};

#endif