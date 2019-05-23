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
	std::vector<int> ageGroupsSpan;
	std::vector<int> ageGroupsStart;
	int scaleCd4Mortality;
	const int hivStepsPerYear;
	Metadata(int tArtStart, double artRelinfect, int eppModel, std::vector<double> projSteps,
	         bool popAdjust, std::vector<int> ageGroupsSp, int scaleCd4Mort, int hivStepsPerYr)
		: hivStepsPerYear(hivStepsPerYr),
		  ageGroupsStart(AGE_GROUPS, 0) {
		timeArtStart = tArtStart;
		relinfectArt = artRelinfect;
		eppMod = eppModel;
		projectionSteps = projSteps;
		populationAdjust = popAdjust;
		ageGroupsSpan = ageGroupsSp;
		scaleCd4Mortality = scaleCd4Mort;

		for (int ageGroup = 1; ageGroup < AGE_GROUPS; ageGroup++) {
			ageGroupsStart[ageGroup] = ageGroupsStart[ageGroup - 1] + ageGroupsSpan[ageGroup - 1];
		}
	}

};

#endif