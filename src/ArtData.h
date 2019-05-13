#ifndef _CPPTEST_ARTDATA_H_
#define _CPPTEST_ARTDATA_H_

#include "read_array.h"
#include "consts.h"

class ArtData {
public:
	array_2d entrantCoverage;
	array_4d paedSurvArtCd4Dist;
	std::vector<int> artCd4EligId;
	std::vector<double> specPopPercentElig;
	std::vector<double> pregnantArtElig;
	double who34PercentElig;
	int everArtEligId;

	ArtData(array_2d entrantArtCov, array_4d paedSurvArtCd4Distrib, std::vector<int> artCd4EligibleId,
	        std::vector<double> specPopulationPercentElig, std::vector<double> pregnantWomenArtElig,
	        double who34PercElig)
		: entrantCoverage(boost::extents[PROJECTION_YEARS][SEXES]),
		  paedSurvArtCd4Dist(boost::extents[PROJECTION_YEARS][SEXES][CD4_STAGES][TREATMENT_STAGES]) {

		entrantCoverage = entrantArtCov;
		paedSurvArtCd4Dist = paedSurvArtCd4Distrib;
		artCd4EligId = artCd4EligibleId;
		specPopPercentElig = specPopulationPercentElig;
		pregnantArtElig = pregnantWomenArtElig;
		who34PercentElig = who34PercElig;
		everArtEligId = CD4_STAGES;
	}


	void calculateEverArtEligId(int t) {
		int cd4EligId = artCd4EligId[t] - 1; // -1 for 0-based indexing vs 1-based in R
		int anyEligId = ((specPopPercentElig[t] > 0) | (pregnantArtElig[t] > 0)) ? 0 : (who34PercentElig > 0) ? hIDX_CD4_350 : cd4EligId;
		everArtEligId = anyEligId < everArtEligId ? anyEligId : everArtEligId;
	}
};



#endif