#ifndef _CPPTEST_INFECTIONDATA_H_
#define _CPPTEST_INFECTIONDATA_H_

#include "read_array.h"
#include "consts.h"

class InfectionData {
public:
	array_3d incrrAge;
	std::vector<double> incrrSex;
	bool incidMod;

	InfectionData(array_3d incrrAges)
		: incrrAge(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]) {
		incrrAge = incrrAges;
		incidMod = false;
	}

	void setIncrrSex(std::vector<double> incrrSexRatio) {
		incidMod = true;
		incrrSex = incrrSexRatio;
	}

};

#endif