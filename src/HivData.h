#ifndef _CPPTEST_HIVDATA_H_
#define _CPPTEST_HIVDATA_H_

#include "read_array.h"
#include "consts.h"

class HivData {
public:
	double netMigrProb;
	std::vector<double> vertTransLag;
	std::vector<double> paedSurveyLag;
	HivData(double netMigrationHivProb, std::vector<double> vertTLag,
	        std::vector<double> paedSurvLag) {
		netMigrProb = netMigrationHivProb;
		vertTransLag = vertTLag;
		paedSurveyLag = paedSurvLag;
	}

};

#endif