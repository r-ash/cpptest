#ifndef _CPPTEST_INTERFACE_H_
#define _CPPTEST_INTERFACE_H_

#include <Rcpp.h>
#include "read_array.h"
#include "model.h"

std::vector<double> runModel(std::vector<double> basePop, std::vector<double> ageGroupsSpan, int timeArtStart,
                             SEXP entrantPrev, std::vector<double> vertTransLag, std::vector<double> paedSurveyLag,
                             bool populationAdjust, std::vector<double> entrantPop, std::vector<double> birthLag,
                             std::vector<double> cumSurv, std::vector<double> cumNetMigr, double netMigrHivProb,
                             std::vector<double> paedSurvCd4Distrib, std::vector<double> paedSurvArtCd4Distrib,
                             std::vector<double> survRate, std::vector<double> netMigr, std::vector<double> asfRate,
                             std::vector<double> sexRatioBirth, int timeSteps);

#endif