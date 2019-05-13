#ifndef _CPPTEST_INTERFACE_H_
#define _CPPTEST_INTERFACE_H_

#include <Rcpp.h>
#include "read_array.h"
#include "model.h"
#include "state.h"
#include "ArtData.h"
#include "Cd4Data.h"
#include "InfectionData.h"
#include "PopulationData.h"
#include "HivData.h"

std::vector<double> runModel(std::vector<double> basePop, std::vector<double> ageGroupsSpan, int timeArtStart,
                             SEXP entrantPrev, std::vector<double> vertTransLag, std::vector<double> paedSurveyLag,
                             bool populationAdjust, std::vector<double> entrantPop, std::vector<double> birthLag,
                             std::vector<double> cumSurv, std::vector<double> cumNetMigr, double netMigrHivProb,
                             std::vector<double> paedSurvCd4Distrib, SEXP entrantArtCoverage, std::vector<double> paedSurvArtCd4Distrib,
                             std::vector<double> survRate, std::vector<double> netMigr, std::vector<double> asfRate,
                             std::vector<double> sexRatioBirth, int hivStepsPerYear, std::vector<double> cd4Prog,
                             std::vector<double> cd4InitDist, std::vector<double> cd4Mort, std::vector<double> incrrAges,
                             double relinfectArt, SEXP iota, std::vector<double> incrrSex, SEXP incidMod, int eppMod,
                             int scaleCd4Mort, std::vector<double> projSteps, SEXP tsEpidemicStart, std::vector<int> artCd4EligId,
                             std::vector<double> specPopPercentElig, std::vector<double> pregnantWomenArtElig,
                             double who34PercentElig, SEXP rSplineRVec, SEXP rTrendBeta, SEXP rTrendTStab,
                             SEXP rTrendR0, int timeSteps);

#endif