#ifndef _CPPTEST_AGEMETHOD_H_
#define _CPPTEST_AGEMETHOD_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "consts.h"
#include "state.h"
#include "ArtData.h"
#include "Cd4Data.h"
#include "read_array.h"
#include "InfectionData.h"
#include "PopulationData.h"
#include "HivData.h"
#include "Metadata.h"

class AgeMethod {
public:
	virtual void agePopulation(int t, State state, HivData hiv, PopulationData pop, Metadata meta,
	                           Cd4Data cd4, ArtData art, bool useEntrantPrev, array_2d entrantPrev,
	                           std::vector<double> previousPregnancyLag, array_2d birthsLag,
	                           array_1d entrantPrevOut);
};

class DefaultAgeMethod: public AgeMethod {
public:
	void agePopulation(int t, State state, HivData hiv, PopulationData pop, Metadata meta,
	                   Cd4Data cd4, ArtData art, bool useEntrantPrev, array_2d entrantPrev,
	                   std::vector<double> previousPregnancyLag, array_2d birthsLag,
	                   array_1d entrantPrevOut) {
		for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
			for (int sex = 0; sex < SEXES; sex++) {
				for (int age = 1; age < MODEL_AGES; age++) {
					// Rcpp::Rcout << "Looping over sex " << sex << " and age " << age << " \n";
					state.population[diseaseStatus][sex][age] = state.previousPopulation[diseaseStatus][sex][age - 1];
				}
				// People do not age out of final open age group (80+)
				state.population[diseaseStatus][sex][MODEL_AGES - 1] += state.previousPopulation[diseaseStatus][sex][MODEL_AGES - 1];
			}
		}

		double hivAgProb[SEXES][AGE_GROUPS];
		for (int sex = 0; sex < SEXES; sex++) {
			int a = 0;
			for (int ageGroup = 0; ageGroup < (AGE_GROUPS - 1); ageGroup++) {
				hivAgProb[sex][ageGroup] = 0;
				for (int i = 0; i < meta.ageGroupsSpan[ageGroup]; i++) {
					hivAgProb[sex][ageGroup] += state.previousPopulation[HIVP][sex][a];
					a++;
				}
				hivAgProb[sex][ageGroup] = (hivAgProb[sex][ageGroup] > 0) ? state.previousPopulation[HIVP][sex][a - 1] / hivAgProb[sex][ageGroup] : 0;
			}
			// People do not age out of final open age group (80+)
			hivAgProb[sex][AGE_GROUPS - 1] = 0.0;
		}

		for (int sex = 0; sex < SEXES; sex++) {
			for (int ageGroup = 1; ageGroup < AGE_GROUPS; ageGroup++) {
				for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
					state.hivPop[sex][ageGroup][cd4Stage] = (1 - hivAgProb[sex][ageGroup]) * state.previousHivPop[sex][ageGroup][cd4Stage] + hivAgProb[sex][ageGroup - 1] * state.previousHivPop[sex][ageGroup - 1][cd4Stage];
					if (t > meta.timeArtStart) {
						for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
							state.artPopulation[sex][ageGroup][cd4Stage][treatmentStage] = (1 - hivAgProb[sex][ageGroup]) * state.previousArtPopulation[sex][ageGroup][cd4Stage][treatmentStage] + hivAgProb[sex][ageGroup - 1] * state.previousArtPopulation[sex][ageGroup - 1][cd4Stage][treatmentStage];
						}
					}
				}
			}
		}

		// add lagged births to youngest age group
		for (int sex = 0; sex < SEXES; sex++) {

			double entrantPrevalence;

			if (useEntrantPrev) {
				entrantPrevalence = entrantPrev[t][sex];
			} else {
				entrantPrevalence = previousPregnancyLag[t - 1] * hiv.vertTransLag[t - 1] * hiv.paedSurveyLag[t - 1];
			}

			if (meta.populationAdjust) {
				state.population[HIVN][sex][0] = pop.entrantPopulation[t - 1][sex] * (1.0 - entrantPrevalence);
				state.population[HIVP][sex][0] = pop.entrantPopulation[t - 1][sex] * entrantPrevalence;
			} else {
				state.population[HIVN][sex][0] = birthsLag[t - 1][sex] * pop.cumulativeSurvey[t - 1][sex] * (1.0 - entrantPrevalence / hiv.paedSurveyLag[t - 1]) + pop.cumulativeNetMigr[t - 1][sex] * (1.0 - previousPregnancyLag[t - 1] * hiv.netMigrProb);
				state.population[HIVP][sex][0] = birthsLag[t - 1][sex] * pop.cumulativeSurvey[t - 1][sex] * entrantPrevalence + pop.cumulativeNetMigr[t - 1][sex] * entrantPrevalence;
			}

			double paedSurvPos = state.population[HIVP][sex][0];

			entrantPrevOut[t] = (state.population[HIVP][MALE][0] + state.population[HIVP][FEMALE][0]) /
			                    (state.population[HIVN][MALE][0] + state.population[HIVN][FEMALE][0] + state.population[HIVP][MALE][0] + state.population[HIVP][FEMALE][0]);

			for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
				state.hivPop[sex][0][cd4Stage] = (1 - hivAgProb[sex][0]) * state.previousHivPop[sex][0][cd4Stage] + paedSurvPos * cd4.paedSurvDist[t][sex][cd4Stage] * (1.0 - art.entrantCoverage[t][sex]);
				if (t > meta.timeArtStart) {
					for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
						state.artPopulation[sex][0][cd4Stage][treatmentStage] = (1 - hivAgProb[sex][0]) * state.previousArtPopulation[sex][0][cd4Stage][treatmentStage];
						state.artPopulation[sex][0][cd4Stage][treatmentStage] += paedSurvPos * art.paedSurvArtCd4Dist[t][sex][cd4Stage][treatmentStage] * art.entrantCoverage[t][sex];
					}
				}
			}
		}
	}
};


#endif