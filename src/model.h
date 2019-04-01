#ifndef _CPPTEST_MODEL_H_
#define _CPPTEST_MODEL_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "read_array.h"
#include "consts.h"


class Model {
private:
	array_4d population;
	array_1d ageGroupsSpan;
	int timeArtStart;
	array_5d artPopulation;
	array_4d hivPopulation;
	bool useEntrantPrev = FALSE;
	array_2d entrantPrev;

public:
	Model(array_2d basePopulation, array_1d ageGroupsSp, int tArtStart)
		: population(boost::extents[PROJECTION_YEARS][DISEASE_STATUS][SEXES][MODEL_AGES]),
		  ageGroupsSpan(boost::extents[AGE_GROUPS]),
		  artPopulation(boost::extents[PROJECTION_YEARS][SEXES][AGE_GROUPS][DISEASE_STATUS][CD4_STAGES]),
		  hivPopulation(boost::extents[PROJECTION_YEARS][SEXES][AGE_GROUPS][DISEASE_STATUS]) {

		timeArtStart = tArtStart;

		for (int sex = 0; sex < SEXES; sex++) {
			for (int age = 0; age < MODEL_AGES; age++) {
				population[0][HIVN][sex][age] = basePopulation[sex][age];
				population[0][HIVP][sex][age] = 0.0;
			}
		}

		for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
			ageGroupsSpan[ageGroup] = ageGroupsSp[ageGroup];
		}

		for (int sex = 0; sex < SEXES; sex++) {
			for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
				for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
					hivPopulation[0][sex][ageGroup][diseaseStatus] = 0.0;
				}
			}
		}

		if (timeArtStart < PROJECTION_YEARS) {
			for (int sex = 0; sex < SEXES; sex++) {
				for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
					for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
						for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
							// Initialise to 0 in year of ART start
							artPopulation[timeArtStart][sex][ageGroup][diseaseStatus][cd4Stage] = 0.0;
						}
					}
				}
			}
		}
	};

	void setEntrantPrev(array_2d entPrev) {
		useEntrantPrev = TRUE;
		entrantPrev = entPrev;
	}

	array_4d getPopulation() {
		return population;
	};

	void agePopulation(int t) {
		for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
			for (int sex = 0; sex < SEXES; sex++) {
				for (int age = 1; age < MODEL_AGES; age++) {
					//Rcpp::Rcout << "Looping over sex " << sex << " and age " << age << " \n";
					population[t][diseaseStatus][sex][age] = population[t - 1][diseaseStatus][sex][age - 1];
				}
				// People do not age out of final open age group (80+)
				population[t][diseaseStatus][sex][MODEL_AGES - 1] += population[t - 1][diseaseStatus][sex][MODEL_AGES - 1];
			}
		}

		double hivAgProb[SEXES][AGE_GROUPS];
		for (int sex = 0; sex < SEXES; sex++) {
			int a = 0;
			for (int ageGroup = 0; ageGroup < (AGE_GROUPS - 1); ageGroup++) {
				hivAgProb[sex][ageGroup] = 0;
				for (int i = 0; i < ageGroupsSpan[ageGroup]; i++) {
					hivAgProb[sex][ageGroup] += population[t - 1][HIVP][sex][a];
					a++;
				}
				hivAgProb[sex][ageGroup] = (hivAgProb[sex][ageGroup] > 0) ? population[t - 1][HIVP][sex][a - 1] / hivAgProb[sex][ageGroup] : 0;
			}
			// People do not age out of final open age group (80+)
			hivAgProb[sex][AGE_GROUPS - 1] = 0.0;
		}

		for (int sex = 0; sex < SEXES; sex++) {
			for (int ageGroup = 1; ageGroup < AGE_GROUPS; ageGroup++) {
				for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
					hivPopulation[t][sex][ageGroup][cd4Stage] = (1 - hivAgProb[sex][ageGroup]) * hivPopulation[t - 1][sex][ageGroup][cd4Stage] + hivAgProb[sex][ageGroup - 1] * hivPopulation[t - 1][sex][ageGroup - 1][cd4Stage];
					if (t > timeArtStart) {
						for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
							artPopulation[t][sex][ageGroup][cd4Stage][treatmentStage] = (1 - hivAgProb[sex][ageGroup]) * artPopulation[t - 1][sex][ageGroup][cd4Stage][treatmentStage] + hivAgProb[sex][ageGroup - 1] * artPopulation[t - 1][sex][ageGroup - 1][cd4Stage][treatmentStage];
						}
					}
				}
			}
		}

		// // add lagged births to youngest age group
		// for (int sex = 0; sex < SEXES; sex++) {

		// 	double paedsurv_g;
		// 	double entrant_prev;

		// 	if (use_entrantprev) {
		// 		entrant_prev = entrantprev[t][sex];
		// 	} else {
		// 		entrant_prev = pregprevlag[t - 1] * verttrans_lag[t - 1] * paedsurv_lag[t - 1];
		// 	}

		// 	if (bin_popadjust) {
		// 		pop[t][HIVN][sex][0] =  entrantpop[t - 1][sex] * (1.0 - entrant_prev);
		// 		paedsurv_g = entrantpop[t - 1][sex] * entrant_prev;
		// 	} else {
		// 		pop[t][HIVN][sex][0] = birthslag[t - 1][sex] * cumsurv[t - 1][sex] * (1.0 - entrant_prev / paedsurv_lag[t - 1]) + cumnetmigr[t - 1][sex] * (1.0 - pregprevlag[t - 1] * netmig_hivprob);
		// 		paedsurv_g = birthslag[t - 1][sex] * cumsurv[t - 1][sex] * entrant_prev + cumnetmigr[t - 1][sex] * entrant_prev;
		// 	}

		// 	pop[t][HIVP][sex][0] = paedsurv_g;

		// 	entrantprev_out[t] = (pop[t][HIVP][MALE][0] + pop[t][HIVP][FEMALE][0]) / (pop[t][HIVN][MALE][0] + pop[t][HIVN][FEMALE][0] + pop[t][HIVP][MALE][0] + pop[t][HIVP][FEMALE][0]);

		// 	for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
		// 		hivpop[t][sex][0][cd4Stage] = (1 - hivAgProb[sex][0]) * hivpop[t - 1][sex][0][cd4Stage] + paedsurv_g * paedsurv_cd4dist[t][sex][cd4Stage] * (1.0 - entrantartcov[t][sex]);
		// 		if (t > timeArtStart) {
		// 			for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
		// 				artpop[t][sex][0][cd4Stage][treatmentStage] = (1 - hivAgProb[sex][0]) * artpop[t - 1][sex][0][cd4Stage][treatmentStage];
		// 				artpop[t][sex][0][cd4Stage][treatmentStage] += paedsurv_g * paedsurv_artcd4dist[t][sex][cd4Stage][treatmentStage] * entrantartcov[t][sex];
		// 			}
		// 		}
		// 	}
		// }
	};
};

std::vector<double> runModel(std::vector<double> basePop, std::vector<double> ageGroupsSpan, int timeArtStart,
                             std::vector<double> entrantPrev, int timeSteps);

std::vector<double> runModel(std::vector<double> basePop, std::vector<double> ageGroupsSpan, int timeArtStart,
                             int timeSteps);

#endif