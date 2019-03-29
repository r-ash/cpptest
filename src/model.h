#ifndef _CPPTEST_MODEL_H_
#define _CPPTEST_MODEL_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "read_array.h"
#include "consts.h"

typedef boost::multi_array<double, 1> array_1d;
typedef boost::multi_array<double, 2> array_2d;
typedef boost::multi_array<double, 3> array_3d;
typedef boost::multi_array<double, 4> array_4d;

class Model {
private:
	array_4d population;

public:
	Model(array_2d basePopulation) : population(boost::extents[PROJECTION_YEARS][DISEASE_STATUS][SEXES][MODEL_AGES]) {
		for (int sex = 0; sex < SEXES; sex++) {
			for (int age = 0; age < MODEL_AGES; age++) {
				population[0][HIVN][sex][age] = basePopulation[sex][age];
				population[0][HIVP][sex][age] = 0.0;
			}
		}
	};

	array_4d getPopulation() {
		return population;
	};

	void agePopulation(int t) {
		for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
			for (int sex = 0; sex < SEXES; sex++) {
				for (int age = 1; age < MODEL_AGES; age++) {
					Rcpp::Rcout << "Looping over sex " << sex << " and age " << age << " \n";
					population[t][diseaseStatus][sex][age] = population[t - 1][diseaseStatus][sex][age - 1];
				}
				// People do not age out of final open age group (80+)
				population[t][diseaseStatus][sex][MODEL_AGES - 1] += population[t - 1][diseaseStatus][sex][MODEL_AGES - 1];
			}
		}

		// double hiv_ag_prob[SEXES][AGE_GROUPS];
		// for (int sex = 0; sex < SEXES; sex++) {
		// 	int a = 0;
		// 	for (int ageGroup = 0; ageGroup < (AGE_GROUPS - 1); ageGroup++) {
		// 		hiv_ag_prob[sex][ageGroup] = 0;
		// 		for (int i = 0; i < AGE_GROUPS_SPAN[ageGroup]; i++) {
		// 			hiv_ag_prob[sex][ageGroup] += pop[t - 1][HIVP][sex][a];
		// 			a++;
		// 		}
		// 		hiv_ag_prob[sex][ageGroup] = (hiv_ag_prob[sex][ageGroup] > 0) ? pop[t - 1][HIVP][sex][a - 1] / hiv_ag_prob[sex][ageGroup] : 0;
		// 	}
		// 	// People do not age out of final open age group (80+)
		// 	hiv_ag_prob[sex][AGE_GROUPS - 1] = 0.0;
		// }

		// for (int sex = 0; sex < SEXES; sex++) {
		// 	for (int ageGroup = 1; ageGroup < AGE_GROUPS; ageGroup++) {
		// 		for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
		// 			hivpop[t][sex][ageGroup][cd4Stage] = (1 - hiv_ag_prob[sex][ageGroup]) * hivpop[t - 1][sex][ageGroup][cd4Stage] + hiv_ag_prob[sex][ageGroup - 1] * hivpop[t - 1][sex][ageGroup - 1][cd4Stage];
		// 			if (t > t_ART_start)
		// 				for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
		// 					artpop[t][sex][ageGroup][cd4Stage][treatmentStage] = (1 - hiv_ag_prob[sex][ageGroup]) * artpop[t - 1][sex][ageGroup][cd4Stage][treatmentStage] + hiv_ag_prob[sex][ageGroup - 1] * artpop[t - 1][sex][ageGroup - 1][cd4Stage][treatmentStage];
		// 				}
		// 		}
		// 	}
		// }

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
		// 		hivpop[t][sex][0][cd4Stage] = (1 - hiv_ag_prob[sex][0]) * hivpop[t - 1][sex][0][cd4Stage] + paedsurv_g * paedsurv_cd4dist[t][sex][cd4Stage] * (1.0 - entrantartcov[t][sex]);
		// 		if (t > t_ART_start) {
		// 			for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
		// 				artpop[t][sex][0][cd4Stage][treatmentStage] = (1 - hiv_ag_prob[sex][0]) * artpop[t - 1][sex][0][cd4Stage][treatmentStage];
		// 				artpop[t][sex][0][cd4Stage][treatmentStage] += paedsurv_g * paedsurv_artcd4dist[t][sex][cd4Stage][treatmentStage] * entrantartcov[t][sex];
		// 			}
		// 		}
		// 	}
		// }
	};
};

std::vector<double> runModel(std::vector<double> basePop, int timeSteps);

#endif