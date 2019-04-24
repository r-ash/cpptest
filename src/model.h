#ifndef _CPPTEST_MODEL_H_
#define _CPPTEST_MODEL_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "read_array.h"
#include "consts.h"


class Model {
private:
	array_4d population;
	std::vector<double> ageGroupsSpan;
	int timeArtStart;
	array_5d artPopulation;
	array_4d hivPop;
	bool useEntrantPrev = FALSE;
	array_2d entrantPrev;
	array_1d entrantPrevOut;
	array_1d previousPregnancyLag;
	std::vector<double> vertTransLag;
	std::vector<double> paedSurveyLag;
	bool populationAdjust;
	array_2d entrantPopulation;
	array_2d birthsLag;
	array_2d cumulativeSurvey;
	array_2d cumulativeNetMigr;
	double netMigrHivProb;
	array_3d paedSurvCd4Dist;
	array_2d entrantArtCoverage;
	array_4d paedSurvArtCd4Dist;
	array_3d survivalRate;
	array_3d naturalDeaths;
	array_3d netMigration;
	array_2d asfr;
	array_2d sexRatioAtBirth;

public:
	Model(array_2d basePopulation, std::vector<double> ageGroupsSp, std::vector<double> vertTLag, std::vector<double> paedSurvLag,
	      bool popAdjust, array_2d entrantPop, array_2d birthLag, array_2d cumSurvey,
	      array_2d cumNetMigr, double netMigrationHivProb, array_3d paedSurvCd4Distrib,
	      array_2d entrantArtCov, array_4d paedSurvArtCd4Distrib, array_3d survRate,
	      array_3d netMigr, array_2d asfRate, array_2d sexRatioBirth, int tArtStart)
		: population(boost::extents[PROJECTION_YEARS][DISEASE_STATUS][SEXES][MODEL_AGES]),
		  artPopulation(boost::extents[PROJECTION_YEARS][SEXES][AGE_GROUPS][DISEASE_STATUS][CD4_STAGES]),
		  hivPop(boost::extents[PROJECTION_YEARS][SEXES][AGE_GROUPS][CD4_STAGES]),
		  entrantPrev(boost::extents[PROJECTION_YEARS][SEXES]),
		  entrantPrevOut(boost::extents[PROJECTION_YEARS]),
		  previousPregnancyLag(boost::extents[PROJECTION_YEARS]),
		  entrantPopulation(boost::extents[PROJECTION_YEARS][SEXES]),
		  birthsLag(boost::extents[PROJECTION_YEARS][SEXES]),
		  cumulativeSurvey(boost::extents[PROJECTION_YEARS][SEXES]),
		  cumulativeNetMigr(boost::extents[PROJECTION_YEARS][SEXES]),
		  paedSurvCd4Dist(boost::extents[PROJECTION_YEARS][SEXES][CD4_STAGES]),
		  entrantArtCoverage(boost::extents[PROJECTION_YEARS][SEXES]),
		  paedSurvArtCd4Dist(boost::extents[PROJECTION_YEARS][SEXES][CD4_STAGES][TREATMENT_STAGES]),
		  survivalRate(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  naturalDeaths(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  netMigration(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  asfr(boost::extents[PROJECTION_YEARS][FERT_AGES]),
		  sexRatioAtBirth(boost::extents[qPROJECTION_YEARS][SEXES]) {

		timeArtStart = tArtStart;
		populationAdjust = popAdjust;
		netMigrHivProb = netMigrationHivProb;
		entrantArtCoverage = entrantArtCov;
		paedSurvArtCd4Dist = paedSurvArtCd4Distrib;
		survivalRate = survRate;
		vertTransLag = vertTLag;
		paedSurveyLag = paedSurvLag;
		entrantPopulation = entrantPop;
		birthsLag = birthLag;
		cumulativeSurvey = cumSurvey;
		cumulativeNetMigr = cumNetMigr;
		paedSurvCd4Dist = paedSurvCd4Distrib;
		netMigration = netMigr;
		asfr = asfRate;
		sexRatioAtBirth = sexRatioBirth;

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
				for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
					hivPop[0][sex][ageGroup][cd4Stage] = 0.0;
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

		for (int year = 0; year < PROJECTION_YEARS; year++) {
			previousPregnancyLag[year] = 0.0;
			for (int sex = 0; sex < SEXES; sex++) {
				for (int modelAge = 0; modelAge < MODEL_AGES; modelAge++) {
					naturalDeaths[year][sex][modelAge] = 0.0;
				}
			}
		}
	};

	void setEntrantPrev(array_2d entPrev) {
		useEntrantPrev = TRUE;
		entrantPrev = entPrev;
	};

	array_4d getPopulation() {
		return population;
	};

	void agePopulation(int t) {
		for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
			for (int sex = 0; sex < SEXES; sex++) {
				for (int age = 1; age < MODEL_AGES; age++) {
					// Rcpp::Rcout << "Looping over sex " << sex << " and age " << age << " \n";
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
					hivPop[t][sex][ageGroup][cd4Stage] = (1 - hivAgProb[sex][ageGroup]) * hivPop[t - 1][sex][ageGroup][cd4Stage] + hivAgProb[sex][ageGroup - 1] * hivPop[t - 1][sex][ageGroup - 1][cd4Stage];
					if (t > timeArtStart) {
						for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
							artPopulation[t][sex][ageGroup][cd4Stage][treatmentStage] = (1 - hivAgProb[sex][ageGroup]) * artPopulation[t - 1][sex][ageGroup][cd4Stage][treatmentStage] + hivAgProb[sex][ageGroup - 1] * artPopulation[t - 1][sex][ageGroup - 1][cd4Stage][treatmentStage];
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
				entrantPrevalence = previousPregnancyLag[t - 1] * vertTransLag[t - 1] * paedSurveyLag[t - 1];
			}

			if (populationAdjust) {
				population[t][HIVN][sex][0] = entrantPopulation[t - 1][sex] * (1.0 - entrantPrevalence);
				population[t][HIVP][sex][0] = entrantPopulation[t - 1][sex] * entrantPrevalence;
			} else {
				population[t][HIVN][sex][0] = birthsLag[t - 1][sex] * cumulativeSurvey[t - 1][sex] * (1.0 - entrantPrevalence / paedSurveyLag[t - 1]) + cumulativeNetMigr[t - 1][sex] * (1.0 - previousPregnancyLag[t - 1] * netMigrHivProb);
				population[t][HIVP][sex][0] = birthsLag[t - 1][sex] * cumulativeSurvey[t - 1][sex] * entrantPrevalence + cumulativeNetMigr[t - 1][sex] * entrantPrevalence;
			}

			double paedSurvPos = population[t][HIVP][sex][0];

			entrantPrevOut[t] = (population[t][HIVP][MALE][0] + population[t][HIVP][FEMALE][0]) /
			                    (population[t][HIVN][MALE][0] + population[t][HIVN][FEMALE][0] + population[t][HIVP][MALE][0] + population[t][HIVP][FEMALE][0]);

			for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
				hivPop[t][sex][0][cd4Stage] = (1 - hivAgProb[sex][0]) * hivPop[t - 1][sex][0][cd4Stage] + paedSurvPos * paedSurvCd4Dist[t][sex][cd4Stage] * (1.0 - entrantArtCoverage[t][sex]);
				if (t > timeArtStart) {
					for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
						artPopulation[t][sex][0][cd4Stage][treatmentStage] = (1 - hivAgProb[sex][0]) * artPopulation[t - 1][sex][0][cd4Stage][treatmentStage];
						artPopulation[t][sex][0][cd4Stage][treatmentStage] += paedSurvPos * paedSurvArtCd4Dist[t][sex][cd4Stage][treatmentStage] * entrantArtCoverage[t][sex];
					}
				}
			}
		}
	};

	void deathsAndMigration(int t) {
		// Non-HIV mortality
		for (int sex = 0; sex < SEXES; sex++) {
			int a = 0;
			for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
				double deathsMigAgeGroup = 0;
				double hivPopAgeGroup = 0;
				for (int i = 0; i < ageGroupsSpan[ageGroup]; i++) {

					hivPopAgeGroup += population[t][HIVP][sex][a];

					// non-HIV mortality
					double deathRate = 1.0 - survivalRate[t][sex][a];
					double hivNegDeaths = population[t][HIVN][sex][a] * deathRate;
					population[t][HIVN][sex][a] -= hivNegDeaths; // survival HIV- population
					double hivPosDeaths = population[t][HIVP][sex][a] * deathRate;
					deathsMigAgeGroup -= hivPosDeaths;
					population[t][HIVP][sex][a] -= hivPosDeaths;   // survival HIV+ population
					naturalDeaths[t][sex][a] = hivNegDeaths + hivPosDeaths;

					// net migration
					// TODO: rename migrate_a and hmig_a - what do these represent?
					double migrate_a = netMigration[t][sex][a] * (1 + survivalRate[t][sex][a]) / 2.0 / (population[t][HIVN][sex][a] + population[t][HIVP][sex][a]);
					population[t][HIVN][sex][a] *= 1 + migrate_a;
					double hmig_a = migrate_a * population[t][HIVP][sex][a];
					deathsMigAgeGroup += hmig_a;
					population[t][HIVP][sex][a] += hmig_a;

					a++;
				}

				// migration and deaths for hivpop
				double hivPopDeathsMigrRate = hivPopAgeGroup > 0 ? deathsMigAgeGroup / hivPopAgeGroup : 0.0;
				for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
					hivPop[t][sex][ageGroup][cd4Stage] *= 1 + hivPopDeathsMigrRate;
					if (t > timeArtStart) {
						for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
							artPopulation[t][sex][ageGroup][cd4Stage][treatmentStage] *= 1 + hivPopDeathsMigrRate;
						}
					}
				}
			}
		}
	};

	void fertility(int t) {
		Rcpp::Rcout << "Calculating fertility \n";
		double births = 0.0;
		double birthsByAgeGroup[FERT_AGE_GROUPS];
		memset(birthsByAgeGroup, 0, FERT_AGE_GROUPS * sizeof(double));
		for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
			int a = pIDX_FERT;
			for (int ha = hIDX_FERT; ha < hIDX_FERT + FERT_AGE_GROUPS; ha++) {
				for (int i = 0; i < ageGroupsSpan[ha]; i++) {
					birthsByAgeGroup[ha - hIDX_FERT] += (population[t - 1][diseaseStatus][FEMALE][a] + population[t][diseaseStatus][FEMALE][a]) / 2 * asfr[t][a];
					a++;
				}
			}
		}
		for (int ha = hIDX_FERT; ha < FERT_AGE_GROUPS; ha++) {
			births += birthsByAgeGroup[ha - hIDX_FERT];
		}

		if (t + AGE_START < PROJECTION_YEARS) {
			for (int sex = 0; sex < SEXES; sex++) {
				birthsLag[t + AGE_START - 1][sex] = sexRatioAtBirth[t][sex] * births;
			}
		}
	};
};

#endif