#ifndef _CPPTEST_MODEL_H_
#define _CPPTEST_MODEL_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "read_array.h"
#include "consts.h"
#include "state.h"

class Model {
private:
	State state;
	array_3d population;
	array_3d previousPopulation;
	std::vector<double> ageGroupsSpan;
	int timeArtStart;
	array_4d artPopulation;
	array_4d previousArtPopulation;
	array_3d hivPop;
	array_3d previousHivPop;
	bool useEntrantPrev;
	array_2d entrantPrev;
	double previousPregnancyLag;
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
	array_2d naturalDeaths;
	array_3d netMigration;
	array_2d asfr;
	array_2d sexRatioAtBirth;
	const int HIVSTEPS_PER_YEAR;
	const int DT;
	array_3d cd4Progression;
	array_3d cd4InitialDist;
	array_3d cd4Mortality;
	array_3d incrrAge;
	double relinfectArt;
	double iota;
	std::vector<double> incrrSex;
	bool incidMod;
	int eppMod;
	int scaleCd4Mortality;
	std::vector<double> projectionSteps;
	double timeEpidemicStart;
	int everArtEligId;
	std::vector<int> artCd4EligId;
	std::vector<double> specPopPercentElig;
	std::vector<double> pregnantArtElig;
	double who34PercentElig;
	double prevalenceCurrent;
	// Store prevalence at last time step for r-trend model
	double prevalenceLast;
	std::vector<int> ageGroupsStart;
	std::vector<double> rSplineRVec;
	std::vector<double> rTrendBeta;
	double rTrendTStab;
	double rTrendR0;
	// Outputs
	array_1d entrantPrevOut;
	array_1d incidence15to49;
	array_1d prevalence15to49;
	array_1d incrate15To49;
	array_1d rVec;
	array_3d infections;

public:
	Model(State &modelState, std::vector<double> ageGroupsSp, std::vector<double> vertTLag, std::vector<double> paedSurvLag,
	      bool popAdjust, array_2d entrantPop, array_2d birthLag, array_2d cumSurvey,
	      array_2d cumNetMigr, double netMigrationHivProb, array_3d paedSurvCd4Distrib,
	      array_2d entrantArtCov, array_4d paedSurvArtCd4Distrib, array_3d survRate,
	      array_3d netMigr, array_2d asfRate, array_2d sexRatioBirth, int hivStepsPerYear,
	      array_3d cd4Prog, array_3d cd4InitDist, array_3d cd4Mort, array_3d incrrAges, int tArtStart,
	      double artRelinfect, int eppModel, int scaleCd4Mort, std::vector<double> projSteps,
	      std::vector<int> artCd4EligibleId, std::vector<double> specPopulationPercentElig,
	      std::vector<double> pregnantWomenArtElig, double who34PercElig)
		: state(modelState),
		  population(boost::extents[DISEASE_STATUS][SEXES][MODEL_AGES]),
		  previousPopulation(boost::extents[DISEASE_STATUS][SEXES][MODEL_AGES]),
		  artPopulation(boost::extents[SEXES][AGE_GROUPS][DISEASE_STATUS][CD4_STAGES]),
		  previousArtPopulation(boost::extents[SEXES][AGE_GROUPS][DISEASE_STATUS][CD4_STAGES]),
		  hivPop(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		  previousHivPop(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		  entrantPrev(boost::extents[PROJECTION_YEARS][SEXES]),
		  entrantPopulation(boost::extents[PROJECTION_YEARS][SEXES]),
		  birthsLag(boost::extents[PROJECTION_YEARS][SEXES]),
		  cumulativeSurvey(boost::extents[PROJECTION_YEARS][SEXES]),
		  cumulativeNetMigr(boost::extents[PROJECTION_YEARS][SEXES]),
		  paedSurvCd4Dist(boost::extents[PROJECTION_YEARS][SEXES][CD4_STAGES]),
		  entrantArtCoverage(boost::extents[PROJECTION_YEARS][SEXES]),
		  paedSurvArtCd4Dist(boost::extents[PROJECTION_YEARS][SEXES][CD4_STAGES][TREATMENT_STAGES]),
		  survivalRate(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  naturalDeaths(boost::extents[SEXES][MODEL_AGES]),
		  netMigration(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  asfr(boost::extents[PROJECTION_YEARS][FERT_AGES]),
		  sexRatioAtBirth(boost::extents[PROJECTION_YEARS][SEXES]),
		  HIVSTEPS_PER_YEAR(hivStepsPerYear),
		  DT(1.0 / HIVSTEPS_PER_YEAR),
		  cd4Progression(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES - 1]),
		  cd4InitialDist(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		  cd4Mortality(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		  incrrAge(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]),
		  ageGroupsStart(AGE_GROUPS, 0),
		  entrantPrevOut(boost::extents[PROJECTION_YEARS]),
		  incidence15to49(boost::extents[PROJECTION_YEARS]),
		  prevalence15to49(boost::extents[(PROJECTION_YEARS - 1) * hivStepsPerYear]),
		  incrate15To49(boost::extents[(PROJECTION_YEARS - 1) * hivStepsPerYear]),
		  rVec(boost::extents[(PROJECTION_YEARS - 1) * hivStepsPerYear]),
		  infections(boost::extents[PROJECTION_YEARS][SEXES][MODEL_AGES]) {

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
		cd4Progression = cd4Prog;
		cd4InitialDist = cd4InitDist;
		cd4Mortality = cd4Mort;
		incrrAge = incrrAges;
		relinfectArt = artRelinfect;
		eppMod = eppModel;
		scaleCd4Mortality = scaleCd4Mort;
		projectionSteps = projSteps;
		artCd4EligId = artCd4EligibleId;
		specPopPercentElig = specPopulationPercentElig;
		pregnantArtElig = pregnantWomenArtElig;
		who34PercentElig = who34PercElig;
		ageGroupsSpan = ageGroupsSp;
		incidMod = FALSE;
		prevalenceCurrent = 0.0;
		useEntrantPrev = FALSE;
		everArtEligId = CD4_STAGES;

		initialiseState();

		for (int ageGroup = 1; ageGroup < AGE_GROUPS; ageGroup++) {
			ageGroupsStart[ageGroup] = ageGroupsStart[ageGroup - 1] + ageGroupsSpan[ageGroup - 1];
		}

		// Prepare outputs
		prevalence15to49[0] = 0.0;

	};

	void initialiseState() {
		population = state.getPopulation();
		previousPopulation = state.getPreviousPopulation();
		artPopulation = state.getArtPopulation();
		previousArtPopulation = state.getPreviousArtPopulation();
		hivPop = state.getHivPopulation();
		previousHivPop = state.getPreviousHivPopulation();
		previousPregnancyLag = state.getPreviousPregnancyLag();
		naturalDeaths = state.getNaturalDeaths();
	};

	void updateModelState() {
		state.updatePopulation(population);
		state.updateArtPopulation(artPopulation);
		state.updateHivPopulation(hivPop);
		state.updateNaturalDeaths(naturalDeaths);
	};

	void setEntrantPrev(array_2d entPrev) {
		useEntrantPrev = TRUE;
		entrantPrev = entPrev;
	};

	void setIncrrSex(std::vector<double> incrrSexRatio) {
		incidMod = TRUE;
		incrrSex = incrrSexRatio;
	}

	void initialiseRSpline(std::vector<double> rSplineVec) {
		rSplineRVec = rSplineVec;
	};

	void initialiseRTrend(std::vector<double> beta, double tStab, double r0) {
		rTrendBeta = beta;
		rTrendTStab = tStab;
		rTrendR0 = r0;
	};

	void intialiseNonDirectIncid(double iotaModel, double tsEpidemicStart) {
		iota = iotaModel;
		timeEpidemicStart = tsEpidemicStart;
	};

	array_3d getPopulation() {
		return population;
	};

	void agePopulation(int t) {
		for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
			for (int sex = 0; sex < SEXES; sex++) {
				for (int age = 1; age < MODEL_AGES; age++) {
					// Rcpp::Rcout << "Looping over sex " << sex << " and age " << age << " \n";
					population[diseaseStatus][sex][age] = previousPopulation[diseaseStatus][sex][age - 1];
				}
				// People do not age out of final open age group (80+)
				population[diseaseStatus][sex][MODEL_AGES - 1] += previousPopulation[diseaseStatus][sex][MODEL_AGES - 1];
			}
		}

		double hivAgProb[SEXES][AGE_GROUPS];
		for (int sex = 0; sex < SEXES; sex++) {
			int a = 0;
			for (int ageGroup = 0; ageGroup < (AGE_GROUPS - 1); ageGroup++) {
				hivAgProb[sex][ageGroup] = 0;
				for (int i = 0; i < ageGroupsSpan[ageGroup]; i++) {
					hivAgProb[sex][ageGroup] += previousPopulation[HIVP][sex][a];
					a++;
				}
				hivAgProb[sex][ageGroup] = (hivAgProb[sex][ageGroup] > 0) ? previousPopulation[HIVP][sex][a - 1] / hivAgProb[sex][ageGroup] : 0;
			}
			// People do not age out of final open age group (80+)
			hivAgProb[sex][AGE_GROUPS - 1] = 0.0;
		}

		for (int sex = 0; sex < SEXES; sex++) {
			for (int ageGroup = 1; ageGroup < AGE_GROUPS; ageGroup++) {
				for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
					hivPop[sex][ageGroup][cd4Stage] = (1 - hivAgProb[sex][ageGroup]) * previousHivPop[sex][ageGroup][cd4Stage] + hivAgProb[sex][ageGroup - 1] * previousHivPop[sex][ageGroup - 1][cd4Stage];
					if (t > timeArtStart) {
						for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
							artPopulation[sex][ageGroup][cd4Stage][treatmentStage] = (1 - hivAgProb[sex][ageGroup]) * previousArtPopulation[sex][ageGroup][cd4Stage][treatmentStage] + hivAgProb[sex][ageGroup - 1] * previousArtPopulation[sex][ageGroup - 1][cd4Stage][treatmentStage];
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
				entrantPrevalence = previousPregnancyLag * vertTransLag[t - 1] * paedSurveyLag[t - 1];
			}

			if (populationAdjust) {
				population[HIVN][sex][0] = entrantPopulation[t - 1][sex] * (1.0 - entrantPrevalence);
				population[HIVP][sex][0] = entrantPopulation[t - 1][sex] * entrantPrevalence;
			} else {
				population[HIVN][sex][0] = birthsLag[t - 1][sex] * cumulativeSurvey[t - 1][sex] * (1.0 - entrantPrevalence / paedSurveyLag[t - 1]) + cumulativeNetMigr[t - 1][sex] * (1.0 - previousPregnancyLag * netMigrHivProb);
				population[HIVP][sex][0] = birthsLag[t - 1][sex] * cumulativeSurvey[t - 1][sex] * entrantPrevalence + cumulativeNetMigr[t - 1][sex] * entrantPrevalence;
			}

			double paedSurvPos = population[HIVP][sex][0];

			entrantPrevOut[t] = (population[HIVP][MALE][0] + population[HIVP][FEMALE][0]) /
			                    (population[HIVN][MALE][0] + population[HIVN][FEMALE][0] + population[HIVP][MALE][0] + population[HIVP][FEMALE][0]);

			for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
				hivPop[sex][0][cd4Stage] = (1 - hivAgProb[sex][0]) * previousHivPop[sex][0][cd4Stage] + paedSurvPos * paedSurvCd4Dist[t][sex][cd4Stage] * (1.0 - entrantArtCoverage[t][sex]);
				if (t > timeArtStart) {
					for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
						artPopulation[sex][0][cd4Stage][treatmentStage] = (1 - hivAgProb[sex][0]) * previousArtPopulation[sex][0][cd4Stage][treatmentStage];
						artPopulation[sex][0][cd4Stage][treatmentStage] += paedSurvPos * paedSurvArtCd4Dist[t][sex][cd4Stage][treatmentStage] * entrantArtCoverage[t][sex];
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

					hivPopAgeGroup += population[HIVP][sex][a];

					// non-HIV mortality
					double deathRate = 1.0 - survivalRate[t][sex][a];
					double hivNegDeaths = population[HIVN][sex][a] * deathRate;
					population[HIVN][sex][a] -= hivNegDeaths; // survival HIV- population
					double hivPosDeaths = population[HIVP][sex][a] * deathRate;
					deathsMigAgeGroup -= hivPosDeaths;
					population[HIVP][sex][a] -= hivPosDeaths;   // survival HIV+ population
					naturalDeaths[sex][a] = hivNegDeaths + hivPosDeaths;

					// net migration
					// TODO: rename migrate_a and hmig_a - what do these represent?
					double migrate_a = netMigration[t][sex][a] * (1 + survivalRate[t][sex][a]) / 2.0 / (population[HIVN][sex][a] + population[HIVP][sex][a]);
					population[HIVN][sex][a] *= 1 + migrate_a;
					double hmig_a = migrate_a * population[HIVP][sex][a];
					deathsMigAgeGroup += hmig_a;
					population[HIVP][sex][a] += hmig_a;

					a++;
				}

				// migration and deaths for hivpop
				double hivPopDeathsMigrRate = hivPopAgeGroup > 0 ? deathsMigAgeGroup / hivPopAgeGroup : 0.0;
				for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
					hivPop[sex][ageGroup][cd4Stage] *= 1 + hivPopDeathsMigrRate;
					if (t > timeArtStart) {
						for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
							artPopulation[sex][ageGroup][cd4Stage][treatmentStage] *= 1 + hivPopDeathsMigrRate;
						}
					}
				}
			}
		}
	};

	void fertility(int t) {
		Rcpp::Rcout << "Calculating fertility \n";
		double births = 0.0;
		std::vector<double> birthsByAgeGroup(FERT_AGE_GROUPS, 0.0);
		for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
			int a = pIDX_FERT;
			for (int ha = hIDX_FERT; ha < hIDX_FERT + FERT_AGE_GROUPS; ha++) {
				for (int i = 0; i < ageGroupsSpan[ha]; i++) {
					birthsByAgeGroup[ha - hIDX_FERT] += (previousPopulation[diseaseStatus][FEMALE][a] + population[diseaseStatus][FEMALE][a]) / 2 * asfr[t][a];
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

	void diseaseProgression(int t) {

		int cd4EligId = artCd4EligId[t] - 1; // -1 for 0-based indexing vs 1-based in R
		int anyEligId = ((specPopPercentElig[t] > 0) | (pregnantArtElig[t] > 0)) ? 0 : (who34PercentElig > 0) ? hIDX_CD4_350 : cd4EligId;
		everArtEligId = anyEligId < everArtEligId ? anyEligId : everArtEligId;

		for (int hivStep = 0; hivStep < HIVSTEPS_PER_YEAR; hivStep++) {
			int ts = (t - 1) * HIVSTEPS_PER_YEAR + hivStep;
			double hivDeathsByAgeGroup[SEXES][AGE_GROUPS];
			double grad[SEXES][AGE_GROUPS][CD4_STAGES];
			for (int sex = 0; sex < SEXES; sex++) {
				for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
					for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {

						double cd4mxScale = 1.0;
						if (scaleCd4Mortality & (t >= timeArtStart) & (cd4Stage >= everArtEligId)) {
							double artPopAgeSex = 0.0;
							for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
								artPopAgeSex += artPopulation[sex][ageGroup][cd4Stage][treatmentStage];
							}
							cd4mxScale = hivPop[sex][ageGroup][cd4Stage] / (hivPop[sex][ageGroup][cd4Stage] + artPopAgeSex);
						}

						double deaths = cd4mxScale * cd4Mortality[sex][ageGroup][cd4Stage] * hivPop[sex][ageGroup][cd4Stage];
						hivDeathsByAgeGroup[sex][ageGroup] += DT * deaths;
						grad[sex][ageGroup][cd4Stage] = -deaths;
					}
					for (int cd4Stage = 1; cd4Stage < CD4_STAGES; cd4Stage++) {
						grad[sex][ageGroup][cd4Stage - 1] -= cd4Progression[sex][ageGroup][cd4Stage - 1] * hivPop[sex][ageGroup][cd4Stage - 1];
						grad[sex][ageGroup][cd4Stage] += cd4Progression[sex][ageGroup][cd4Stage - 1] * hivPop[sex][ageGroup][cd4Stage - 1];
					}
				}
			}

			if (eppMod != EPP_DIRECTINCID) {
				// incidence

				// calculate r(t)
				if (eppMod == EPP_RSPLINE) {
					rVec[ts] = rSplineRVec[ts];
				} else {
					rVec[ts] = calcRtrendRt(hivStep, ts);
				}

				// calculate new infections by sex and age
				double infectionsBySexAge[SEXES][MODEL_AGES];
				if (incidMod) {
					calcInfectionsEppSpectrum((projectionSteps[ts] == timeEpidemicStart) ? iota : 0.0,
					                          t, hivStep, ts, infectionsBySexAge);
				}

				prevalence15to49[ts] = prevalenceCurrent;

				// add new infections to HIV population
				for (int sex = 0; sex < SEXES; sex++) {
					int a = 0;
					for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
						double infectionsPerAgeGroup = 0.0;
						for (int i = 0; i < ageGroupsSpan[ageGroup]; i++) {
							infectionsPerAgeGroup += infectionsBySexAge[sex][a];
							infections[t][sex][a] += DT * infectionsBySexAge[sex][a];
							population[HIVN][sex][a] -= DT * infectionsBySexAge[sex][a];
							population[HIVP][sex][a] += DT * infectionsBySexAge[sex][a];
							a++;
						}
						if (ageGroup < hIDX_15TO49 + hAG_15TO49 ) {
							incidence15to49[t] += DT * infectionsPerAgeGroup;
						}

						// add infections to grad hivpop
						for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
							grad[sex][ageGroup][cd4Stage] += infectionsPerAgeGroup * cd4InitialDist[sex][ageGroup][cd4Stage];
						}
					}
				}
			}
		}
	}

	double calcRtrendRt(int hivStep, int ts) {
		// sum population sizes
		double Xhivn = 0.0, Xhivp = 0.0;
		for (int sex = 0; sex < SEXES; sex++)
			for (int a = pIDX_15TO49; a < pIDX_15TO49 + pAG_15TO49; a++) {
				Xhivn += population[HIVN][sex][a];
				Xhivp += population[HIVP][sex][a];
			}

		// adjust HIV population for partial year time step
		for (int sex = 0; sex < SEXES; sex++) {
			Xhivn -= population[HIVN][sex][pIDX_15TO49] * (1.0 - DT * hivStep);
			Xhivp -= population[HIVP][sex][pIDX_15TO49] * (1.0 - DT * hivStep);
			Xhivn += population[HIVN][sex][pIDX_15TO49 + pAG_15TO49] * (1.0 - DT * hivStep);
			Xhivp += population[HIVP][sex][pIDX_15TO49 + pAG_15TO49] * (1.0 - DT * hivStep);
		}

		double Xtot = Xhivn + Xhivp;

		prevalenceLast = prevalenceCurrent;
		prevalenceCurrent = Xhivp / Xtot;

		// calculate r(t)
		double projStep = projectionSteps[ts];
		if (projStep > timeEpidemicStart) {
			double rVecLast = rVec[ts - 1];
			double gammaTs = (projStep < rTrendTStab) ? 0.0 : (prevalenceCurrent - prevalenceLast) * (projStep - rTrendTStab) / (DT * (prevalenceLast));
			double logrDiff = rTrendBeta[1] * (rTrendBeta[0] - rVecLast) + rTrendBeta[2] * (prevalenceLast) + rTrendBeta[3] * gammaTs;
			return exp(log(rVecLast) + logrDiff);
		} else {
			return rTrendR0;
		}
	}

	void calcInfectionsEppSpectrum(double iota, int t, int hts, int ts, double infectionsBySexAge[SEXES][MODEL_AGES]) {

		// sum population sizes
		double Xhivn_g[SEXES], Xhivn_incagerr[SEXES], Xhivp_noart = 0.0, Xart = 0.0;
		for (int sex = 0; sex < SEXES; sex++) {
			Xhivn_g[sex] = 0.0;
			Xhivn_incagerr[sex] = 0.0;
			for (int a = pIDX_15TO49; a < pIDX_15TO49 + pAG_15TO49; a++) {
				Xhivn_g[sex] += population[HIVN][sex][a];
				Xhivn_incagerr[sex] += incrrAge[t][sex][a] * population[HIVN][sex][a];
			}

			for (int ha = hIDX_15TO49; ha < hIDX_15TO49 + hAG_15TO49 + 1; ha++) {

				// adjustment to first and last age group for partial year time step
				// calculation proportion of HIV population to include / exclude based on hivpop in single-year ages.
				double prop_include;
				if (ha == hIDX_15TO49) {
					double hivp_ha = 0.0;
					for (int a = ageGroupsStart[ha]; a < ageGroupsStart[ha] + ageGroupsSpan[ha]; a++)
						hivp_ha += population[HIVP][sex][a];
					prop_include = (hivp_ha > 0) ? 1.0 - population[HIVP][sex][ageGroupsStart[ha]] / hivp_ha * (1.0 - DT * hts) : 1.0;
				} else if (ha == hIDX_15TO49 + hAG_15TO49) {
					double hivp_ha = 0.0;
					for (int a = ageGroupsStart[ha]; a < ageGroupsStart[ha] + ageGroupsSpan[ha]; a++)
						hivp_ha += population[HIVP][sex][a];
					prop_include = (hivp_ha > 0) ? population[HIVP][sex][ageGroupsStart[ha]] / hivp_ha * (1.0 - DT * hts) : 1.0;
				} else
					prop_include = 1.0;

				for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
					Xhivp_noart += hivPop[sex][ha][cd4Stage] * prop_include;
					if (t >= timeArtStart)
						for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++)
							Xart += artPopulation[sex][ha][cd4Stage][treatmentStage] * prop_include;
				}
			}
		}
		double Xhivn = Xhivn_g[MALE] + Xhivn_g[FEMALE];

		// adjust HIV negative population for partial year time step
		for (int sex = 0; sex < SEXES; sex++) {
			Xhivn -= population[HIVN][sex][pIDX_15TO49] * (1.0 - DT * hts);
			Xhivn += population[HIVN][sex][pIDX_15TO49 + pAG_15TO49] * (1.0 - DT * hts);
		}

		double Xtot = Xhivn + Xhivp_noart + Xart;
		prevalenceCurrent = (Xhivp_noart + Xart) / Xtot;

		incrate15To49[ts] = rVec[ts] * (Xhivp_noart + relinfectArt * Xart) / Xtot + iota;

		// incidence by sex
		double incrate15To49Sex[SEXES];
		incrate15To49Sex[MALE] = incrate15To49[ts] * (Xhivn_g[MALE] + Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + incrrSex[t] * Xhivn_g[FEMALE]);
		incrate15To49Sex[FEMALE] = incrate15To49[ts] * incrrSex[t] * (Xhivn_g[MALE] + Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + incrrSex[t] * Xhivn_g[FEMALE]);

		// annualized infections by age and sex
		for (int sex = 0; sex < SEXES; sex++)
			for (int age = 0; age < MODEL_AGES; age++) {
				infectionsBySexAge[sex][age] = population[HIVN][sex][age] * incrate15To49Sex[sex] * incrrAge[t][sex][age] * Xhivn_g[sex] / Xhivn_incagerr[sex];
			}

		return;
	}
};

#endif