#include "model.h"

void Model::agePopulation(int t) {
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
			state.population[HIVN][sex][0] = pop.birthsLag[t - 1][sex] * pop.cumulativeSurvey[t - 1][sex] * (1.0 - entrantPrevalence / hiv.paedSurveyLag[t - 1]) + pop.cumulativeNetMigr[t - 1][sex] * (1.0 - previousPregnancyLag[t - 1] * hiv.netMigrProb);
			state.population[HIVP][sex][0] = pop.birthsLag[t - 1][sex] * pop.cumulativeSurvey[t - 1][sex] * entrantPrevalence + pop.cumulativeNetMigr[t - 1][sex] * entrantPrevalence;
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

void Model::deathsAndMigration(int t) {
	// Non-HIV mortality
	for (int sex = 0; sex < SEXES; sex++) {
		int a = 0;
		for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
			double deathsMigAgeGroup = 0;
			double hivPopAgeGroup = 0;
			for (int i = 0; i < meta.ageGroupsSpan[ageGroup]; i++) {

				hivPopAgeGroup += state.population[HIVP][sex][a];

				// non-HIV mortality
				double deathRate = 1.0 - pop.survivalRate[t][sex][a];
				double hivNegDeaths = state.population[HIVN][sex][a] * deathRate;
				state.population[HIVN][sex][a] -= hivNegDeaths; // survival HIV- population
				double hivPosDeaths = state.population[HIVP][sex][a] * deathRate;
				deathsMigAgeGroup -= hivPosDeaths;
				state.population[HIVP][sex][a] -= hivPosDeaths;   // survival HIV+ population
				state.naturalDeaths[sex][a] = hivNegDeaths + hivPosDeaths;

				// net migration
				// TODO: rename migrate_a and hmig_a - what do these represent?
				double migrate_a = pop.netMigration[t][sex][a] * (1 + pop.survivalRate[t][sex][a]) / 2.0 / (state.population[HIVN][sex][a] + state.population[HIVP][sex][a]);
				state.population[HIVN][sex][a] *= 1 + migrate_a;
				double hmig_a = migrate_a * state.population[HIVP][sex][a];
				deathsMigAgeGroup += hmig_a;
				state.population[HIVP][sex][a] += hmig_a;

				a++;
			}

			// migration and deaths for hivpop
			double hivPopDeathsMigrRate = hivPopAgeGroup > 0 ? deathsMigAgeGroup / hivPopAgeGroup : 0.0;
			for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
				state.hivPop[sex][ageGroup][cd4Stage] *= 1 + hivPopDeathsMigrRate;
				if (t > meta.timeArtStart) {
					for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
						state.artPopulation[sex][ageGroup][cd4Stage][treatmentStage] *= 1 + hivPopDeathsMigrRate;
					}
				}
			}
		}
	}
}

void Model::fertility(int t) {
	Rcpp::Rcout << "Calculating fertility \n";
	double births = 0.0;
	std::vector<double> birthsByAgeGroup(FERT_AGE_GROUPS, 0.0);
	for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
		int a = pIDX_FERT;
		for (int ha = hIDX_FERT; ha < hIDX_FERT + FERT_AGE_GROUPS; ha++) {
			for (int i = 0; i < meta.ageGroupsSpan[ha]; i++) {
				birthsByAgeGroup[ha - hIDX_FERT] += (state.previousPopulation[diseaseStatus][FEMALE][a] + state.population[diseaseStatus][FEMALE][a]) / 2 * pop.asfr[t][a];
				a++;
			}
		}
	}
	for (int ha = hIDX_FERT; ha < FERT_AGE_GROUPS; ha++) {
		births += birthsByAgeGroup[ha - hIDX_FERT];
	}

	if (t + AGE_START < PROJECTION_YEARS) {
		for (int sex = 0; sex < SEXES; sex++) {
			pop.birthsLag[t + AGE_START - 1][sex] = pop.sexRatioAtBirth[t][sex] * births;
		}
	}
}

void Model::diseaseProgression(int t) {
	for (int hivStep = 0; hivStep < meta.hivStepsPerYear; hivStep++) {
		int ts = (t - 1) * meta.hivStepsPerYear + hivStep;
		double hivDeathsByAgeGroup[SEXES][AGE_GROUPS];
		double grad[SEXES][AGE_GROUPS][CD4_STAGES];
		for (int sex = 0; sex < SEXES; sex++) {
			for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
				for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {

					double cd4mxScale = 1.0;
					if (meta.scaleCd4Mortality & (t >= meta.timeArtStart) & (cd4Stage >= art.everArtEligId[t - 1])) {
						double artPopAgeSex = 0.0;
						for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
							artPopAgeSex += state.artPopulation[sex][ageGroup][cd4Stage][treatmentStage];
						}
						cd4mxScale = state.hivPop[sex][ageGroup][cd4Stage] / (state.hivPop[sex][ageGroup][cd4Stage] + artPopAgeSex);
					}

					double deaths = cd4mxScale * cd4.mortality[sex][ageGroup][cd4Stage] * state.hivPop[sex][ageGroup][cd4Stage];
					hivDeathsByAgeGroup[sex][ageGroup] += DT * deaths;
					grad[sex][ageGroup][cd4Stage] = -deaths;
				}
				for (int cd4Stage = 1; cd4Stage < CD4_STAGES; cd4Stage++) {
					grad[sex][ageGroup][cd4Stage - 1] -= cd4.progression[sex][ageGroup][cd4Stage - 1] * state.hivPop[sex][ageGroup][cd4Stage - 1];
					grad[sex][ageGroup][cd4Stage] += cd4.progression[sex][ageGroup][cd4Stage - 1] * state.hivPop[sex][ageGroup][cd4Stage - 1];
				}
			}
		}

		if (meta.eppMod != EPP_DIRECTINCID) {
			// incidence

			// calculate r(t)
			if (meta.eppMod == EPP_RSPLINE) {
				rVec[ts] = rSplineRVec[ts];
			} else {
				rVec[ts] = calcRtrendRt(hivStep, ts);
			}

			// calculate new infections by sex and age
			double infectionsBySexAge[SEXES][MODEL_AGES];
			if (infection.incidMod) {
				calcInfectionsEppSpectrum((meta.projectionSteps[ts] == timeEpidemicStart) ? iota : 0.0,
				                          t, hivStep, ts, infectionsBySexAge);
			}

			prevalence15to49[ts] = prevalenceCurrent;

			// add new infections to HIV population
			for (int sex = 0; sex < SEXES; sex++) {
				int a = 0;
				for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
					double infectionsPerAgeGroup = 0.0;
					for (int i = 0; i < meta.ageGroupsSpan[ageGroup]; i++) {
						infectionsPerAgeGroup += infectionsBySexAge[sex][a];
						infections[sex][a] += DT * infectionsBySexAge[sex][a];
						state.population[HIVN][sex][a] -= DT * infectionsBySexAge[sex][a];
						state.population[HIVP][sex][a] += DT * infectionsBySexAge[sex][a];
						a++;
					}
					if (ageGroup < hIDX_15TO49 + hAG_15TO49 ) {
						incidence15to49[t] += DT * infectionsPerAgeGroup;
					}

					// add infections to grad hivpop
					for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
						grad[sex][ageGroup][cd4Stage] += infectionsPerAgeGroup * cd4.initialDist[sex][ageGroup][cd4Stage];
					}
				}
			}
		}
	}
}

double Model::calcRtrendRt(int hivStep, int ts) {
	// sum population sizes
	double Xhivn = 0.0, Xhivp = 0.0;
	for (int sex = 0; sex < SEXES; sex++)
		for (int a = pIDX_15TO49; a < pIDX_15TO49 + pAG_15TO49; a++) {
			Xhivn += state.population[HIVN][sex][a];
			Xhivp += state.population[HIVP][sex][a];
		}

	// adjust HIV population for partial year time step
	for (int sex = 0; sex < SEXES; sex++) {
		Xhivn -= state.population[HIVN][sex][pIDX_15TO49] * (1.0 - DT * hivStep);
		Xhivp -= state.population[HIVP][sex][pIDX_15TO49] * (1.0 - DT * hivStep);
		Xhivn += state.population[HIVN][sex][pIDX_15TO49 + pAG_15TO49] * (1.0 - DT * hivStep);
		Xhivp += state.population[HIVP][sex][pIDX_15TO49 + pAG_15TO49] * (1.0 - DT * hivStep);
	}

	double Xtot = Xhivn + Xhivp;

	prevalenceLast = prevalenceCurrent;
	prevalenceCurrent = Xhivp / Xtot;

	// calculate r(t)
	double projStep = meta.projectionSteps[ts];
	if (projStep > timeEpidemicStart) {
		double rVecLast = rVec[ts - 1];
		double gammaTs = (projStep < rTrendTStab) ? 0.0 : (prevalenceCurrent - prevalenceLast) * (projStep - rTrendTStab) / (DT * (prevalenceLast));
		double logrDiff = rTrendBeta[1] * (rTrendBeta[0] - rVecLast) + rTrendBeta[2] * (prevalenceLast) + rTrendBeta[3] * gammaTs;
		return exp(log(rVecLast) + logrDiff);
	} else {
		return rTrendR0;
	}
}

void Model::calcInfectionsEppSpectrum(double iota, int t, int hts, int ts, double infectionsBySexAge[SEXES][MODEL_AGES]) {

	// sum population sizes
	double Xhivn_g[SEXES], Xhivn_incagerr[SEXES], Xhivp_noart = 0.0, Xart = 0.0;
	for (int sex = 0; sex < SEXES; sex++) {
		Xhivn_g[sex] = 0.0;
		Xhivn_incagerr[sex] = 0.0;
		for (int a = pIDX_15TO49; a < pIDX_15TO49 + pAG_15TO49; a++) {
			Xhivn_g[sex] += state.population[HIVN][sex][a];
			Xhivn_incagerr[sex] += infection.incrrAge[t][sex][a] * state.population[HIVN][sex][a];
		}

		for (int ha = hIDX_15TO49; ha < hIDX_15TO49 + hAG_15TO49 + 1; ha++) {

			// adjustment to first and last age group for partial year time step
			// calculation proportion of HIV population to include / exclude based on hivpop in single-year ages.
			double prop_include;
			if (ha == hIDX_15TO49) {
				double hivp_ha = 0.0;
				for (int a = ageGroupsStart[ha]; a < ageGroupsStart[ha] + meta.ageGroupsSpan[ha]; a++)
					hivp_ha += state.population[HIVP][sex][a];
				prop_include = (hivp_ha > 0) ? 1.0 - state.population[HIVP][sex][ageGroupsStart[ha]] / hivp_ha * (1.0 - DT * hts) : 1.0;
			} else if (ha == hIDX_15TO49 + hAG_15TO49) {
				double hivp_ha = 0.0;
				for (int a = ageGroupsStart[ha]; a < ageGroupsStart[ha] + meta.ageGroupsSpan[ha]; a++)
					hivp_ha += state.population[HIVP][sex][a];
				prop_include = (hivp_ha > 0) ? state.population[HIVP][sex][ageGroupsStart[ha]] / hivp_ha * (1.0 - DT * hts) : 1.0;
			} else
				prop_include = 1.0;

			for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
				Xhivp_noart += state.hivPop[sex][ha][cd4Stage] * prop_include;
				if (t >= meta.timeArtStart)
					for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++)
						Xart += state.artPopulation[sex][ha][cd4Stage][treatmentStage] * prop_include;
			}
		}
	}
	double Xhivn = Xhivn_g[MALE] + Xhivn_g[FEMALE];

	// adjust HIV negative population for partial year time step
	for (int sex = 0; sex < SEXES; sex++) {
		Xhivn -= state.population[HIVN][sex][pIDX_15TO49] * (1.0 - DT * hts);
		Xhivn += state.population[HIVN][sex][pIDX_15TO49 + pAG_15TO49] * (1.0 - DT * hts);
	}

	double Xtot = Xhivn + Xhivp_noart + Xart;
	prevalenceCurrent = (Xhivp_noart + Xart) / Xtot;

	incrate15To49[ts] = rVec[ts] * (Xhivp_noart + meta.relinfectArt * Xart) / Xtot + iota;

	// incidence by sex
	double incrate15To49Sex[SEXES];
	incrate15To49Sex[MALE] = incrate15To49[ts] * (Xhivn_g[MALE] + Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + infection.incrrSex[t] * Xhivn_g[FEMALE]);
	incrate15To49Sex[FEMALE] = incrate15To49[ts] * infection.incrrSex[t] * (Xhivn_g[MALE] + Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + infection.incrrSex[t] * Xhivn_g[FEMALE]);

	// annualized infections by age and sex
	for (int sex = 0; sex < SEXES; sex++)
		for (int age = 0; age < MODEL_AGES; age++) {
			infectionsBySexAge[sex][age] = state.population[HIVN][sex][age] * incrate15To49Sex[sex] * infection.incrrAge[t][sex][age] * Xhivn_g[sex] / Xhivn_incagerr[sex];
		}

	return;
}