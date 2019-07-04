#ifndef _CPPTEST_DISEASEPROGRESSION_H_
#define _CPPTEST_DISEASEPROGRESSION_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "consts.h"
#include "AgeMethod.h"
#include "state.h"
#include "ArtData.h"
#include "Cd4Data.h"
#include "InfectionData.h"
#include "PopulationData.h"
#include "HivData.h"
#include "Metadata.h"

class DiseaseProgressionMethod {
public:
	virtual void diseaseProgression();
}

class DefaultDiseaseProgressionMethod: public DiseaseProgressionMethod {
	public diseaseProgression() {
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
		}
	}
}

class RVecMethod() {
public:
	virtual void rVecMethod(int hivStep, int ts);
}


class RSplineRVecMethod: public RVecMethod {
	public rVecMethod(int hivStep, int ts) {
		rVec[ts] = rSplineRVec[ts];
	}
}

class DefaultRVecMethod: public RVecMethod {
	public rVecMethod(int hivStep, int ts) {
		// sum population sizes
		double Xhivn = 0.0, Xhivp = 0.0;
		for (int sex = 0; sex < SEXES; sex++) {
			for (int a = pIDX_15TO49; a < pIDX_15TO49 + pAG_15TO49; a++) {
				Xhivn += state.population[HIVN][sex][a];
				Xhivp += state.population[HIVP][sex][a];
			}
		}

		// adjust HIV population for partial year time step
		for (int sex = 0; sex < SEXES; sex++) {
			Xhivn -= state.population[HIVN][sex][pIDX_15TO49] * (1.0 - DT * hivStep);
			Xhivp -= state.population[HIVP][sex][pIDX_15TO49] * (1.0 - DT * hivStep);
			Xhivn += state.population[HIVN][sex][pIDX_15TO49 + pAG_15TO49] * (1.0 - DT * hivStep);
			Xhivp += state.population[HIVP][sex][pIDX_15TO49 + pAG_15TO49] * (1.0 - DT * hivStep);
		}

		double Xtot = Xhivn + Xhivp;

		state.prevalence = Xhivp / Xtot;

		// calculate r(t)
		double projStep = meta.projectionSteps[ts];
		if (projStep > timeEpidemicStart) {
			double rVecLast = rVec[ts - 1];
			double gammaTs = (projStep < rTrendTStab) ? 0.0 : (state.prevalence - state.previousPrevalence) * (projStep - rTrendTStab) / (DT * (state.previousPrevalence));
			double logrDiff = rTrendBeta[1] * (rTrendBeta[0] - rVecLast) + rTrendBeta[2] * (state.previousPrevalence) + rTrendBeta[3] * gammaTs;
			return exp(log(rVecLast) + logrDiff);
		} else {
			return rTrendR0;
		}
	}
}

class DirectIncidDiseaseProgressionMethod: public DefaultDiseaseProgressionMethod,
	public RVecMethod {
public:
	diseaseProgression() {
		DefaultDiseaseProgressionMethod::diseaseProgression();
		RVecMethod::rVecMethod();
		// calculate new infections by sex and age
		double infectionsBySexAge[SEXES][MODEL_AGES];
		if (infection.incidMod) {
			calcInfectionsEppSpectrum((meta.projectionSteps[ts] == timeEpidemicStart) ? iota : 0.0,
			                          t, hivStep, ts, infectionsBySexAge);
		}

		prevalence15to49[ts] = state.prevalence;

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

#endif