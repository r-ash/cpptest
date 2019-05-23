#ifndef _CPPTEST_STATE_H_
#define _CPPTEST_STATE_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "read_array.h"
#include "consts.h"

typedef array_3d::array_view<2>::type array_3d_view_2d;
typedef array_4d::array_view<3>::type array_4d_view_3d;
typedef boost::multi_array_types::index_range range;

class State {
public:
	array_3d population;
	array_3d previousPopulation;
	array_4d outputPopulation;
	array_4d artPopulation;
	array_4d previousArtPopulation;
	array_3d hivPop;
	array_3d previousHivPop;
	array_2d naturalDeaths;
	array_2d infections;
	double prevalence;
	double previousPrevalence;

	State(array_2d basePopulation) :
		population(boost::extents[DISEASE_STATUS][SEXES][MODEL_AGES]),
		previousPopulation(boost::extents[DISEASE_STATUS][SEXES][MODEL_AGES]),
		outputPopulation(boost::extents[PROJECTION_YEARS][DISEASE_STATUS][SEXES][MODEL_AGES]),
		artPopulation(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES][TREATMENT_STAGES]),
		previousArtPopulation(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES][TREATMENT_STAGES]),
		hivPop(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		previousHivPop(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		naturalDeaths(boost::extents[SEXES][MODEL_AGES]),
		infections(boost::extents[SEXES][MODEL_AGES]) {

		array_3d_view_2d hivNegPopulation = population[ boost::indices[HIVN][range(0, SEXES)][range(0, MODEL_AGES)] ];
		array_3d_view_2d hivPosPopulation = population[ boost::indices[HIVP][range(0, SEXES)][range(0, MODEL_AGES)] ];
		hivNegPopulation = basePopulation;

		// TODO: Can this be done outside of the for loops? Tried using fill_n like we do in read_array readEntrantArtCoverage
		// but it looks like this approach isn't valid for array view.
		for (int sex = 0; sex < SEXES; sex++) {
			for (int age = 0; age < MODEL_AGES; age++) {
				hivPosPopulation[sex][age] = 0.0;
				naturalDeaths[sex][age] = 0.0;
			}
		}
		previousPopulation = population;
		array_4d_view_3d firstYearPopulation = outputPopulation[ boost::indices[0][range(0, DISEASE_STATUS)][range(0, SEXES)][range(0, MODEL_AGES)] ];
		firstYearPopulation = population;

		for (int sex = 0; sex < SEXES; sex++) {
			for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
				for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
					hivPop[sex][ageGroup][cd4Stage] = 0.0;
					for (int treatmentStage = 0; treatmentStage < TREATMENT_STAGES; treatmentStage++) {
						// Initialise to 0 in year of ART start
						artPopulation[sex][ageGroup][cd4Stage][treatmentStage] = 0.0;
					}
				}
			}
		}
		previousArtPopulation = artPopulation;
		previousHivPop = hivPop;
	}

	array_2d getNaturalDeaths() {
		return naturalDeaths;
	}

	void updatePopulation(int t) {
		Rcpp::Rcout << "Updating population with time t = " << t << "\n";
		previousPopulation = population;
		array_4d_view_3d tPopulation = outputPopulation[ boost::indices[t][range(0, DISEASE_STATUS)][range(0, SEXES)][range(0, MODEL_AGES)] ];
		tPopulation = population;
	}

	void updateArtPopulation() {
		previousArtPopulation = artPopulation;
	}

	void updateHivPopulation() {
		previousHivPop = hivPop;
	}

	void updateNaturalDeaths() {
		// Record for output?
	}

	void updateInfections() {
		// Record for output?
	}
};

#endif