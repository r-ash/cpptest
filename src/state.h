#ifndef _CPPTEST_STATE_H_
#define _CPPTEST_STATE_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "read_array.h"
#include "consts.h"

typedef array_3d::array_view<2>::type array_3d_view_2d;
typedef boost::multi_array_types::index_range range;

class State {
private:
	array_3d population;
	array_3d previousPopulation;
	array_4d artPopulation;
	array_4d previousArtPopulation;
	array_4d yearArtStartPopulation;
	array_3d hivPop;
	array_3d previousHivPop;
	double previousPregnancyLag;
	array_2d naturalDeaths;

public:
	State(array_2d basePopulation) :
		population(boost::extents[DISEASE_STATUS][SEXES][MODEL_AGES]),
		previousPopulation(boost::extents[DISEASE_STATUS][SEXES][MODEL_AGES]),
		artPopulation(boost::extents[SEXES][AGE_GROUPS][DISEASE_STATUS][CD4_STAGES]),
		previousArtPopulation(boost::extents[SEXES][AGE_GROUPS][DISEASE_STATUS][CD4_STAGES]),
		hivPop(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		previousHivPop(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		naturalDeaths(boost::extents[SEXES][MODEL_AGES]) {

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

		for (int sex = 0; sex < SEXES; sex++) {
			for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
				for (int diseaseStatus = 0; diseaseStatus < DISEASE_STATUS; diseaseStatus++) {
					for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
						// Initialise to 0 in year of ART start
						artPopulation[sex][ageGroup][diseaseStatus][cd4Stage] = 0.0;
					}
				}
			}
		}
		previousArtPopulation = artPopulation;

		for (int sex = 0; sex < SEXES; sex++) {
			for (int ageGroup = 0; ageGroup < AGE_GROUPS; ageGroup++) {
				for (int cd4Stage = 0; cd4Stage < CD4_STAGES; cd4Stage++) {
					hivPop[sex][ageGroup][cd4Stage] = 0.0;
				}
			}
		}
		previousHivPop = hivPop;

		previousPregnancyLag = 0.0;
	}

	array_3d getPopulation() {
		return population;
	}

	array_3d getPreviousPopulation() {
		return previousPopulation;
	}

	array_4d getArtPopulation() {
		return artPopulation;
	}

	array_4d getPreviousArtPopulation() {
		return previousArtPopulation;
	}

	array_3d getHivPopulation() {
		return hivPop;
	}

	array_3d getPreviousHivPopulation() {
		return previousHivPop;
	}

	double getPreviousPregnancyLag() {
		return previousPregnancyLag;
	}

	array_2d getNaturalDeaths() {
		return naturalDeaths;
	}

	void updatePopulation(array_3d newPopulation) {
		previousPopulation = population;
		population = newPopulation;
	}

	void updateArtPopulation(array_4d newArtPopulation) {
		previousArtPopulation = artPopulation;
		artPopulation = newArtPopulation;
	}

	void updateHivPopulation(array_3d newHivPopulation) {
		previousHivPop = hivPop;
		hivPop = newHivPopulation;
	}

	void updateNaturalDeaths(array_2d newNaturalDeaths) {
		naturalDeaths = newNaturalDeaths;
	}
};

#endif