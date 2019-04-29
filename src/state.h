#ifndef _CPPTEST_STATE_H_
#define _CPPTEST_STATE_H_

#include <Rcpp.h>
#include <boost/multi_array.hpp>
#include "read_array.h"
#include "consts.h"

class State {
private:
	array_3d population;
	array_3d previousPopulation;

public:
	State(array_2d basePopulation) :
		population(boost::extents[DISEASE_STATUS][SEXES][MODEL_AGES]),
		previousPopulation(boost::extents[DISEASE_STATUS][SEXES][MODEL_AGES]) {

		array_3d_view_2d hivNegPopulation = population[ boost::indices[HIVN][range(0, SEXES)][range(0, MODEL_AGES)] ];
		array_3d_view_2d hivPosPopulation = population[ boost::indices[HIVP][range(0, SEXES)][range(0, MODEL_AGES)] ];
		hivNegPopulation = basePopulation;

		// TODO: Can this be done outside of the for loops? Tried using fill_n like we do in read_array readEntrantArtCoverage
		// but it looks like this approach isn't valid for array view.
		for (int sex = 0; sex < SEXES; sex++) {
			for (int age = 0; age < MODEL_AGES; age++) {
				hivPosPopulation[sex][age] = 0.0;
			}
		}
	}

	array_3d getPopulation() {
		return population;
	}

	array_3d getPreviousPopulation() {
		return previousPopulation;
	}

	void updatePopulation(array_3d newPopulation) {
		previousPopulation = population;
		population = newPopulation;
	}
}