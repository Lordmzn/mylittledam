/*
 * MyLittleDamSimulator.h
 *
 *  Created on: 26/ott/2012
 *      Author: lordmzn
 */

#ifndef MYLITTLEDAMSIMULATOR_H_
#define MYLITTLEDAMSIMULATOR_H_
#include <vector>

struct DamSimParam {
	std::vector<double>::size_type horizon;
	double s0;
	double water_demand;
	double liv_crit;
	double inflow;
	int RBF_N;
	int nvars;
	int nobjs;
};

class MyLittleDamSimulator {
public:
	MyLittleDamSimulator(DamSimParam p);
	void runSimulation(double vars[], double objs[]);

private:
	double RBFs_policy(std::vector<double> input,
			std::vector<double> parameters, unsigned int N);
	DamSimParam p;
};

#endif /* MYLITTLEDAMSIMULATOR_H_ */
