/*
 * MIT License
 * 
 * Copyright (c) 2016 Emanuele Mason
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
 * MyLittleDam.cpp
 *
 *  Created on: 26/ott/2012
 *      Author: lordmzn
 */


#include "MyLittleDamSimulator.h"
#include "moeaframework.h"
#include <cmath>

MyLittleDamSimulator::MyLittleDamSimulator(DamSimParam params) :
	p(params) {
	// ciaone
}

double MyLittleDamSimulator::RBFs_policy(std::vector<double> input,
		std::vector<double> parameters, unsigned int N) {
	// setting
	unsigned int D = input.size(); // number of states
	const double u_min = 0.0;
	const double u_max = 155;
	// reshape of parameters
	double c[N][D]; // center = (N-RBFs, D-state variables)
	double b[N][D]; // radius = (N-RBFs, D-state variables)
	unsigned int count = 0;
	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < D; j++) {
			c[i][j] = parameters[count];
			b[i][j] = parameters[count + 1];
			count = count + 2;
		}
	}
	std::vector<double> w; // weights of the RBFs
	for (unsigned int i = 0; i < N; i++) {
		w.push_back(parameters[count]);
		count = count + 1;
	}
	double bias = 0.0;
	bias = parameters[count];

	// Normalization of the input
	double input_norm[D][2];
	input_norm[0][0] = 0;
	input_norm[0][1] = 155;
	//input_norm[1][0] = 0 ;
	//input_norm[1][1] = 200 ;
	std::vector<double> input1;
	for (unsigned int i = 0; i < D; i++) {
		input1.push_back(
				(input[i] - input_norm[i][0])
						/ (input_norm[i][1] - input_norm[i][0]));
	}
	// RBFs
	std::vector<double> phi;
	double bf = 0;
	for (unsigned int j = 0; j < N; j++) {
		bf = 0;
		for (unsigned int i = 0; i < D; i++) {
			bf = bf
					+ (input1[i] - c[j][i]) * (input1[i] - c[j][i])
							/ (b[j][i] * b[j][i]);
		}
		phi.push_back(exp(-bf));
	}
	// compute decision
	double out = bias;
	for (unsigned int i = 0; i < N; i++) {
		out = out + w[i] * phi[i];
	}
	return out * (u_max - u_min) + u_min;
}

void MyLittleDamSimulator::runSimulation(double vars[], double objs[]) {
	std::vector<double> storage(p.horizon + 1, 0.0);
	std::vector<double> irr_stepcost(p.horizon + 1, 0.0);
	std::vector<double> flo_stepcost(p.horizon + 1, 0.0);
	std::vector<double> hyd_stepcost(p.horizon + 1, 0.0);
	std::vector<double> rel_stepcost(p.horizon + 1, 0.0);
	std::vector<double> rbf_params(p.nvars, 0.0);
	storage[0] = p.s0;
	double ril_SUP = 0.0;
	double ril_INF = 0.0;
	double release;
	double decision = 0.0;
	double hp;
	for (int i = 0; i < p.nvars; i++) {
		rbf_params[i] = vars[i];
	}

	for (std::vector<double>::size_type t = 0; t < p.horizon; t++) {
		decision = RBFs_policy(std::vector<double>(1, storage[t]), rbf_params, p.RBF_N);
		ril_SUP = storage[t];
		ril_INF = std::max(storage[t] - 100, 0.0);
		release = std::min(ril_SUP, std::max(ril_INF, decision));
		storage[t + 1] = storage[t] + p.inflow - release;
		irr_stepcost[t + 1] = std::max(50 - release, 0.0);
		flo_stepcost[t + 1] = std::max(storage[t + 1] - 50, 0.0);
		hyd_stepcost[t + 1] = std::max(4.36 - hp, 0.0);
		hp = 1 * 9.81 * 1000 / 3600000 * storage[t + 1] * std::max(release - 0.0, 0.0);
		rel_stepcost[t + 1] = std::max(release - 30.0, 0.0);
	}

	objs[0] = 0.0;
	objs[1] = 0.0;
	objs[2] = 0.0;
	objs[3] = 0.0;
	for (int t = 0; t < (int) irr_stepcost.size(); t++) {
		objs[0] = objs[0] + irr_stepcost[t];
		objs[1] = objs[1] + flo_stepcost[t];
		objs[2] = objs[2] + hyd_stepcost[t];
		objs[3] = objs[3] + rel_stepcost[t];
	}
	objs[0] = objs[0] / (double) irr_stepcost.size();
	objs[1] = objs[1] / (double) flo_stepcost.size();
	objs[2] = objs[2] / (double) hyd_stepcost.size();
	objs[3] = objs[3] / (double) rel_stepcost.size();
}

int main(int argc, char **argv) {
	DamSimParam param;
	param.horizon = 2401;
	param.s0 = 100.0;
	param.inflow = 40.0;
	param.RBF_N = 3;
	param.nvars = 10;
	param.nobjs = 4;
	double vars[param.nvars];
	double objs[param.nobjs];

	MyLittleDamSimulator dam(param);

	MOEA_Init(param.nobjs, 0);

	while (MOEA_Next_solution() == MOEA_SUCCESS) {
		MOEA_Read_doubles(param.nvars, vars);
		dam.runSimulation(vars, objs);
		MOEA_Write(objs, NULL);
	}

	MOEA_Terminate();

	return MOEA_SUCCESS;
}

