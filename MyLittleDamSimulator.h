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
