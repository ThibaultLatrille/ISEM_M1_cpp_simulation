/*
 * Simulations.h
 *
 *  Created on: Mar 13, 2014
 *      Author: Thibault
 */

#ifndef SIMULATIONS_H_
#define SIMULATIONS_H_

using namespace std;

string IntToStr(int);
string DoubleToStr( double);

int conditional_expectation( void );
int exponential_growth( void );
int simulate_infection(bool,bool,int);

#endif /* SIMULATIONS_H_ */
