// By Thibault Latrille
// thibault.latrille@ens-lyon.fr
// 11/03/2014
// No copyright (Why on hell should there be one anyway!)

// This simulation intend verify our analytical results by a Monte Carlo algorithm, refer to the paper (also in the github repository) for more information.
// This piece of code is specifically aimed at testing the following formula (for 2 families):
// \mathbb{E}\left[\dfrac{1}{N_+(t)-1} \right]=\dfrac{e^{-\rho t}}{r_+-1}
// Write a file of the form r=100_tmax=20_rho=0.1_replicate=1000.txt in the current directory.
// The first line is the header.
// Each lines contains the estimation, the 95% confidence interval as long as the true value for one fixed time.
// Note that two estimation on different lines are not independent of one another.
// Must be compiled by cpp11 (-std=c++11)
// Must allocate dynamically the array for improved performance

#include <iostream>		// setbuf, NULL
#include <cstdio>		// stdout
#include <cstdlib>      // standard library
#include <cmath> 		// sqrt, pow
#include <random>
#include <algorithm>
#include <chrono>
#include <fstream>		// ofstream
#include <unistd.h>     // chdir
#include "Simulations.h"

int main() {

	//setbuf(stdout, NULL);
	// Specifically for eclipse IDE, don't buffer output

	chdir("/cygdrive/c/Users/Thibault/Desktop/ISEM_M1_Cpp_Simulations/simulated_data/expectation_time_14_03");
	// set current  directory

	//////////////////////////////////

	const double lambda(0.005);
	const double tau(0.2);

	const int number_simulation=1;
	const int number_nematode=2;
	const int r=100;

	double runi;

	double time_growth[number_nematode];
	double family_size[number_nematode];

	const bool stoch_time=true;
	const bool stoch_growth=true;

	fill(family_size, family_size+number_nematode, double(r) );

	vector<double> vector_replicate( number_simulation );

	////////////////////////////////////

	//unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	unsigned seed=123;
	mt19937 generator(seed);
	uniform_real_distribution<double> u_distribution(0.0,1.0);

	negative_binomial_distribution<int> nb_distribution(10000,0.5);

	cout << nb_distribution(generator) << endl;

	for (int ith_simu=0; ith_simu<number_simulation; ith_simu++){

		fill(family_size, family_size+number_nematode, double(r) );
		time_growth[number_nematode-1]=0.0;

		for (int ith_nematode=(number_nematode-1) ; ith_nematode>0 ; ith_nematode-- ){

			if (stoch_time)
				{
				runi = u_distribution(generator);
				time_growth[ith_nematode-1] = time_growth[ith_nematode]-log(runi)/(tau);
				}
				else time_growth[ith_nematode-1] = time_growth[ith_nematode]+1./(tau);
		}

		for (int ith_nematode=0 ; ith_nematode<number_nematode ; ith_nematode++ ){
			if (stoch_growth){
				double p=exp(-lambda*(time_growth[ith_nematode]));
				negative_binomial_distribution<int> nb_distribution(family_size[ith_nematode],p);
				family_size[ith_nematode]+=nb_distribution(generator);
			}
			else family_size[ith_nematode]=family_size[ith_nematode]*exp(lambda*time_growth[ith_nematode]);
		}

		double sum_r=0.;
		double sum_rsquared=0.;

		cout << "For the " << ith_simu << "th simulation:" << endl;
		for (int ith_nematode=0 ; ith_nematode<number_nematode ; ith_nematode++ ){
			cout << "for the " << ith_nematode << "th nematode, the size is " << family_size[ith_nematode] << " and growth time is" << time_growth[ith_nematode] <<endl;
			sum_r+=family_size[ith_nematode];
			sum_rsquared+=pow(family_size[ith_nematode],2);
		}

		vector_replicate[ith_simu]=sum_rsquared/pow(sum_r,2);

	}

	double s1;
	double mean;
	s1 = accumulate( vector_replicate.begin() , vector_replicate.end() , 0.0 );
	mean = s1 / number_simulation;

	double s2(0.0);
	double stdev;
	double upper_bound_stdev;
	double lower_bound_stdev;
	s2 = 0.0;
	for ( vector<double>::iterator it = vector_replicate.begin(); it != vector_replicate.end(); it++ ){
		s2+=pow(*it,2.0);
	}

	stdev = 1.96*sqrt( (s2/number_simulation-pow( mean,2.) ) / (number_simulation-1) );
	lower_bound_stdev=mean-stdev;
	upper_bound_stdev=mean+stdev;
	cout << "lower bound stdev = " << lower_bound_stdev <<endl;
	cout << "mean = " << mean <<endl;
	cout << "upper bound  stdev= " << upper_bound_stdev <<endl;

	double p=0.025;
	int rank=floor(p*number_simulation);
	sort (vector_replicate.begin() , vector_replicate.end());
	cout << "rank = " << rank <<endl;

	double lower_bound(vector_replicate[rank]);
	double upper_bound(vector_replicate[number_simulation-rank]);

	cout << "lower bound = " << lower_bound <<endl;
	cout << "mean = " << mean <<endl;
	cout << "upper bound = " << upper_bound <<endl;
	double approximat_relatedness;
	double q=exp(lambda/tau);

	cout << "q = " << q <<endl;

	approximat_relatedness=(q-1)*(pow(q,number_nematode)+1)/((q+1)*(pow(q,number_nematode)-1));

	cout << "approximat_relatedness = " << approximat_relatedness <<endl;
    return 0;
}

