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
#include <chrono>
#include <fstream>		// ofstream
#include <unistd.h>     // chdir
#include "Simulations.h"

int main() {

	setbuf(stdout, NULL); // Specifically for eclipse IDE, don't buffer output
	chdir("/cygdrive/c/Users/Thibault/Desktop/ISEM_M1_Cpp_Simulations/simulated_data/expectation_time_14_03");
	const double rho=0.2;
	double tmax=20.;
	const int r=100;

	unsigned seed = chrono::system_clock::now().time_since_epoch().count(); // The seed
	mt19937 generator(seed); // The mersenne_twister generator
	uniform_real_distribution<double> distribution(0.0,1.0); // Uniform [0,1] distribution

	const int number_of_pop = 1000;
	// The number of populations we are computing, influence the range of the confidence interval
	const int number_of_time_steps = 250;
	// The number of time steps we are intersted in

	double processus[number_of_time_steps][number_of_pop];
	// The array containing all trajectories at each time steps.

	double output_array[number_of_time_steps][5]; // The array containing all relevants informations

	for (int index_time=0; index_time<number_of_time_steps; index_time++) // Initialize time steps.
	{
		output_array[index_time][0]=index_time*tmax/(number_of_time_steps-1);
	}

	for (int pop_index=0 ; pop_index<number_of_pop ; pop_index++ ) // For each independent simulation
	{
		int n=r;
		double t=0.0;
		double runi;
		double T;
		int index_time=0;
		double time_step=output_array[index_time][0];
		while (t<tmax) // Until tmax is not reached
		{
			runi = distribution(generator);
			T = -log(runi)/(rho*n); // T is exponentialy distributed
			t+=T;
			while (t>time_step)  // Write 1/(population size-1) at each time step
			{
				processus[index_time][pop_index]=1.0/(double(n)-1);
				index_time++;
				if (index_time==number_of_time_steps)
					break;
				time_step=output_array[index_time][0];
			}
			n++;
		}
	}

	string file_name = "r="+IntToStr(r)+"_tmax="+IntToStr(tmax)+"_rho="+DoubleToStr(rho)+"_replicate="+IntToStr(number_of_pop)+".txt";
	// The file name
	cout << file_name << endl;

	ofstream myfile;
	myfile.open( file_name.c_str() ); // Open the file
	myfile << "//time estimation lower_bound upper_bound expectation" << endl;
	// The first line of the file, the header

	for (int index_time=0; index_time<number_of_time_steps; index_time++)
	{
		double s1=0.; // The sum of 1/(population size-1) size at each time step
		for (int pop_index=0 ; pop_index<number_of_pop ; pop_index++ )
		{
			s1+=processus[index_time][pop_index];
		}
		output_array[index_time][1]=s1/number_of_pop;

		double s2=0.; // The sum of square of 1/(population size-1) at each time step
		for (int pop_index=0 ; pop_index<number_of_pop ; pop_index++ )
		{
			s2+=pow(processus[index_time][pop_index],2.0);
		}

		double stdev = 1.96*sqrt( (s2/number_of_pop-pow(output_array[index_time][1],2.0) ) / (number_of_pop-1) );
		// The standard deviation of 1/(population size-1) at each time step

		output_array[index_time][2] = output_array[index_time][1]-stdev; // Upper bound for the 95% confidence interval
		output_array[index_time][3] = output_array[index_time][1]+stdev; // Lower bound for the 95% confidence interval

		output_array[index_time][4]=exp(-rho*output_array[index_time][0])/(double (r)-1);
		// The expected value of 1/(population size-1)
		myfile << output_array[index_time][0] << " " << output_array[index_time][1] << " " << output_array[index_time][2] << " " << output_array[index_time][3] << " " << output_array[index_time][4] << endl;
	}
	myfile.close();
	cout << "work finish" <<endl;
    return 0;
}
