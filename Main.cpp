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

    chdir("/home/thibault/Desktop/ISEM_M1_Cpp_Simulations/simulated_data/");
    for (int nematode=2; nematode<15; nematode++)
    {
    simulate_infection(true, true,nematode);
    simulate_infection(true, false,nematode);
    simulate_infection(false, true,nematode);
    simulate_infection(false, false,nematode);
    }
}
