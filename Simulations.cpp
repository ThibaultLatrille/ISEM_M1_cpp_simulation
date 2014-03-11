// By Thibault Latrille
// thibault.latrille@ens-lyon.fr
// 11/03/2014
// No copyright (Why on hell should there be one anyway!)

// This simulation intend verify our analytical results by a Monte Carlo algorithm, refer to the paper for more information.
// This piece of code is specifically aimed at testing the following formula (for 2 families):
// \mathbb{E}\left[ \sum_{i=1}^2 \left. \dfrac{N_i(t)(N_i(t)-1)}{N_+(t)( N_+(t)-1 ) } \right\vert N_+(t) \right]
// = \dfrac{r_1(r_2+1)+r_2(r_2+1)}{r_+ (r_+ +1 )}-\dfrac{2}{ (r_+ +1 )}  \dfrac{1}{N_+(t)-1}


#include <iostream>		// setbuf, NULL
#include <cstdio>		// stdout
#include <cstdlib>
#include <ctime>	    // time as seed for pseudo-random number generator
#include <vector> 		// vector
#include <cmath> 		// sqrt, pow
#include <numeric>		// inner_product, accumulate

using namespace std;

int main() {

	setbuf(stdout, NULL); // Specifically for eclipse IDE, don't buffer output

	int r1 = 50;  		// Initial size of family 1
	int n1;				// Size of family 1
	int r2 = 100; 		// Initial size of family 2
	int n2;				// Size of family 2
	int rplus = r1 + r2;  // Initial size of population
	int nplus;			// Size of population
	int nfinal = 1000; 	// Size of population at the end of computation
	int random;
	srand(time(NULL)); // random seed
	cout << "r1=" << r1 << endl;
	cout << "r2=" << r2 << endl;

	int number_of_pop = 1000; // The number of populations we are computing
	vector<double> relatedness_replicate;
	// The vector that will contain the estimated relatedness for each independent population
	relatedness_replicate.reserve(number_of_pop); // Allocate memory for the vector

	for (int i = 0; i<number_of_pop; i++)  // For each replicated population
	{
		n1 = r1;
		n2 = r2;
		nplus = n1 + n2;

		while (nplus < nfinal)  // Monte Carlo algorithm for generating one population
		{
			random = rand() % nplus + 1;
			if (random <= n1)
				++n1;
			else
				++n2;
			++nplus;
		}
		relatedness_replicate.push_back(double(n1*(n1 - 1) + n2*(n2 - 1)) / double(nplus*(nplus - 1)));
		// Append the estimated relatedness to the vector
	}

	double s1 = accumulate(relatedness_replicate.begin(), relatedness_replicate.end(), 0.0);
	double mean_relatedness_replicate = s1 / number_of_pop; // The mean of estimated relatedness

	double s2 = 0.0; // The sum of square of estimated relatedness
	for (vector<double>::iterator it = relatedness_replicate.begin(); it != relatedness_replicate.end(); it++)
	{
		s2 += pow(*it, 2.0);
	}
	double stdev = 1.96*sqrt((s2 / number_of_pop - pow(mean_relatedness_replicate, 2.)) / (number_of_pop - 1));
	// The standard deviation of estimated relatedness

	double relatedness_upper_estimate = mean_relatedness_replicate + stdev; // Upper bound for the 95% confidence interval
	double relatedness_lower_estimate = mean_relatedness_replicate - stdev; // Lower bound for the 95% confidence interval

	double relatedness = double(r1*(r1 + 1) + r2*(r2 + 1)) / double(rplus*(rplus + 1)) - 2. / double((rplus + 1)*(nplus - 1));
	// The relatedness given by analytic formula

	cout << "n1=" << n1 << endl;
	cout << "n2=" << n2 << endl;
	cout << "ntotal=" << nplus << endl;
	cout << "relatedness=" << relatedness << endl;
	cout << "estimated relatedness=" << mean_relatedness_replicate << endl;
	cout << "with confidence interval [" << relatedness_lower_estimate << "," << relatedness_upper_estimate << "]" << endl; // display results

	return 0;
}
