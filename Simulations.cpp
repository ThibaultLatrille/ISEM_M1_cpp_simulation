#include <iostream>		// setbuf, NULL
#include <cstdio>		// stdout
#include <cstdlib>		//
#include <ctime>	    // time
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
	vector<double> relatedness_replicate; // The vector for each computed pop
	relatedness_replicate.reserve(number_of_pop); // allocate memory for the vector

	for (int i = 0; i<number_of_pop; i++)
	{
		n1 = r1;
		n2 = r2;
		nplus = n1 + n2;

		while (nplus < nfinal)
		{
			random = rand() % nplus + 1;
			if (random <= n1)
				++n1;
			else
				++n2;
			++nplus;
		}
		relatedness_replicate.push_back(double(n1*(n1 - 1) + n2*(n2 - 1)) / double(nplus*(nplus - 1)));
	}

	double s1 = accumulate(relatedness_replicate.begin(), relatedness_replicate.end(), 0.0);
	double mean_relatedness_replicate = s1 / number_of_pop;

	double s2 = 0.0;
	for (vector<double>::iterator it = relatedness_replicate.begin(); it != relatedness_replicate.end(); it++)
	{
		s2 += pow(*it, 2.0);
	}
	double stdev = 1.96*sqrt((s2 / number_of_pop - pow(mean_relatedness_replicate, 2.)) / (number_of_pop - 1));

	double relatedness_upper_estimate = mean_relatedness_replicate + stdev;
	double relatedness_lower_estimate = mean_relatedness_replicate - stdev;

	double relatedness = double(r1*(r1 + 1) + r2*(r2 + 1)) / double(rplus*(rplus + 1)) - 2. / double((rplus + 1)*(nplus - 1));

	cout << "n1=" << n1 << endl;
	cout << "n2=" << n2 << endl;
	cout << "ntotal=" << nplus << endl;
	cout << "relatedness=" << relatedness << endl;
	cout << "estimated relatedness=" << mean_relatedness_replicate << endl;
	cout << "with confidence interval [" << relatedness_lower_estimate << "," << relatedness_upper_estimate << "]" << endl;

	return 0;
}
