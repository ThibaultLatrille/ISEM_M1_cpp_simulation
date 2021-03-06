// By Thibault Latrille
// thibault.latrille@ens-lyon.fr
// 11/03/2014
// No copyright (Why on hell should there be one anyway!)

// This simulation intend verify our analytical results by a Monte Carlo algorithm, refer to the paper (also in the github repository) for more information.
// This piece of code is specifically aimed at testing the following formula (for 2 families):
//  \mathbb{E}\left[ \left. \dfrac{N_i(t)(N_i(t)-1)}{N_+(t)( N_+(t)-1 ) } \right\vert N_+(t)=n_+ \right]
//  =\dfrac{r_i(r_i+1)}{r_+ (r_+ +1)}- \dfrac{2 r_i (r_+ -r_{i})}{r_+ (r_+ +1)} \dfrac{1}{ n_+ -1  }.
//
// Write severals file of the form r1=20_r2=500_replicate=10000.txt in the current directory.
// The first line is the header.
// Each lines contains the estimated relatedness, the 95% confidence interval as long as the true value for one fixed population size.
// Note that two estimated relatedness on different lines are not independent of one another.

#include <iostream>		// setbuf, NULL
#include <cstdio>		// stdout
#include <cstdlib>      // standard library
#include <ctime>	    // time as seed for pseudo-random number generator
#include <vector> 		// vector
#include <cmath> 		// sqrt, pow
#include <numeric>		// inner_product, accumulate
#include <fstream>		// ofstream
#include <sstream>		// ostringstream
#include <unistd.h>     // chdir
#include <random>
#include <chrono>
#include "Simulations.h"

int exponential_growth( void ){
	const double rho=0.2;
	double tmax=20.;
	const int r=100;

	unsigned seed = chrono::system_clock::now().time_since_epoch().count(); // The seed
	mt19937 generator(seed); // The mersenne_twister generator
	uniform_real_distribution<double> distribution(0.0,1.0); // Uniform [0,1] distribution

	const int number_of_pop = 1000;
	// The number of populations we are computing, influence the range of the confidence interval
	const int number_of_time_steps = 250;
	// The number of time steps we are interested in

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


int conditional_expectation( void ){
	int r_list[9]={10,20,30,40,50,100,200,500,1000};

	for (int r1_index=0; r1_index<9; r1_index++) // For every r1 in r_list
	{
		for (int r2_index=0; r2_index<9; r2_index++) // And the for every r2 in r_list
		{
			int r1 = r_list[r1_index];  		// Initial size of family 1
			int n1;				// Size of family 1
			int r2 = r_list[r2_index]; 		// Initial size of family 2
			int n2;				// Size of family 2
			int rplus = r1+r2;  // Initial size of population
			int nplus;			// Size of population
			int random;
			srand (time(NULL)); // random seed

			int number_of_steps=1000; // The number of estimation produced
			// The number of lines in the files
			double output_array[number_of_steps][5]; // The array containing all the relevant information
			for (int i=0; i<number_of_steps; i++)    // Initiate all the population size at which we will compute relatedness
			{
				output_array[i][0] = rplus+i*rplus/50; // the max population size is (rplus+ number_of_steps-1)*rplus/100
			}

			int number_of_pop = 1000;
			// The number of populations we are computing, influence the range of the confidence interval

			vector<double> relatedness_replicate( number_of_pop );
			// The vector that will contain the estimated relatedness for each independent population

			int n1_replicate[number_of_pop]; // Array to keep track of the size of family 1
			fill(n1_replicate, n1_replicate+number_of_pop, r1); // Initiated at r1
			int n2_replicate[number_of_pop]; // Array to keep track of the size of family 2
			fill(n2_replicate, n2_replicate+number_of_pop, r2); // Initiated at r1

			string file_name = "r1="+IntToStr(r1)+"_r2="+IntToStr(r2)+"_replicate="+IntToStr(number_of_pop)+".txt";
			// The file name
			cout << file_name << endl;

			ofstream myfile;
			myfile.open( file_name.c_str() ); // Open the file
			myfile << "//total_population estimated_relatedness lower_bound upper_bound relatedness" << endl;
			// The first line of the file, the header

			double s1;
			double stdev;
			double s2;

			for (int i=0; i<number_of_steps; i++)
			{
				for (int j=0; j<number_of_pop; j++)  // For each replicated population
				{
					n1 = n1_replicate[j]; // Initiate the size of the jth family (replicate) 1 to the previous state
					n2 = n2_replicate[j]; // Initiate the size of the jth family (replicate) 2 to the previous state
					nplus=n1+n2; // Replace

					while ( nplus < output_array[i][0] )  // Monte Carlo algorithm for generating one population
					{
						random = rand() % nplus + 1;// Pick a random number between 1 and nplus
						if (random <= n1) 			// If the random number is lower or equal to n1, add an individual to n1
							++n1;
						else						// Otherwise add an individual to n2
							++n2;
						++nplus;					// Add an individual to nplus
					}

					n1_replicate[j] = n1;
					n2_replicate[j] = n2;
					relatedness_replicate[j] = ( double (n1) * double (n1-1) + double(n2)*double(n2-1)  ) / ( double(nplus) *double (nplus-1) ) ;
					// Append the estimated relatedness to the vector
				}

				s1 = accumulate( relatedness_replicate.begin() , relatedness_replicate.end() , 0.0 );
				output_array[i][1] = s1 / number_of_pop; // The mean of estimated relatedness

				s2 = 0.0; // The sum of square of estimated relatedness
				for ( vector<double>::iterator it = relatedness_replicate.begin(); it != relatedness_replicate.end(); it++ )
				{
					s2+=pow(*it,2.0);
				}
				stdev = 1.96*sqrt( (s2/number_of_pop-pow( output_array[i][1],2.) ) / (number_of_pop-1) );
				// The standard deviation of estimated relatedness

				output_array[i][2] = output_array[i][1]-stdev; // Upper bound for the 95% confidence interval
				output_array[i][3] = output_array[i][1]+stdev; // Lower bound for the 95% confidence interval

				output_array[i][4] = double ( double (r1) *(r1+1)+double (r2)*(r2+1) ) / double ( rplus*double(rplus+1) ) - 2. * double ( r1*r2+r2*r1 ) /  double ( rplus*(rplus+1)*double(nplus-1) );
				// The relatedness given by analytic formula
				myfile << output_array[i][0] << " " << output_array[i][1] << " " << output_array[i][2] << " " << output_array[i][3] << " " << output_array[i][4] << endl;
				// Write relevant informations in the file
			}
			myfile.close(); // Write the file
		}
	}
	cout << "work done captain" << endl;
    return 0;
}

int simulate_infection( bool stoch_time, bool stoch_growth, int number_nematode ) {

	const int number_simulation=1000;
	const int r=100;

	double runi;

	////////////////////////////////////
    string file_name=IntToStr(number_nematode)+"nematodes_";

    if (stoch_time) file_name += "randtime_";
    else file_name += "dettime_";
    if (stoch_growth) file_name += "randgrowth";
    else file_name += "detgrowth";
    file_name += ".txt";

	cout << file_name << endl;

	ofstream myfile;
	myfile.open( file_name.c_str() ); // Open the file
    myfile << "//tau lambda mean lower_bound upper_bound" << endl;
			// The first line of the file, the header

	//unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	unsigned seed=123456;
	mt19937 generator(seed);
	uniform_real_distribution<double> u_distrib(0.0,1.0);

	const double lambda=0.01;
	for (int tau_rank=0; tau_rank<1000; tau_rank++){
        vector<double> vector_replicate;
        vector<int> number_infecting;

        double tau=pow(10,-2+double(tau_rank)*4/1000)*lambda;


        for (int ith_simu=0; ith_simu<number_simulation; ith_simu++){

            double interval_time=0;
            double sum_lineage=double(r);
            vector<double> family_size;
            vector<double> time_growth;
            time_growth.push_back(0.);

            for (int ith_nematode=0 ; ith_nematode<(number_nematode-1) ; ith_nematode++ ){
                if (stoch_time){
                    runi = u_distrib(generator);
                    interval_time = -log(runi)/(tau);
                }
                else interval_time = 1./(tau);

                if ((sum_lineage*exp(lambda*interval_time)+r)<pow(10,9)){
                    for ( vector<double>::iterator it = time_growth.begin(); it != time_growth.end(); it++ ){
                        *it+=interval_time;
                    }
                    time_growth.push_back(0.);
                    sum_lineage*=exp(lambda*interval_time);
                    sum_lineage+=r;
                }
                else break;
            }
            interval_time=log(pow(10,9)/sum_lineage)/lambda;
            for ( vector<double>::iterator it = time_growth.begin(); it != time_growth.end(); it++ ){
                *it+=interval_time;
            }
            int infection=0;
            for ( vector<double>::iterator it = time_growth.begin(); it != time_growth.end(); it++ ){
                infection++;
                if (stoch_growth){
                    double p=exp(-(*it)*lambda);
                    negative_binomial_distribution<unsigned long long> nb_distrib(r,p);
                    family_size.push_back(r+double(nb_distrib(generator)));
                }
                else family_size.push_back(r*exp((*it)*lambda));
            }
            double sum_r=0.;
            double sum_rsquared=0.;

            for ( vector<double>::iterator it = family_size.begin(); it != family_size.end(); it++ ){
                sum_r+=(*it);
                sum_rsquared+=pow(*it,2);
            }

            number_infecting.push_back(infection);
            vector_replicate.push_back(sum_rsquared/pow(sum_r,2));
        }

        double s1_infection;
        s1_infection = accumulate( number_infecting.begin() , number_infecting.end() , 0.0 );
        double mean_infection = s1_infection / number_simulation;

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

        //double p=0.025;
        //int rank=floor(p*number_simulation);
        //sort (vector_replicate.begin() , vector_replicate.end());

        //double lower_bound(vector_replicate[rank]);
        //double upper_bound(vector_replicate[number_simulation-rank]);

        //double approximat_relatedness;
        //double q=exp(lambda/tau);

        //approximat_relatedness=(q-1)*(pow(q,number_nematode)+1)/((q+1)*(pow(q,number_nematode)-1));

        if (tau_rank*100%1000==0){
        cout << tau_rank*100/1000 << "% done" << endl;
        }
        myfile << tau/lambda << " " << mean_infection << " " << mean << " " << lower_bound_stdev << " " << upper_bound_stdev << endl;
	}
	myfile.close(); // Write the file
    return 0;
}

string IntToStr( int integer)
// Convert integer to string
{
 ostringstream mystring;
 mystring << integer;
 return mystring.str();
}

string DoubleToStr( double doublenumber)
// Convert integer to string
{
 ostringstream mystring;
 mystring << doublenumber;
 return mystring.str();
}

