#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>

using namespace std;

int main() {

	setbuf(stdout, NULL);

	int r1 = 50;
	int n1 = r1;
	int r2 = 100;
	int n2 = r2;
	int rplus = r1 + r2;
	int nplus;
	int nfinal=1000;
	int random;

    cout << "r1=" << r1 << endl;
    cout << "r2=" << r2 << endl;

	srand (time(NULL));

	nplus = n1 + n2;

	while ( nplus < nfinal )
	{
		random = rand() % nplus + 1;
		if (random <= n1)
			++n1;
		else
			++n2;
		++nplus;
	}

	double relatedness_estimated;
	double relatedness;

	relatedness = double ( r1 * (r1 + 1 ) ) / double ( rplus * ( rplus + 1) ) - 2. * (rplus - r1 ) /  double ( rplus * ( rplus + 1) * ( nplus -1) );
	relatedness_estimated = double ( n1 * ( n1 - 1 ) ) / double ( nplus * ( nplus -1 ) );

    cout << "n1=" << n1 << endl;
    cout << "n2=" << n2 << endl;
    cout << "ntotal=" << nplus << endl;
    cout << "relatedness=" << relatedness << endl;
    cout << "estimated relatedness=" << relatedness_estimated << endl;

    return 0;
}
