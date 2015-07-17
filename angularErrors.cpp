/*
 *   Use a Monte-Carlo method to determine the angular errors
 *
 *   Compile with g++ -o angularErrors angularErrors.cpp -lgsl
 *
 *   Morag I. Scrimgeour, 16/07/2015
 */

#define pi 3.1415926535897932384626433832795

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <gsl/gsl_rng.h>
#include </Users/Morag/Cosmology_codes/CPP_functions/usefulfunctions.h>
#include </Users/Morag/Cosmology_codes/CPP_functions/gslrandfunctions.h>

using namespace std;

int main() {
	
	double* u = new double[3];
	double* du = new double[3];
	double RA_l, Dec_b, mag, dmag, dRA_l, dDec_b;
	
	//u[0] = 0.;    // Equatorial coords
	//u[1] = 0.;
	//u[2] = -100.;
	
	//u[0] = -208.;    // Equatorial coords, RI = 50
	//u[1] = -99.;
	//u[2] = -91.;
	
	//u[0] = -203.;    // Equatorial coords, RI = 70
	//u[1] = -97.;
	//u[2] = -90.;
	
	u[0] = -212.;    // Equatorial coords, MLE
	u[1] = -125.;
	u[2] = 162.;
	
	
	
	du[0] = 3.;
	du[1] = 8.;
	du[2] = 128.;
	
	mag = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
	
	//-------------------------------------------------------------------------------------------
	//------------------------------ Set up random number generator -----------------------------
	//-------------------------------------------------------------------------------------------
	
	gsl_rng * rr;
	long seed;
	
	gsl_rng_env_setup();
	
	rr = gsl_rng_alloc(gsl_rng_rand48);     // pick random number generator
	seed = time(NULL); // * getpid();
	gsl_rng_set(rr, seed);
	
	//-------------------------------------------------------------------------------------------
	//------------------------------ Convert BF -----------------------------
	//-------------------------------------------------------------------------------------------
	
	angular_errors_MC(u, du, RA_l, Dec_b, dmag, dRA_l, dDec_b, rr, 1, 1);
	// iopt: input option. 1: Equatorial, 2: Galactic Cartesian coords.
	// oopt: output option. 1: Equatorial (RA,Dec), 2: Galactic (l,b)
	
	cout << "mag = " << mag << " +/- " << dmag << endl;
	cout << "RA = " << RA_l << " +/- " << dRA_l << endl;
	cout << "Dec = " << Dec_b << " +/- " << dDec_b << endl;
	
	angular_errors_MC(u, du, RA_l, Dec_b, dmag, dRA_l, dDec_b, rr, 1, 2);
	cout << "l = " << RA_l << " +/- " << dRA_l << endl;
	cout << "b = " << Dec_b << " +/- " << dDec_b << endl;
	

}