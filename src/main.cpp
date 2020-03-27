#include "blackscholes.hpp"
//#include "interval_bisection.h"
#include "newton.hpp"
//#include "brents_method.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace bsc;

int main() {

        double S = 100.0;
  	double K = 100.0;
  	double C_M = 10.5;
  	double r = 0.06;   
	double T = 1.0;  

	double epsilon = 0.001;
	double init = 0.3;
	  

	BlackScholes bs(S, K, r, T);
  	double sigma_ns = Newton(S,K,r,T,C_M, init, epsilon); 
	  	
	cout << "Implied Vol using Newton-Raphson: " << sigma_ns << std::endl;

        return 0;
}