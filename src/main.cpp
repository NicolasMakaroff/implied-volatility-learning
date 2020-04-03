#include "blackscholes.hpp"
//#include "interval_bisection.h"
#include "newton.hpp"
#include "brent.hpp"
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
        double low_vol = 0.01;
	double high_vol = 0.70;

	double epsilon = 0.001;
	double init = 0.3;
	  

	BlackScholes bs(S, K, r, T);
  	double sigma_ns = Newton(S,K,r,T,C_M, init, epsilon); 
	  	
	// << "Implied Vol using Newton-Raphson: " << sigma_ns << std::endl;
        
        BlackScholes brent(S,K,r,T);
        double sigma_brent = brents_method(S, K, r, T, C_M, low_vol, high_vol, epsilon);
        //cout << "Implied Vol using Brent: " << sigma_brent << std::endl;

        for(int i=0;i<100000;i++){
		S += pow(-1,i)*10;  
	  	K = S;  
	  	C_M += 0.1*pow(-1,i);
	  
	  	//Black-Scholes functor
	  	cout<<"Before BSC "<<i<<endl;
	  	BlackScholes bsc(S, K, r, T);
  		double sigma_bm = brents_method(S, K, r, T, C_M, low_vol, high_vol, epsilon);

  		low_vol = 0.90*sigma_bm;
  		high_vol = 1.10*sigma_bm;
	  	
		cout << "Implied Vol using Brent: " << sigma_bm << std::endl;
	}	

        return 0;
}