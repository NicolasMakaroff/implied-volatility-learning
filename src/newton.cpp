#include "newton.hpp"
#include "blackscholes.hpp"

#include<cmath>

double Newton(double S_, double K_, double r_, double T_, double C_M, double init, double epsilon){

        bsc::BlackScholes bs(S_, K_, r_, T_);



	double y = bs.OptionEurop(init);
	double x = init;

	while(fabs(y-C_M) > epsilon) {
		double d_x = (bs.OptionVega)(x);
		x+=(C_M - y)/d_x;
		y = (bs.OptionEurop)(x);
	}
	return x;

}