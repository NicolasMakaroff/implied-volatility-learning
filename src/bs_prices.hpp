#ifndef __BS_PRICES_H
#define __BS_PRICES_H

#include <iostream>
#include <cmath>
#include <armadillo>

double d_j(const int j, const double S, const double K, const double r, const double sigma, const double T) {
  return (log(S/K) + (r + (pow(-1,j-1))*0.5*sigma*sigma)*T)/(sigma*(pow(T,0.5)));
}

double call_price(const double S, const double K, const double r, const double sigma, const double T) {
  return S * arma::normcdf(d_j(1, S, K, r, sigma, T))-K*exp(-r*T) * arma::normcdf(d_j(2, S, K, r, sigma, T));
}

double call_vega(const double S, const double K, const double r, const double sigma, const double T){
	return S*sqrt(T)*arma::normpdf(d_j(1, S, K, r, sigma, T));
}





#endif 