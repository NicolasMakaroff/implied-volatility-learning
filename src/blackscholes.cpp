#include "blackscholes.hpp"
#include "bs_prices.hpp"
#include <cmath>

bsc::BlackScholes::BlackScholes(double S_, double K_, double r_, double T_)
:S_(S_), K_(K_), r_(r_), T_(T_){

}

double bsc::BlackScholes::OptionEurop(double sigma_){
        return call_price(S_, K_, r_, sigma_, T_);
}

double bsc::BlackScholes::OptionVega(double sigma_){
        return call_vega(S_, K_, r_, sigma_, T_);
}