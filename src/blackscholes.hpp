#ifndef __BLACKSCHOLES_H
#define __BLACKSCHOLES_H

#include <iostream>
#include <cmath>
#include <armadillo>

class BlackScholes{

        private:
                double S_;
                double K_;
                double r_;
                double T_;

        public:
                BlackScholes(double S_,double K_,double r_,double T_);
                double OptionEurop(double);
                double OptionVega(double);

        
};


#endif 