#ifndef __BRENTS_METHOD_HPP
#define __BRENTS_METHOD_HPP

#include <cmath>
#include "blackscholes.hpp"

double SIGN(const double a, const double b){
	if(b>=0.0){
		return fabs(a);
	}
	else{
		return -fabs(a);
	}
}

double brents_method(const double s, const double k, const double r, 
					const double t, const double call_option_price, 
					const double x1, const double x2, const double tol)
{

	int ITMAX = 100;
	double epsilon = 3.0e-8;

	bsc::BlackScholes bsc(s, k, r, t);

	double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2;
	double fa = call_option_price - bsc.OptionEurop(a);
	double fb = call_option_price - bsc.OptionEurop(b);

	double fc, p, q, r_, s_, tol1, xm;

	if((fa>0.0 && fb<0.0) || (fa<0.0 && fb<0.0))
		return -400.0;

	fc = fb;

	for(int i=0;i<ITMAX;i++){
		if((fb>0.0 && fc>0.0)||(fb<0.0 && fc<0.0)){
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if(fabs(fc)<fabs(fb)){
			a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
		}
		tol1 = 2.0*epsilon*fabs(b)+0.5*tol;
		xm=0.5*(c-b);

		if(fabs(xm)<=tol1||fb==0.0)
			return b;

		if(fabs(e)>=tol1 && fabs(fa) > fabs(fb)) {
			s_ = fb/fa;
			if(a==c){
				p=2.0*xm*s;
				q=1.0-s_;
			}
			else{
				q=fa/fc;
				r_=fb/fc;
				p=s_*(2.0*xm*q*(q-r_)-(b-a)*(r_-1.0));
				q=(q-1.0)*(r_-1.0)*(s_-1.0);
			}
			if (p > 0.0) 
				q = -q; // Check whether in bounds.
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {                     
                e=d;             
                b += d;
            } else {
                b += SIGN(tol1,xm);
            }
            fb = call_option_price - bsc.OptionEurop(b);
		}
	}
	return -200.00;
}

#endif