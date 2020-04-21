import numpy as np
import scipy.stats as scp

def d(S_, K_, r_,sigma_ ,tau_):
    d1_ = 1 / (sigma_ * np.sqrt(tau_)) * ( np.log(S_/K_) + (r_ + sigma_**2/2) * tau_)
    d2_ = d1_ - sigma_ * np.sqrt(tau_)
    return d1_, d2_

def call_option(S_, K_, r_,sigma_ ,tau_ ,d1_,d2_):
    return scp.norm.cdf(d1_) * S_ - scp.norm.cdf(d2_) * K_ * np.exp(-r_ *tau_)
    
def call_vega(S_, tau_ , d1_):
    return S_ * scp.norm.pdf(d1_) * np.sqrt(tau_)
    
def Newton_Raphson(S_,K_,r_,tau_, sigma0_ ,price_, epsilon_):
    d1_,d2_ = d(S_, K_, r_,sigma0_ ,tau_)
    g = call_option(S_, K_, r_,sigma0_ ,tau_,d1_,d2_)-price_
    sigma_ = sigma0_ - g/call_vega(S_,tau_, d1_)
    while np.abs(sigma_-sigma0_)/sigma0_ > epsilon_:
        sigma0_ = sigma_
        d1_,d2_ = d(S_, K_, r_,sigma0_ ,tau_)
        g = call_option(S_, K_, r_,sigma0_ ,tau_,d1_,d2_)-price_
        sigma_ = sigma0_ - g/call_vega(S_,tau_, d1_)
    return sigma_