import numpy as np
import scipy.stats as scp

def Heston_discret1(kappa_,theta_,sigma_,rho_,r_,T_,L_,V0_,S0_):
    V = [V0_]
    X = [np.log(S0_)]
    delta_t_ = T_/L_
    delta_W_ = np.array([scp.norm.rvs(loc = 0, scale = np.sqrt(delta_t_), size = L_)]).T
    delta_W_orth_ = np.array([scp.norm.rvs(loc = 0, scale = np.sqrt(delta_t_), size = L_)]).T
    for i in range(L_):
        Vi_ = V[i]+kappa_*(theta_-V[i])*delta_t_ + sigma_*np.sqrt(max(V[i],0))*(rho_*delta_W_[i][0] + np.sqrt(1-rho_**2)*delta_W_orth_[i][0])
        V.append(Vi_)
        Xi_ = X[i] + (r_-1/2*V[i])*delta_t_ + np.sqrt(max(V[i],0)) * delta_W_[i][0]
        X.append(Xi_)
    return X[L_]

def f(X_,K_):
    return max(np.exp(X_)-K_,0)
    
def monte_carlo(kappa_,theta_,sigma_,rho_,r_,T_,L_,V0_,S0_,K0_,N_):
    X = []
    for i in range(N_):
        X.append(Heston_discret1(kappa_,theta_,sigma_,rho_,r_,T_,L_,V0_,S0_))
    s = 0
    for i in range(N_):
        s += f(X[i],K0_)
    return s/N_
    
def call(esp_,r_,T_):
    return np.exp(-r_*T_)*esp_