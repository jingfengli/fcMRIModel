# -*- coding: utf-8 -*-
#
# Simulate BOLD signal based on local neural activity
#
################################################################# 	
#################################################################
# balloon-windkessel model in Friston KJ, Harrison L, & Penny W
# (2003) Dynamic causal modelling. NeuroImage 19(4):1273-1302
#################################################################

import numpy 
import scipy
from scipy.io import loadmat
from Sim_variables import kappa,gamma,tau,alpha,rho

def boldodes(x,t,z):
    # Haemodynamic model embedding the Balloon-Windkessel model
    # x(0) = s = vasodilatory signal
    # x(1) = f = inflow
    # x(2) = v = blood volume
    # x(3) = q = deoxyhaemoglobin content
    # z = neuronal activity
    
    tind = numpy.round(t*1000)
        
    zt = z[tind]

    
    y = numpy.empty(shape=(4,))
    
    y[0] = zt - kappa * x[0] - gamma * (x[1] - 1)
    y[1] = x[0]
    y[2] = (x[1] - (x[2]) ** (1.0/alpha)) / tau
    y[3] = (x[1] * (1-(1-rho)**(1/x[1]))/rho - (x[2])**((1-alpha)/alpha)*x[3])/tau
    return y


