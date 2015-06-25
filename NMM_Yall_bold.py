# -*- coding: utf-8 -*-
import numpy
import scipy
from scipy.io import savemat
from Sim_variables import rho,V0
from Ode23 import rk4
from boldodes import boldodes

def NMM_Yall_bold(Yall):
    tfrom = 0;
    tto = numpy.size(Yall,0)

    N = numpy.size(Yall,1)
    
    # Solve BOLD ODEs
    # x(0) = s = vasodilatory signal
    # x(1) = f = inflow
    # x(2) = v = blood volume
    # x(3) = q = deoxyhaemoglobin content
    
    # initialize bold signal
    Ybold = numpy.empty(shape=(tto,N))
    
    # loop over time series
    for n in xrange(N):

        # get time series
        z = Yall[tfrom:tto,n]
        
        # get first 30 seconds and append in front in order to get rid of initial transiet
        ttrns = 30000
        zt = Yall[tfrom:ttrns,n]
        
        # get abs(diff(glu)) 
        # Maybe try reverse padding

        z = numpy.abs(numpy.diff(numpy.concatenate((zt,z),axis=0),axis=0))

        # the shape of z is weird now.
        z = numpy.concatenate((z, numpy.zeros(shape=(1,))),axis = 0)
        
        # ICs
        ics = numpy.array([0,1,1,1])
        
        
        # solve
        tzero = 0
        tend = (tto + ttrns - 1) /1000.0
        
        # arange function does not include the end, thus, we added a 0.001 to theoretic end
        x = rk4(boldodes, ics,numpy.arange(tzero,tend+0.001,0.001), z = z); 
#        print numpy.arange(tzero,0.001,tend)
        # BOLD signal
        k1 = 7.0 * rho
        k2 = 2
        k3 = 2.0 * rho - 2
        
        y = V0*(k1*(1-x[:,3])+k2*(1-x[:,2]/x[:,2])+k3*(1-x[:,2]))
        
        # remove transient
#        print numpy.shape(x)
        Ybold[:,n] = y[ttrns:]
        
    return Ybold

