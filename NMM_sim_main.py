# -*- coding: utf-8 -*-
import time
import os
import numpy as np
from numpy.random import rand
from scipy.io import loadmat, savemat
import matplotlib.pyplot as plt
from Sim_variables import c
from simvec import simvec
from Ode23 import rk2
from NMM_Yall_bold import NMM_Yall_bold
from Support import rowswap

##############
# runname
init = 'randm'
#init = 'saved'
rn = 'PY_01'

##############
# CONNECTION MATRIX
CIJ = loadmat('./macaque47.mat')
CM = CIJ['CIJ']

N = np.size(CM,0)

###############
# TIME PARAMS ###
tseg = 1;       # number of segments used in the intial transient
lseg = 8;       # number of segments used in the actual run
sseg = 4;       # the positon of the segment for which fast neural activity is saved 
llen = 6000;    # length of each segment, in milliseconds
tres = 0.2;     # time resolution of model output, in milliseconds

##############
# Error Tolerances
# NOT SETUP YET

##############
# INITIAL CONDITION
# Initial condition - random

#######
if init == 'randm':
    ics = np.empty(shape=(N*3,1));
    for i in xrange(N):
        ics[i*3:i*3+3] = np.array([(rand(1)-0.5)*0.8-0.2 , \
                                   (rand(1)-0.5)*0.6+0.3, \
                                   (rand(1)-0.5)*0.16+0.05])
elif init == 'saved':
    ICS = loadmat('/Users/jingfengli/Work/fcmodel/ics.mat')
    ics = ICS['ics']
    
##############
# START SIMULATION
# TRANSIENT
print '\n\nconnection strength is set to %s \n' %c
print 'beginning dynamics (transient) ...'

for SEGMENT in xrange(tseg):
    tic = time.time()
    # arange function does not include the end, thus, we added a tres to theoretic end
    y = rk2(simvec,ics,np.arange(0,llen+tres,tres));
    yics = y[-1,:]
    print 'finished segment %s' % SEGMENT
    ics = yics
    toc = time.time()
    print '%.2f' %(toc-tic)

print 'finished transient'
#END TRANSIENT

# save model parameters and initial condition
#savemat('./init_' + rn + '.mat',mdict={'ics':ics,'CM':CM, 'c':c})

#################
# RUN
# loop over 'lseg' segments of equal length
for SEGMENT in xrange(lseg):
    tic = time.time()
    # arange function does not include the end, thus, we added a tres to theoretic end
    y = rk2(simvec,ics,np.arange(0,llen+tres,tres));
    # keep only excitatory variable and downsample to 1 ms resolution
    Y = y[::5, ::3]
    # save last time step as initial condition for next time segment
    yics = y[-1]
    # save downsampled time series of excitatory variable, plus parameters
    savemat('./' + rn + '_part' + str(SEGMENT) + '.mat',mdict={'yics':yics,'Y':Y, 'c':c, 'tres':tres})
    print 'saved segment %d' % SEGMENT
    # swap initial condition
    ics = yics
    toc = time.time()
    print '%.2f' %(toc-tic)
#    print len(Y)
#################
# END OF RUN

# SAVE ONE RUN SEGMENT (fast activity)
savemat('./' + rn + '_Yfast.mat',mdict={'Y':Y, 'c':c, 'sseg': sseg})

# CONCATENATE OUTPUT FILES
Yall = np.empty((0,47,1))
# concatenate
for s in xrange(lseg):
    Yfile = loadmat('./' + rn + '_part' + str(s) + '.mat')
    Yall = np.concatenate((Yall,Yfile['Y'][:-1]),axis=0)

Yall = np.squeeze(Yall)
# Delete the small segements
#for s in xrange(lseg):
#    os.remove('./' + rn + '_part' + str(s) + '.mat')

####################
# COMPUTE BOLD SIGNAL
# using the nonlinear ballon-windkessel model
print "beginning bold calculation ..."
tic = time.time()
Ybold = NMM_Yall_bold(Yall)
toc = time.time()
print '%.2f' %(toc-tic)

# COMPUTE SOME BASIC BOLD SIGNAL ANALYSES ==========
# settings for bold averages
T = np.size(Ybold,axis = 0);
xsec = 2000 
xgap = 500       # in msec, window size and spacing
t0 = np.arange(0,T-xsec+xgap,xgap)
#len(t0)
te = np.arange(xsec,T+xgap,xgap)
#len(te)
# initialize...
Ybold_w = np.empty(shape=(len(t0),N));
# compute bold averages
for w in xrange(len(t0)):
    Ybold_w[w] = np.mean(Ybold[t0[w]:te[w]],axis=0)

# remove NaNs, get average bold signal over whole brain, and regress out

Ybold_w[np.isnan(Ybold_w)] = 0
Ybold_w_mean = np.mean(Ybold_w,axis = 1)

Ybold_w_reg = np.empty(shape=(len(t0),N))
for i in xrange(N):
    (lina,linb) = np.polyfit(Ybold_w_mean,Ybold_w[:,i],1)
    Ybold_w_reg[:,i] = Ybold_w[:,i] - lina*Ybold_w_mean - linb

# get BOLD cross-correlations
C = np.corrcoef(Ybold_w_reg.transpose());

# save processed BOLD data
#savemat('./' + rn + '_Yfast.mat',mdict={'Y':Y, 'c':c, 'sseg': sseg})

savemat('./' + rn + '_Ybold_proc.mat',mdict={'Y':Y, 'CIJ':CIJ, 'Yall': Yall , 'Ybold_w':Ybold_w, \
        'Ybold_w_mean':Ybold_w_mean,'Ybold_w_reg':Ybold_w_reg,'C':C})

print '... all done ...'

NC = rowswap(C)






