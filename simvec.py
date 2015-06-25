# -*- coding: utf-8 -*-
# Setting ODES for simulating local neural activity

import numpy
from Sim_variables import *
from scipy.io import loadmat
from Support import gain

CIJ = loadmat('./macaque47.mat')
CM = CIJ['CIJ']
k_in = numpy.sum(CM*1.0,axis=0)    

def simvec(y,t):
#    print V1, V2, V3, V4, V5, V6, V7, gCa, gK, gL, VK, VL, VCa, I, b, ani, aei 
#    print aie, aee, phi, V8, V9, gNa, VNa, ane, nse, rnmda, N, CM, vs, c, k_in
    
    out = numpy.empty(shape=(N*3,1))
    v   = y[::3]
    w   = y[1::3]
    z   = y[2::3]
    
    gainv = gain(v,V5,V6,0.5);    
    
    vgainmean = (numpy.dot(gainv.transpose(),CM))/k_in

    vgainmean = vgainmean.transpose()
    
    vgainmean[numpy.isnan(vgainmean)] = 0

    out[::3] = -(gCa+rnmda*(1.0-c)*aee*gainv+rnmda*c*aee*vgainmean)*gain(v,V1,V2,0.5)*(v-VCa) \
        -gK*w*(v-VK) \
        -gL*(v-VL)  \
        -(gNa*gain(v,V9,V8,0.5)+(1.0-c)*aee*gainv+c*aee*vgainmean)*(v-VNa) \
        +ane*I \
        -gain(z,V7,V6,aie*0.5)*z
                
    out[1::3] = phi*(gain(v,V3,V4,0.5)-w);

    out[2::3] = b*(ani*I+gain(v,V5,V6,aei*0.5)*v);
    
    return out

def sintext(y,t):
    return numpy.sin(t*1.0)

