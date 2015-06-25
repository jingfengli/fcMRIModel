# -*- coding: utf-8 -*-
# Reorder results for visualization
#
import scipy 
import numpy
from Sim_variables import OrderIndex

def gain(VAR,C1,C2,C3):
    return C3*(1+scipy.tanh((VAR-C1)/C2))*1.0
    
def rowswap(OLDOrder):
    NEWOrder = numpy.empty(shape=numpy.shape(OLDOrder))
    for ai in xrange(47):
        for bi in xrange(47):
            NEWOrder[OrderIndex.index(ai + 1),OrderIndex.index(bi + 1)] = OLDOrder[ai,bi]
#    NEWOrder[OrderIndex-1==ai,OrderIndex -1 ==bi,:] = OLDOrder[ai,bi,:]
    return NEWOrder

