# -*- coding: utf-8 -*-
# Setting up environmental variables
# Setting up the connectivity strength 'c'
#
import scipy.io as spy

#CIJ = spy.loadmat('/Users/jingfengli/Work/fcmodel/macaque47.mat')
#nCIJ = CIJ('CIJ')

V1 = -0.01; V2 = 0.15; V3 = 0.0; V4 = 0.3; V5 = 0.0; V6 = 0.65; V7 = 0.0; V9 = 0.3; V8 = 0.15;
gCa = 1.0; gK = 2.0; gL = 0.5; gNa = 6.7;
VK = -0.7; VL = -0.5; I = 0.3; b = 0.1; phi = 0.7; VNa = 0.53; VCa = 1.0;
ani = 0.4; vs = 1.0; aei = 2.0; aie = 2.0; aee = 0.36; ane = 1.0; rnmda = 0.25;

nse = 0.0;
N = 47.0;

cij = 0.15

c = cij

# Balloon model parameters
kappa = 0.65;   # rate of signal decay, s^-1                0.65
gamma = 0.41;   # rate of flow-dependent elimination, s^-1  0.41
tau = 0.98;     # haemodynamic transit time, s              0.98
alpha = 0.32;   # Grubb's exponent                          0.32
rho = 0.34;     # resting oxygen extraction fraction        0.34
V0 = 0.02;      # resting blood volume fraction             0.02

# row swap index

OrderIndex =[36,39,38,37,15,31,16,18,21,17,46,44,30,47,8,23,41,45,43,40,42,28,27,29,32,\
    26,20,19,13,5,25,24,12,14,6,33,2,35,34,22,11,1,10,9,7,4,3];
