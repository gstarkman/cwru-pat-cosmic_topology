#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 07:36:34 2017

@author: jiggly
"""

import healpy as hp
import healpy.visufunc as hpv
import healpy.pixelfunc as hpp
import matplotlib.pyplot as plt
import numpy as np

RCpp = np.loadtxt("/home/jiggly/Documents/Research/CorrelationMatrix/BayesianModComp/R_Cpp.csv")
ICpp = np.loadtxt("/home/jiggly/Documents/Research/CorrelationMatrix/BayesianModComp/I_Cpp.csv")
Csize = np.size(RCpp)

RCvv = np.loadtxt("/home/jiggly/Documents/Research/CorrelationMatrix/BayesianModComp/R_Cvv.csv")
ICvv = np.loadtxt("/home/jiggly/Documents/Research/CorrelationMatrix/BayesianModComp/I_Cvv.csv")
corrsz = np.size(RCvv)

res = 3 #res from cpp code
nside = 2**res
pnum = hpp.nside2npix(nside)
xCpp = np.zeros((pnum,pnum))
yCpp = np.zeros((pnum,pnum))
zCpp = np.zeros((pnum,pnum))

modeno = int(np.sqrt(corrsz))
xCvv = np.zeros((modeno,modeno))
yCvv = np.zeros((modeno,modeno))
zCvv = np.zeros((modeno,modeno))

for i in range(1, modeno):
    xCvv[i] = RCvv[(i-1)*modeno:i*modeno] #/ 10e40

plt.pcolor(xCvv)
plt.colorbar()

# loop over RCpp and store in each above
for i in range(1, Csize/pnum):
#    if (RCpp[(i-1)*pnum:i*pnum] < 10e50 and RCpp[(i-1)*pnum:i*pnum] < 10e-50):
    xCpp[i] = RCpp[(i-1)*pnum:i*pnum] / 10e110
    yCpp[i] = ICpp[(i-1)*pnum:i*pnum] / 10e110

#Cpad = np.lib.pad(ICpp, (pnum-Csize,0), 'constant', constant_values=(0,0))

hpv.mollview(np.abs(xCpp[3]), rot=90, norm='log', cbar=False, hold=False)
hpv.mollview(np.abs(xCpp[50]), rot=90, norm='log', cbar=False, hold=False)
hpv.mollview(np.abs(xCpp[200]), rot=90, norm='log', cbar=False, hold=False)
hpv.mollview(np.abs(xCpp[600]), rot=90, norm='log', cbar=False, hold=False)
#hpv.mollview(yCpp[3], rot=90, norm='hist')
