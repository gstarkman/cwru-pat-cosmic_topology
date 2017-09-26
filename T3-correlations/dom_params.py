# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 14:06:48 2016

@author: jiggly
"""

import numpy as np

#fundamental domain lengths
lx_min, lx_max, lx_grid = 10, 25, 15
ly_min, ly_max, ly_grid = 10, 25, 15
lz_min, lz_max, lz_grid = 10, 25, 15
lengArr = [[] for _ in range(2)]
#for lx in range(lx_min, lx_max+lx_grid, lx_grid):
#    for ly in range(ly_min, ly_max+ly_grid, ly_grid):
#        for lz in range(lz_min, lz_max+lz_grid, lz_grid):
#            domLeng = [lx, ly, lz]
#            lengArr += [domLeng]
lengArr[0] += [10, 10, 10]
lengArr[1] += [25, 25, 25]
leng_fname = 'lengthData.csv'
np.savetxt(leng_fname, lengArr)

#fundamental domain tilt angles
theta_min, theta_max, theta_grid = 0, np.pi/2, 0.1*np.pi
phi_min, phi_max, phi_grid = 0, np.pi/2, 0.1*np.pi
psi_min, psi_max, psi_grid = 0, np.pi/2, 0.1*np.pi
angArr = [[]]
#for theta in np.arange(theta_min, theta_max+theta_grid, theta_grid):
#    for phi in np.arange(phi_min, phi_max+phi_grid, phi_grid):
#        for psi in np.arange(psi_min, psi_max+psi_grid, psi_grid):
#            domAngs = [theta, phi, psi]
#            angArr += [domAngs]
angArr[0] += [np.pi/2, np.pi/2, np.pi/2]
ang_fname = 'angleData.csv'
np.savetxt(ang_fname, angArr)

#orientation angles of topology w.r.t. observer
alpha_min, alpha_max, alpha_grid = 0, np.pi/2, 0.1*np.pi
beta_min, beta_max, beta_grid = 0, np.pi/2, 0.1*np.pi
gamma_min, gamma_max, gamma_grid = 0, np.pi/2, 0.1*np.pi
eulerArr = [[]]
#for alpha in np.arange(alpha_min, alpha_max+alpha_grid, alpha_grid):
#    for beta in np.arange(beta_min, beta_max+beta_grid, beta_grid):
#        for gamma in np.arange(gamma_min, gamma_max+gamma_grid, gamma_grid):
#            eulerAngs = [alpha, beta, gamma]
#            eulerArr += [eulerAngs]
eulerArr[0] += [np.pi/2, np.pi/2, np.pi/2]          
eul_fname = 'eulerAngData.csv'
np.savetxt(eul_fname, eulerArr)

#position w/in topology
rTx_min, rTx_max, rTx_grid = lx_min, lx_max, 1
rTy_min, rTy_max, rTy_grid = ly_min, ly_max, 1
rTz_min, rTz_max, rTz_grid = lz_min, lz_max, 1
rTarr = [[]]
#for rTx in range(rTx_min, rTx_max, rTx_grid):
#    for rTy in range(rTy_min, rTy_max, rTy_grid):
#        for rTz in range(rTz_min, rTz_max, rTz_grid):
#            rT = [rTx, rTy, rTz]
#            rTarr += [rT]
rTarr[0] += [1, 1, 1]
pos_fname = 'positionData.csv'
np.savetxt(pos_fname, rTarr)
