
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import camb
#from camb import model, initialpower
#import csv
import numpy as np

l_max = 28
Boost = 3.0

pars = camb.model.CAMBparams()
powpars = camb.initialpower.InitialPowerParams
transpars = camb.TransferParams

#=============================================================================================================
#set parameter values for CAMB
#   document stating which to vary for Planck: https://wiki.cosmos.esa.int/planckpla/index.php/Simulation_data
#   values for model taken from: http://xxx.lanl.gov/pdf/1502.01589v2
#       note: delta_nnu = DeltaN_eff = N_eff(=3.13 +\- 0.64\0.63) - 3.046 = 0.084 (ignoring +\-)
#=============================================================================================================

pars.set_accuracy(AccuracyBoost = Boost, 
                  lSampleBoost = Boost,
                  lAccuracyBoost = Boost,
                  HighAccuracyDefault = True, 
                  DoLateRadTruncation = True)
                  
pars.set_bbn_helium(ombh2 = 0.02222, 
                    delta_nnu = 0.084, 
                    tau_neutron = 880.3)

pars.set_cosmology(H0 = 67.31,
                   cosmomc_theta = None,#1.04085e-2, 
                   ombh2 = 0.02222, 
                   omch2 = 0.1197, 
                   omk = -0.040, 
                   num_massive_neutrinos = 1, 
                   mnu = 0.7, 
                   nnu = 3.13, 
                   YHe = None, 
                   meffsterile = 0, 
                   standard_neutrino_neff = 3.13, 
                   TCMB = 2.7255, 
                   tau = 0.078, 
                   tau_neutron = 880.3)
                   
pars.set_dark_energy(w = -1.54, 
                     sound_speed = 1.0, 
                     dark_energy_model = 'fluid')
                     
pars.set_for_lmax(lmax = l_max, max_eta_k = None, 
                  lens_potential_accuracy = 0, 
                  lens_margin = 150, 
                  k_eta_fac = 2.5, 
                  lens_k_eta_reference = 18000.0)
                  
pars.InitPower.set_params(As = 2.198e-09, 
                          ns = 0.9655, 
                          nrun = 1, 
                          nrunrun = 0, 
                          r = 0, 
                          nt = None, 
                          ntrun = 0, 
                          pivot_scalar = 0.05, 
                          pivot_tensor = 0.05, 
                          parameterization = 2)
                          
pars.set_matter_power(redshifts = [0.0], 
                      kmax = 1.2, 
                      k_per_logint = None, 
                      silent = False)
                      
#set parameters for transfer functions
transpars(AccuracyBoost = Boost, high_precision = True)

#get transfer function from CAMB
data = camb.get_transfer_functions(pars)
transfer = data.get_cmb_transfer_data()
transfer_k = np.reshape(transfer.q, (len(transfer.q), 1))
transfer_data = []
transfer_data += [transfer_k[0:np.size(transfer_k)]]
for ell in range(0, 28):
    transfer_deltaplk = np.reshape(transfer.delta_p_l_k[0][ell], (len(transfer.delta_p_l_k[0][ell]), 1))
    transfer_data += [transfer_deltaplk[0:np.size(transfer_k)]]
#check validity of parameter choices
if pars.validate() == False:
    print('invalid parameters')
    quit()
transData = np.resize(transfer_data, (29, np.size(transfer_k)))
np.savetxt('transFuncData.csv', transData)
