#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 10:43:45 2017

@author: jiggly
"""

import numpy as np
import astropy.io.fits as fits

mask_file_1024 = 'COM_Mask_CMB-confidence-Tmask-IQU-commander_1024_R2.02_full.fits'
#mask_file_2048 = 'COM_Mask_CMB-confidence-Tmask-IQU-commander-field-Int_2048_R2.01_full.fits'

f = fits.open(mask_file_1024)
masks = f[1].data

masks_arr = np.zeros((np.size(masks),np.size(masks[0])))
for i in range(0,np.size(masks)):
    for j in range(0,np.size(masks[0])):
        masks_arr[i][j] = masks[i][0][j]

np.savetx('dataMask.csv', masks_arr)

