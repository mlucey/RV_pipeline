from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join
from pyraf import iraf

#get skyflat rvs

obsrun = 'Feb2016'
outpath = '/Volumes/MADDIE/'+obsrun+'/final_stuff/donefxcor/'
#where you have the skyflat image
imagepath = '/Volumes/MADDIE/'+obsrun+'/images/'
#template your using
temp = '/Volumes/MADDIE/RV_std/new_solar.fits'

fits = fits.open(imagepath+'RED.ms.fits')
head = fits[0].header
numstars = head['NAXIS2']

aps = []
for i in range(numstars):
    splits = head['APNUM'+str(i+1)].split(' ')
    ap = splits[0]
    aps.append(ap)
for j in range(len(aps)):
    iraf.rv.fxcor(imagepath+'RED.ms.fits' ,temp, apertures=aps[j] ,output= outpath+'Reds'+str(aps[j]) , verbose = 'txtonly', interactive = 'no')
