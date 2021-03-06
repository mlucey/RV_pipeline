from astropy.io import fits
import numpy as np
#import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join
from pyraf import iraf
iraf.rv()

#get skyflat rvs

obsrun = 'Mar2018'
outpath = '/Users/Natalie/mann/fxcor_rv/'+obsrun+'/results/'
#where you have the skyflat image
imagepath = '/Users/Natalie/mann/fxcor_rv/'+obsrun+'/'
#template your using
temp = '/Users/Natalie/mann/templates/new_sunspec.fits'

fits = fits.open(imagepath+'RED.ms.fits')
head = fits[0].header
numstars = head['NAXIS2']

aps = []
for i in range(numstars):
    splits = head['APNUM'+str(i+1)].split(' ')
    ap = splits[0]
    aps.append(ap)
for j in range(len(aps)):
    iraf.rv.fxcor(imagepath+'RED.ms.fits' ,temp, apertures=aps[j] ,output= outpath+'Reds'+str(aps[j]) , function='gaussian', background='INDEF', verbose = 'txtonly', interactive = 'no')

