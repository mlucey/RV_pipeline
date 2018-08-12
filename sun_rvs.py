from astropy.io import fits
import numpy as np
#import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join
from pyraf import iraf
iraf.rv()

#get skyflat rvs

#**need to edit obsrun and the paths**

obsrun = ['Apr2018','Mar2018','Dec2016','Feb2017','Feb2016','Mar2016']

for t in obsrun:
    outpath = '/Volumes/MADDIE/'+t+'/results/donefxcor/'
    #where you have the skyflat image
    imagepath = '/Volumes/MADDIE/'+t+'/images/'
    #template your using
    temp = '/Volumes/MADDIE/newtemplates/new_sunspec.fits'

    fitsfile = fits.open(imagepath+'RED.ms.fits')
    head = fitsfile[0].header
    numstars = head['NAXIS2']

    aps = []
    for i in range(numstars):
        splits = head['APNUM'+str(i+1)].split(' ')
        ap = splits[0]
        aps.append(ap)
    for j in range(len(aps)):
        iraf.rv.fxcor(imagepath+'RED.ms.fits' ,temp, apertures=aps[j] ,output= outpath+'Reds'+str(aps[j]) , function='gaussian', background='INDEF', verbose = 'txtonly', interactive = 'no')

