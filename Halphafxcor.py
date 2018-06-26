from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join
import itertools
import pandas
from scipy import optimize as op
import scipy.integrate as integrate
import math

from pyraf import iraf
iraf.rv()

#needs to be run in pyraf
#for getting halpha fxcor rvs
#assumes you already have an halpha template made

obsrun = 'Apr2018'
# where you have the object spectra
inpath = '/Volumes/MADDIE/'+obsrun+'/images/'
#the halpha tempfile
temp = '/Volumes/MADDIE/Forfxcor/1647952iplushalpha.fits'
#where you want the fxcor outputs saved
outpath = '/Volumes/MADDIE/'+obsrun+'/final_stuff/donefxcor/
#where you want the list of halpha stars saved
filepath = '/Volumes/MADDIE/'+obsrun+'/final_stuff/Halphalegend.txt'
#where you have non-halpha rvcorrect.txt
rvcorrect = '/Volumes/MADDIE/'+obsrun+'/final_stuff/rvcorrect.txt'

halphafiles = []
halphaaps = []
halpharas = []


hght = np.loadtxt(rvcorrect,usecols = (12,),dtype = str)
files = np.loadtxt(rvcorrect,usecols = (13,),dtype = str)
aps = np.loadtxt(rvcorrect,usecols = (10,))
ra = np.loadtxt(rvcorrect,usecols = (9,),dtype = str)

#grabs spectra that failed to get an rv when matched with a template without halpha emission
for i in range(len(hght)):
    if hght[i] == 'INDEF' or float(hght[i]) < .5:
        halphafiles.append(files[i])
        halphaaps.append(aps[i])
        halpharas.append(ra[i])
        
        
#runs fxcor and writes a file with the ras of each star used 
file = open(filepath,'w')
for k in range(len(halphafiles)):
    iraf.rv.fxcor(inpath+halphafiles[k]+'_nosky.fits',temp,apertures = str(int(halphaaps[k])),output = outpath+"halpha_"+halphafiles[k]+"_"+str(halpharas[k]),verbose='txtonly')
    file.write("%8s %8s %8.3f \n" % (halpharas[k],'1647952i',-71.084))
file.close()
