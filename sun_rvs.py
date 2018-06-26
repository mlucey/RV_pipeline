from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join
from pyraf import iraf

obsrun = 'Feb2016'
inpath = '/Volumes/MADDIE/'+obsrun+'/final_stuff/donefxcor/'

def removeb(input):
    if input[0] == 'b':
        split = input.split("'")
        return split[1]
    else:
        return input

 #getting the suns rv.
 #CHANGE TEMPLATE
fits = fits.open('/Volumes/MADDIE/'+obsrun+'/images/RED.ms.fits')
head = fits[0].header
numstars = head['NAXIS2']

aps = []
for i in range(numstars):
    splits = head['APNUM'+str(i+1)].split(' ')
    ap = splits[0]
    aps.append(ap)
for j in range(len(aps)):
    iraf.rv.fxcor('/Volumes/MADDIE/'+obsrun+'/images/RED.ms.fits' ,'/Volumes/MADDIE/RV_std/new_solar.fits', apertures=aps[j] ,output= '/Volumes/MADDIE/'+obsrun+'/final_stuff/donefxcor/phoenix_Reds'+str(aps[j]) , verbose = 'txtonly', interactive = 'no')
"""
fitsfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]

reds = []
aps = []
for i in range(len(fitsfiles)):
    if fitsfiles[i][0:4] == 'Reds':
        reds.append(fitsfiles[i])


temprvs = []
for j in reds:
    temprv = np.loadtxt(inpath+str(j), usecols = (11,))
    temprvs.append(temprv)
    reds1 = j.split("s")
    reds2 = reds1[1].split(".")
    ap = reds2[0]
    aps.append(ap)

rvs = []

for l in range(len(temprvs)):
    rv = -22.016+float(temprvs[l])
    rvs.append(rv)

fiber_rvs = open(inpath+obsrun+'_fiber_rvs.txt','w')
for i in range(len(rvs)):
    fiber_rvs.write("%2s %4f \n" % (aps[i],rvs[i]))
fiber_rvs.close()


"""

"""

badbadrv = np.loadtxt('/Volumes/MADDIE/Halpha/Dec2016/final_stuff/one_rv.dat', usecols = (2,),dtype = str)
ap = np.loadtxt('/Volumes/MADDIE/Halpha/Dec2016/final_stuff/one_rvcorrect.txt', usecols = (10,))
fiberrvs = np.loadtxt('/Volumes/MADDIE/Dec2016/final_stuff/donefxcor/Dec2016_fiber_rvs.txt', usecols =(1,))
fibers = np.loadtxt('/Volumes/MADDIE/Dec2016/final_stuff/donefxcor/Dec2016_fiber_rvs.txt',usecols =(0,))

badrv = []
for i in badbadrv:
    adbadrv = removeb(i)
    if adbadrv != 'INDEF':
        badrv.append(float(adbadrv))



rvs = []
for i in range(len(badrv)):
    index = (fibers.tolist()).index(ap[i])
    rvs.append(badrv[i]-fiberrvs[index])

praervs = []
for i in rvs:
    if 30< i <40:
        praervs.append(i)

print(np.mean(praervs))
"""
