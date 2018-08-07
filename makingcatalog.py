from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join, isdir
import itertools
import csv

# making catalog
#make sure that you renamed the old catalog file!!!

#list of each obsrun
obsrun = ['Apr2018','Mar2018','Dec2016','Feb2017','Feb2016','Mar2016']
#cluster name for catalog
cluster = 'prae'
#first digit of ra (for seperating clusters)
first = '8'
#path to location of each obs run directory
inpath = '/Volumes/MADDIE/'
#path to the rvcorrect and rv.dat file once in the obsrun directory
rvpath = '/results/'
#path to the fits files once in obsrun directory
imagepath = '/images/'
#path to the final fxcor outputs once in obsrun directory
fxcorpath = '/results/donefxcor'


#creating lists to hold all info
allras = []
alldecs = []
allids = []
allrs = []
allr_js = []
allimages = []
allaps = []
allhjds = []
allrvs = []
allerrs = []
allhghts = []
allruns = []
halpharas = []
halphadecs = []
halphaids = []
halphars = []
halphar_js = []
halphaimages = []
halphaaps = []
halphahjds = []
halpharvs = []
halphaerrs = []
halphahghts = []
halpharuns = []



for i in obsrun:
    #loading in rvs and aps from rv.dat and rvcorrect.txt
    baps = np.loadtxt(inpath+i+rvpath+'rvcorrect.txt', usecols = (10,))
    rvs = np.loadtxt(inpath+i+rvpath+'rv.dat', usecols = (2,),dtype =str)
    files = np.loadtxt(inpath+i+'/files.txt', usecols = (1,), dtype = str)
    origfiles = np.loadtxt(inpath+i+'/files.txt', usecols = (0,), dtype = str)
    #getting images,ra and heights from rvcorrect.txt
    bras = np.loadtxt(inpath+i+rvpath+'rvcorrect.txt', usecols = (8,),dtype = str)
    bimages = np.loadtxt(inpath+i+rvpath+'rvcorrect.txt', usecols = (13,),dtype = str)
    bhghts = np.loadtxt(inpath+i+rvpath+'rvcorrect.txt', usecols = (12,),dtype = str)
    #get hjds, and verrs from rv.dat and rvcorrect files
    bhjds = np.loadtxt(inpath+i+rvpath+'rv.dat', usecols = (0,))
    bverrs = np.loadtxt(inpath+i+rvpath+'rvcorrect.txt', usecols = (11,),dtype = str)
    #get ra, dec, id, rs, r-js, and aps from original images
    ogras = []
    checkras = []
    ogdecs = []
    ogrs = []
    ogr_js = []
    ogids = []
    for j in range(len(files)):
        file = fits.open(inpath+i+imagepath+files[j]+'_nosky.fits')
        origfile = fits.open(inpath+i+imagepath+origfiles[j]+'.fits')
        head = origfile[0].header
        #get info out of original headers
        for k in range(1,101):
            splits = head['SLFIB'+str(k)].split(" ")
            if splits[1] == '1':
                ra = splits[2]
                dec = splits[3]
                r = splits[5]
                r_j = splits[6]
                id = splits[4]
                ogras.append(ra)
                splits = ra.split(":")
                checkras.append(splits[0]+splits[1]+splits[2])
                ogdecs.append(dec)
                ogrs.append(r)
                ogr_js.append(r_j)
                ogids.append(id)
    #makes a list of all the fxcor output files
    fitsfiles = [f for f in listdir(inpath+i+fxcorpath) if isfile(join(inpath+i+fxcorpath', f))]

    reds = []
    redaps = []
    #grabs all the fxcor outputs for the skyflats
    for k in range(len(fitsfiles)):
        if fitsfiles[k][0] == 'R':
            reds.append(fitsfiles[k])
    temprvs = []
    for j in reds:
        #grabs rv and ap info for all the skyflats
        temprv = np.loadtxt(inpath+i+fxcorpath+'/'+str(j), usecols = (11,))
        temprvs.append(temprv)
        reds1 = j.split("s")
        reds2 = reds1[1].split(".")
        ap = reds2[0]
        redaps.append(ap)
    donerv = []
    aps = []
    images = []
    ras = []
    hghts = []
    hjds = []
    verrs = []
    for l in range(len(rvs)):
        #makes sure there is a skyflat for each aperture used in the object pointings
        if str(int(baps[l])) in redaps:
            #matches the aperutre from the skyflat to the object
            index = redaps.index(str(int(baps[l])))
            aps.append(baps[l])
            images.append(bimages[l])
            ras.append(bras[l])
            hghts.append(bhghts[l])
            hjds.append(bhjds[l])
            verrs.append(bverrs[l])
            if rvs[l] != 'INDEF':
                #subtracts the velocity of the skyflat specturm of the aperture that matches the object spectum from that spectrums heliocentrically corrected rv
                donerv.append(str("{0:.3f}".format(float(rvs[l])-float(temprvs[index]))))
            else:
                donerv.append('INDEF')
        else:
            #prings the image, ap and observational run of spectrum that do not have a skyflat for that aperture.
            print(bimages[l],baps[l],i)
    for k in range(len(ras)):
        #matches orig header info and rvcorrect info
        index = checkras.index(str(ras[k]))
        #writing all info into one big list for each colomn in final catalog
        allras.append(ogras[index])
        alldecs.append(ogdecs[index])
        allrs.append(ogrs[index])
        allr_js.append(ogr_js[index])
        allruns.append(i)
        allimages.append(images[k])
        allaps.append(aps[k])
        allhjds.append(hjds[k])
        allrvs.append(donerv[k])
        allerrs.append(verrs[k])
        allhghts.append(hghts[k])


#final catalog with every spectrum as it's own row
with open(inpath+'halphaobs_'+cluster+'_catalog.csv','w') as f:
    writer = csv.writer(f)
    writer.writerow(['RA','DEC','r','r-J','obsrun','image','ap','hjd','rv','err','hght'])
    for i in range(len(allras)):
        writer.writerow([allras[i],alldecs[i],allrs[i],allr_js[i],allruns[i],allimages[i],allaps[i],"{0:.3f}".format(allhjds[i]-2400000),allrvs[i],allerrs[i],allhghts[i]])

