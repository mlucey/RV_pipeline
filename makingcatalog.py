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
#run in python 2 to avoid formatting issues

obsrun = ['Apr2018','Mar2018','Dec2016','Feb2017','Feb2016','Mar2016']
cluster = 'prae'
first = '8'
#put in 8 for praesepe and 3 for pleiades

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

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
allhalpha = []


#loading in rvs and aps from rv.dat and rvcorrect.txt
for i in obsrun:
    halphaimaps = np.loadtxt('/Volumes/MADDIE/'+i+'/results/Halphalegend.txt', usecols = (4,),dtype=str)
    ohalpharas = np.loadtxt('/Volumes/MADDIE/'+i+'/results/Halphalegend.txt', usecols = (1,),dtype=str)
    onehalpharas = [s.replace(':', '') for s in ohalpharas]
    allhalpha.append(onehalpharas)
    baps = np.loadtxt('/Volumes/MADDIE/'+i+'/results/rvcorrect.txt', usecols = (10,))
    rvs = np.loadtxt('/Volumes/MADDIE/'+i+'/results/rv.dat', usecols = (2,),dtype =str)
    files = np.loadtxt('/Volumes/MADDIE/'+i+'/files.txt', usecols = (1,), dtype = str)
    origfiles = np.loadtxt('/Volumes/MADDIE/'+i+'/files.txt', usecols = (0,), dtype = str)
    #getting images,ra and heights from rvcorrect.txt
    bras = np.loadtxt('/Volumes/MADDIE/'+i+'/results/rvcorrect.txt', usecols = (8,),dtype = str)
    bimages = np.loadtxt('/Volumes/MADDIE/'+i+'/results/rvcorrect.txt', usecols = (13,),dtype = str)
    bhghts = np.loadtxt('/Volumes/MADDIE/'+i+'/results/rvcorrect.txt', usecols = (12,),dtype = str)
    #get hjds, and verrs from rv.dat and rvcorrect files
    bhjds = np.loadtxt('/Volumes/MADDIE/'+i+'/results/rv.dat', usecols = (0,))
    bverrs = np.loadtxt('/Volumes/MADDIE/'+i+'/results/rvcorrect.txt', usecols = (11,),dtype = str)
    #get ra, dec, id, rs, r-js, and aps from original images
    ogras = []
    checkras = []
    ogdecs = []
    ogrs = []
    ogr_js = []
    ogids = []
    for j in range(len(files)):
        file = fits.open('/Volumes/MADDIE/'+i+'/images/'+files[j]+'_nosky.fits')
        origfile = fits.open('/Volumes/MADDIE/'+i+'/images/'+origfiles[j]+'.fits')
        head = origfile[0].header
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
    fitsfiles = [f for f in listdir('/Volumes/MADDIE/'+i+'/results/donefxcor') if isfile(join('/Volumes/MADDIE/'+i+'/results/donefxcor', f))]

    reds = []
    redaps = []
    for k in range(len(fitsfiles)):
        if fitsfiles[k][0] == 'R':
            reds.append(fitsfiles[k])
    temprvs = []
    for j in reds:
        temprv = np.loadtxt('/Volumes/MADDIE/'+i+'/results/donefxcor/'+str(j), usecols = (11,))
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
        if str(int(baps[l])) in redaps:
            index = redaps.index(str(int(baps[l])))
            aps.append(baps[l])
            images.append(bimages[l])
            ras.append(bras[l])
            hghts.append(bhghts[l])
            hjds.append(bhjds[l])
            verrs.append(bverrs[l])
            if rvs[l] != 'INDEF':
                #donerv.append(str("{0:.3f}".format(float(rvs[l]))))
                donerv.append(str("{0:.3f}".format(float(rvs[l])-float(temprvs[index]))))
            else:
                donerv.append('INDEF')
        else:
            print(bimages[l],baps[l],i)
    for k in range(len(ras)):
        index = checkras.index(str(ras[k]))
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
        if str(images[k])+'_'+str(int(aps[k])) in halphaimaps:
            halpharas.append(ogras[index])
            halphadecs.append(ogdecs[index])
            halphars.append(ogrs[index])
            halphar_js.append(ogr_js[index])
            halpharuns.append(i)
            halphaimages.append(images[k])
            halphaaps.append(aps[k])
            halphahjds.append(hjds[k])
            halpharvs.append(donerv[k])
            halphaerrs.append(verrs[k])
            halphahghts.append(hghts[k])

with open('/Volumes/MADDIE/halphaobs_'+cluster+'_catalog.csv','w') as f:
    writer = csv.writer(f)
    writer.writerow(['RA','DEC','r','r-J','obsrun','image','ap','hjd','rv','err','hght'])
    for i in range(len(allras)):
        writer.writerow([allras[i],alldecs[i],allrs[i],allr_js[i],allruns[i],allimages[i],allaps[i],"{0:.3f}".format(allhjds[i]-2400000),allrvs[i],allerrs[i],allhghts[i]])



with open('/Volumes/MADDIE/halphaonly_'+cluster+'_catalog.csv','w') as f:
    writer = csv.writer(f)
    writer.writerow(['RA','DEC','r','r-J','obsrun','image','ap','hjd','rv','err','hght'])
    for i in range(len(halpharas)):
        writer.writerow([halpharas[i],halphadecs[i],halphars[i],halphar_js[i],halpharuns[i],halphaimages[i],halphaaps[i],"{0:.3f}".format(halphahjds[i]-2400000),halpharvs[i],halphaerrs[i],halphahghts[i]])




