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

obsrun = ['Apr2018','Mar2018','Dec2016','Feb2017','Feb2016','Mar2016']
cluster = 'prae'
first = '8'
#put in 8 for praesepe and 3 for pleiades
#location of obsrun folder
inpath = '/Volumes/MADDIE/'
#rvcorrect and rv.dat file location inside the obsrun folder
rvcorrect = '/final_stuff/'
#folder inside of obsrun folder that holds the raw data and the reduced data
imagefolder = '/images/'
#location of fxcor outputs (specifaically for the skyflats)
fxcorpath = '/final_stuff/donefxcor'
#first letter of skyflatfxcor ouputs (should be different than object ouputs)
firstskyflat = 'R'
#rv of the temp used for the skyflats (in km/s)
skyflattemp = 0
#where and what you want save the csv that prints one row per spectra as
outpathspec = inpath+'obs_'+cluster+'_catalog.csv'
#where and what you want save the csv that prints one row per star as
outpathstar = inpath+cluster+'_catalog.csv'

def number1(x):
    a = (ras[x],decs[x],int(ids[x]),rs[x],r_js[x])
    return list(a)

def onewrite(x):
    a = (obsrun+images[x],int(aps[x]),'{:.3f}'.format(float(hjds[x]-2400000)),'{:.3f}'.format(rvs[x]),'{:.3f}'.format(verrs[x]),'{:.3f}'.format(heights[x]))
    return list(a)

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

#loading in rvs and aps from rv.dat and rvcorrect.txt
for i in obsrun:
    baps = np.loadtxt(inpath+i+rvcorect+'rvcorrect.txt', usecols = (9,))
    #bhaps = np.loadtxt(inpath+i+rvcorrect+'halpharvcorrect.txt', usecols = (9,))
    rvs = np.loadtxt(inpath+i+rvcorrect+'rv.dat', usecols = (2,),dtype =str)
    #halpharvs = np.loadtxt(inpath+i+rvcorrect+'halpharv.dat', usecols = (2,),dtype = str)
    files = np.loadtxt(inpath+i+'/files.txt', usecols = (1,), dtype = str)
    origfiles = np.loadtxt(inpath+i+'/files.txt', usecols = (0,), dtype = str)
    #getting images,ra and heights from rvcorrect.txt
    bras = np.loadtxt(inpath+i+rvcorrect+'rvcorrect.txt', usecols = (8,),dtype = str)
    #bhras = np.loadtxt(inpath+i+rvcorrect+'halpharvcorrect.txt', usecols = (8,),dtype = str)
    bimages = np.loadtxt(inpath+i+rvcorrect+'rvcorrect.txt', usecols = (12,),dtype = str)
    print(bimages,i)
    #bhimages = np.loadtxt(inpath+i+rvcorrect+'halpharvcorrect.txt',usecols = (12,),dtype = str)
    bhghts = np.loadtxt(inpath+i+rvcorrect+'rvcorrect.txt', usecols = (11,),dtype = str)
    #bhhghts = np.loadtxt(inpath+i+rvcorrect+'halpharvcorrect.txt',usecols = (11,), dtype = str)
    #get hjds, and verrs from rv.dat and rvcorrect files
    bhjds = np.loadtxt(inpath+i+rvcorrect+'rv.dat', usecols = (0,))
    #bhhjds = np.loadtxt(inpath+i+rvcorrect+'halpharv.dat',usecols = (0,))
    bverrs = np.loadtxt(inpath+i+rvcorrect+'rvcorrect.txt', usecols = (10,),dtype = str)
    #bhverrs = np.loadtxt(inpath+i+rvcorrect+'halpharvcorrect.txt',usecols = (10,),dtype = str)
    #get ra, dec, id, rs, r-js, and aps from original images
    ogras = []
    checkras = []
    ogdecs = []
    ogrs = []
    ogr_js = []
    ogids = []
    for j in range(len(files)):
        file = fits.open(inpath+i+imagefolder+files[j]+'_nosky.fits')
        origfile = fits.open(inpath+i+imagefolder+origfiles[j]+'.fits')
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
    #grabbing all fxcor ouputs
    fitsfiles = [f for f in listdir(inpath+i+fxcorpath) if isfile(join(inpath+i+fxcorpath, f))]

    reds = []
    redaps = []
    #grabbing only skyflat outputs
    for k in range(len(fitsfiles)):
        if fitsfiles[k][0] == firstskyflat:
            reds.append(fitsfiles[k])
    temprvs = []
    #grabbing rv and ap information from fxcor skyflat outputs
    for j in reds:
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
    #removing info for stars that do not have a skyflat for their aperture and prints the ones that get cut
    for l in range(len(rvs)):
        if str(int(baps[l])) in redaps:
            index = redaps.index(str(int(baps[l])))
            aps.append(baps[l])
            images.append(bimages[l])
            ras.append(bras[l])
            hghts.append(bhghts[l])
            hjds.append(bhjds[l])
            verrs.append(bverrs[l])
            #subtracting skyflat rv from rv
            if rvs[l] != 'INDEF':
                #donerv.append(str("{0:.3f}".format(float(rvs[l]))))
                donerv.append(str("{0:.3f}".format(float(rvs[l])-float(temprvs[index])-float(skyflattemp)))
            else:
                donerv.append('INDEF')
        else:
            print(bimages[l],baps[l],i)

    donehrv = []
    halphaaps = []
    halphaimages = []
    halpharas = []
    halphahghts = []
    halphahjds = []
    halphaverrs = []
    """
    #does the same as above but for the halpha spectra
    for l in range(len(halpharvs)):
        if str(int(bhaps[l])) in redaps:
            index = redaps.index(str(int(bhaps[l])))
            halphaaps.append(bhaps[l])
            halphaimages.append(bhimages[l])
            halpharas.append(bhras[l])
            halphahghts.append(bhhghts[l])
            halphahjds.append(bhhjds[l])
            halphaverrs.append(bhverrs[l])
            if halpharvs[l] != 'INDEF':
                #donehrv.append(str("{0:.3f}".format(float(halpharvs[l]))))
                donehrv.append(str("{0:.3f}".format(float(halpharvs[l])-float(temprvs[index]-float(skyflattemp))))
            else:
                donehrv.append('INDEF')
    """
      
    nopeats = list(set(halpharas))
    hindexs = []
    seenonce = []
    seentwice = []
    seenthrice = []
    #gets indexs of repeats in halphras                          
    for m in nopeats:
        hindexs.append(list_duplicates_of(halpharas,m))
    #matching fxcor and info from raw data header                          
    for k in range(len(ras)):
        index = checkras.index(str(ras[k]))
        #using data from halpha if it was run in the halpha pipeline
        #still has indexing issue problem is there are repeats in ras[k] and in halpha, so having trouble matching them
        if ras[k] in halpharas:
            nindex = nopeats.index(ras[k])
            if ras[k] in seenonce and ras[k] not in seentwice and len(hindexs[nindex]) >1:
                allras.append(ogras[index])
                alldecs.append(ogdecs[index])
                allrs.append(ogrs[index])
                allr_js.append(ogr_js[index])
                allruns.append(i)
                allimages.append(halphaimages[hindexs[nindex][1]])
                allaps.append(halphaaps[hindexs[nindex][1]])
                allhjds.append(halphahjds[hindexs[nindex][1]])
                allrvs.append(donehrv[hindexs[nindex][1]])
                allerrs.append(halphaverrs[hindexs[nindex][1]])
                allhghts.append(halphahghts[hindexs[nindex][1]])
                seentwice.append(ras[k])
            if ras[k] in seentwice and ras[k] not in seenthrice and len(hindexs[nindex]) > 2:
                allras.append(ogras[index])
                alldecs.append(ogdecs[index])
                allrs.append(ogrs[index])
                allr_js.append(ogr_js[index])
                allruns.append(i)
                allimages.append(halphaimages[hindexs[nindex][2]])
                allaps.append(halphaaps[hindexs[nindex][2]])
                allhjds.append(halphahjds[hindexs[nindex][2]])
                allrvs.append(donehrv[hindexs[nindex][2]])
                allerrs.append(halphaverrs[hindexs[nindex][2]])
                allhghts.append(halphahghts[hindexs[nindex][2]])
                seenthrice.append(ras[k])
            if ras[k] in seenthrice and len(hindexs[nindex]) >3 and len(hindexs[nindex]) > 3:
                allras.append(ogras[index])
                alldecs.append(ogdecs[index])
                allrs.append(ogrs[index])
                allr_js.append(ogr_js[index])
                allruns.append(i)
                allimages.append(halphaimages[hindexs[nindex][3]])
                allaps.append(halphaaps[hindexs[nindex][3]])
                allhjds.append(halphahjds[hindexs[nindex][3]])
                allrvs.append(donehrv[hindexs[nindex][3]])
                allerrs.append(halphaverrs[hindexs[nindex][3]])
                allhghts.append(halphahghts[hindexs[nindex][3]])
            if ras[k] not in seenonce:
                allras.append(ogras[index])
                alldecs.append(ogdecs[index])
                allrs.append(ogrs[index])
                allr_js.append(ogr_js[index])
                allruns.append(i)
                allimages.append(halphaimages[hindexs[nindex][0]])
                allaps.append(halphaaps[hindexs[nindex][0]])
                allhjds.append(halphahjds[hindexs[nindex][0]])
                allrvs.append(donehrv[hindexs[nindex][0]])
                allerrs.append(halphaverrs[hindexs[nindex][0]])
                allhghts.append(halphahghts[hindexs[nindex][0]])
                seenonce.append(ras[k])
        else:
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

#write to csv file this does one row per spectra                              
with open(outpathspec,'w') as f:
    writer = csv.writer(f)
    writer.writerow(['RA','DEC','r','r-J','obsrun','image','ap','hjd','rv','err','hght'])
    for i in range(len(allras)):
        writer.writerow([allras[i],alldecs[i],allrs[i],allr_js[i],allruns[i],allimages[i],allaps[i],"{0:.3f}".format(allhjds[i]-2400000),allrvs[i],allerrs[i],allhghts[i]])

norepeats = list(set(allras))



#write to csv file, one row per star
with open(outpathstar,'w') as f:
    writer = csv.writer(f)
    for i in range(len(norepeats)):
        indexs = list_duplicates_of(allras,norepeats[i])
        row = [allras[indexs[0]],alldecs[indexs[0]],allrs[indexs[0]],allr_js[indexs[0]],allruns[indexs[0]],allimages[indexs[0]],allaps[indexs[0]],"{0:.3f}".format(allhjds[indexs[0]]-2400000),allrvs[indexs[0]],allerrs[indexs[0]],allhghts[indexs[0]]]
        for j in range(1,len(indexs)-1):
            row = row+[allruns[indexs[j]],allimages[indexs[j]],allaps[indexs[j]],"{0:.3f}".format(allhjds[indexs[j]]-2400000),allrvs[indexs[j]],allerrs[indexs[j]],allhghts[indexs[j]]]
        writer.writerow(row)

