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
#what you want to name the csv that prints one row per spe

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
                donehrv.append(str("{0:.3f}".format(float(halpharvs[l])-float(temprvs[index]))))
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
with open('/Volumes/MADDIE/allphoenixobs_'+cluster+'_catalog.csv','w') as f:
    writer = csv.writer(f)
    writer.writerow(['RA','DEC','r','r-J','obsrun','image','ap','hjd','rv','err','hght'])
    for i in range(len(allras)):
        writer.writerow([allras[i],alldecs[i],allrs[i],allr_js[i],allruns[i],allimages[i],allaps[i],"{0:.3f}".format(allhjds[i]-2400000),allrvs[i],allerrs[i],allhghts[i]])

norepeats = list(set(allras))




with open('/Volumes/MADDIE/allphoenix_'+cluster+'_catalog.csv','w') as f:
    writer = csv.writer(f)
    for i in range(len(norepeats)):
        indexs = list_duplicates_of(allras,norepeats[i])
        row = [allras[indexs[0]],alldecs[indexs[0]],allrs[indexs[0]],allr_js[indexs[0]],allruns[indexs[0]],allimages[indexs[0]],allaps[indexs[0]],"{0:.3f}".format(allhjds[indexs[0]]-2400000),allrvs[indexs[0]],allerrs[indexs[0]],allhghts[indexs[0]]]
        for j in range(1,len(indexs)-1):
            row = row+[allruns[indexs[j]],allimages[indexs[j]],allaps[indexs[j]],"{0:.3f}".format(allhjds[indexs[j]]-2400000),allrvs[indexs[j]],allerrs[indexs[j]],allhghts[indexs[j]]]
        writer.writerow(row)




"""

# make histos
#update avg formula
boldpraeavgras = np.loadtxt('/Volumes/MADDIE/new1avgpraervs.txt',usecols = (0,),dtype = str)
boldpraeavgrvs = np.loadtxt('/Volumes/MADDIE/new1avgpraervs.txt',usecols =(1,),dtype = str)

boldprae2avgras = np.loadtxt('/Volumes/MADDIE/new2avgpraervs.txt',usecols = (0,),dtype = str)
boldprae2avgrvs  = np.loadtxt('/Volumes/MADDIE/new2avgpraervs.txt',usecols =(1,),dtype = str)

boldprae3avgras = np.loadtxt('/Volumes/MADDIE/new3avgpraervs.txt', usecols = (0,),dtype = str)
boldprae3avgrvs = np.loadtxt('/Volumes/MADDIE/new3avgpraervs.txt', usecols = (1,),dtype = str)

boldpraeras = np.loadtxt('/Volumes/MADDIE/newnoavgpraervs.txt',usecols = (0,),dtype = str)
boldpraervs = np.loadtxt('/Volumes/MADDIE/newnoavgpraervs.txt',usecols =(1,),dtype = str)

newnoavgpraervs = open('/Volumes/MADDIE/allnoavgpraervs.txt','w')
avgpraervs = open('/Volumes/MADDIE/all1avgpraervs.txt','w')
doubavgpraervs = open('/Volumes/MADDIE/all2avgpraervs.txt','w')
tripavgpraervs = open('/Volumes/MADDIE/all3avgpraervs.txt','w')
quadavgpraervs = open('/Volumes/MADDIE/all4avgpraervs.txt','w')
fifavgpraervs = open('/Volumes/MADDIE/all5avgpraervs.txt','w')
sixavgpraervs = open('/Volumes/MADDIE/all6avgpraervs.txt','w')
sevavgpraervs = open('/Volumes/MADDIE/all7avgpraervs.txt','w')
eigavgpraervs = open('/Volumes/MADDIE/all8avgpraervs.txt','w')
ninavgpraervs = open('/Volumes/MADDIE/all9avgpraervs.txt','w')
oldpraeras = []
oldpraervs = []
oldavgpraeras = []
oldavgpraervs = []
praedata = []
praestd = []
old2avgpraervs = []
old2avgpraeras = []
old3avgpraervs = []
old3avgpraeras = []

for i in range(len(boldpraeras)):
    if boldpraeras[i][0] == 'b':
        poldpraeras = boldpraeras[i].split("'")
        oldpraeras.append((poldpraeras[1]))
    else:
        oldpraeras.append(boldpraeras[i])
    if boldpraervs[i][0] == 'b':
        poldpraervs = boldpraervs[i].split("'")
        oldpraervs.append((poldpraervs[1]))
    else:
        oldpraervs.append(boldpraervs[i])

for i in range(len(boldprae2avgras)):
    if boldprae2avgras[i][0] == 'b':
        poldprae2avgras = boldprae2avgras[i].split("'")
        old2avgpraeras.append((poldprae2avgras[1]))
    else:
        old2avgpraeras.append(boldprae2avgras[i])
    if boldprae2avgrvs[i][0] == 'b':
        poldprae2avgrvs = boldprae2avgrvs[i].split("'")
        old2avgpraervs.append((poldprae2avgrvs[1]))
    else:
        old2avgpraervs.append(boldprae2avgervs[i])

for i in range(len(boldpraeavgras)):
    if boldpraeavgras[i][0] == 'b':
        poldpraeavgras = boldpraeavgras[i].split("'")
        oldavgpraeras.append((poldpraeavgras[1]))
    else:
        oldavgpraeras.append(boldpraeavgras[i])
    if boldpraeavgrvs[i][0] == 'b':
        poldpraeavgrvs = boldpraeavgrvs[i].split("'")
        oldavgpraervs.append((poldpraeavgrvs[1]))
    else:
        oldavgpraervs.append(boldpraevgervs[i])

for i in range(len(boldprae3avgras)):
    if boldprae3avgras[i][0] == 'b':
        poldprae3avgras = boldprae3avgras[i].split("'")
        old3avgpraeras.append((poldprae3avgras[1]))
    else:
        old3avgpraeras.append(boldprae3avgras[i])
    if boldprae3avgrvs[i][0] == 'b':
        poldprae3avgrvs = boldprae3avgrvs[i].split("'")
        old3avgpraervs.append((poldprae3avgrvs[1]))
    else:
        old3avgpraervs.append(boldprae3avgervs[i])
#1avg- first avergae happening now - adding averages to previously unaveraged
for i in range(len(ras)):
    for j in range(len(oldpraeras)):
        if ras[i] == oldpraeras[j]:
            for x in range(len(peats4)):
                if ras[i] in peats4:
                    o= peats3.index(ras[i])
                    p = peats2.index(ras[i])
                    q = repeats.index(ras[i])
                    r = uniq.index(ras[i])
                    avg = np.mean([float(oldpraervs[j]),finished[index4[x]],finished[index3[o]],finished[index2[p]],finished[rindex[q]],finished[uindex[r]]])
                    std = np.std([float(oldpraervs[j]),finished[index4[x]],finished[index3[o]],finished[index2[p]],finished[rindex[q]],finished[uindex[r]]])
                    praedata.append(avg)
                    praestd.append(std)
                    fifavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                if ras[i] not in peats4:
                    if ras[i] in peats3:
                        o= peats3.index(ras[i])
                        p = peats2.index(ras[i])
                        q = repeats.index(ras[i])
                        r = uniq.index(ras[i])
                        avg = np.mean([float(oldpraervs[j]),finished[index3[o]],finished[index2[p]],finished[rindex[q]],finished[uindex[r]]])
                        std = np.std([float(oldpraervs[j]),finished[index3[o]],finished[index2[p]],finished[rindex[q]],finished[uindex[r]]])
                        praedata.append(avg)
                        praestd.append(std)
                        quadavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                    if ras[i] not in peats3:
                        if ras[i] in peats2:
                            p = peats2.index(ras[i])
                            q = repeats.index(ras[i])
                            r = uniq.index(ras[i])
                            avg = np.mean([float(oldpraervs[j]),finished[index2[p]],finished[rindex[q]],finished[uindex[r]]])
                            std = np.std([float(oldpraervs[j]),finished[index2[p]],finished[rindex[q]],finished[uindex[r]]])
                            praedata.append(avg)
                            praestd.append(std)
                            tripavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                        if ras[i] not in peats2:
                            if ras[i] in repeats:
                                q = repeats.index(ras[i])
                                r = uniq.index(ras[i])
                                avg = np.mean([float(oldpraervs[j]),finished[rindex[q]],finished[uindex[r]]])
                                std = np.std([float(oldpraervs[j]),finished[rindex[q]],finished[uindex[r]]])
                                praedata.append(avg)
                                praestd.append(std)
                                doubavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                            if ras[i] not in repeats:
                                if ras[i] in uniq:
                                    r = uniq.index(ras[i])
                                    avg = np.mean([float(oldpraervs[j]),finished[uindex[r]]])
                                    std = np.std([float(oldpraervs[j]),finished[uindex[r]]])
                                    praedata.append(avg)
                                    praestd.append(std)
                                    avgpraervs.write("%10s %4s \n" % (ras[i],avg))
                                if ras[i] not in uniq:
                                    print(ras[i])


#1avg - no new avg
for p in range(len(oldavgpraeras)):
    if oldavgpraeras[p] not in ras:
        if oldavgpraervs[p] != 'INDEF':
            praedata.append(oldavgpraervs[p])
            avgpraervs.write("%10s %4s \n" % (oldavgpraeras[p],oldavgpraervs[p]))
avgpraervs.close()

#2 avg - second average happening now - adding average to previously averaged once
for i in range(len(ras)):
    for j in range(len(oldavgpraeras)):
        if ras[i] == oldavgpraeras[j]:
            for x in range(len(peats4)):
                if ras[i] in peats4:
                    o= peats3.index(ras[i])
                    p = peats2.index(ras[i])
                    q = repeats.index(ras[i])
                    r = uniq.index(ras[i])
                    avg = (2*float(oldavgpraervs[j])+finsihed[x]+finished[index3[o]]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/7
                    praedata.append(avg)
                    sixavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                if ras[i] not in peats4:
                    if ras[i] in peats3:
                        o= peats3.index(ras[i])
                        p = peats2.index(ras[i])
                        q = repeats.index(ras[i])
                        r = uniq.index(ras[i])
                        avg = (2*oldavgpraervs[j]+finished[index3[o]]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/6
                        praedata.append(avg)
                        fifavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                    if ras[i] not in peats3:
                        if ras[i] in peats2:
                            p = peats2.index(ras[i])
                            q = repeats.index(ras[i])
                            r = uniq.index(ras[i])
                            avg = (2*oldavgpraervs[j]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/5
                            praedata.append(avg)
                            quadavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                        if ras[i] not in peats2:
                            if ras[i] in repeats:
                                q = repeats.index(ras[i])
                                r = uniq.index(ras[i])
                                avg = (2*oldavgpraervs[j]+finished[rindex[q]]+finished[uindex[r]])/4
                                praedata.append(avg)
                                trippavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                            if ras[i] not in repeats:
                                if ras[i] in uniq:
                                        r = uniq.index(ras[i])
                                        if finished[uindex[r]] != 'INDEF':
                                            avg = (2*oldavgpraervs[j]+finished[uindex[r]])/3
                                            praedata.append(avg)
                                            doubavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                                        else:
                                            praedata.append(oldpraervs[j])
                                            avgpraervs.write("%10s %4s \n" % (oldavgpraeras[j],oldavgpraervs[j]))
                                if ras[i] not in uniq:
                                    print(ras[i])

#2 avg - no new avg
for j in range(len(old2avgpraeras)):
    if old2avgpraeras[j] not in ras:
        if old2avgpraervs[j] != 'INDEF':
            praedata.append(old2avgpraervs[j])
            doubavgpraervs.write("%10s %4s \n" % (old2avgpraeras[j],old2avgpraervs[j]))
doubavgpraervs.close()
#3 avg - happening now
for i in range(len(ras)):
    for j in range(len(old2avgpraeras)):
        if ras[i] == old2avgpraeras[j]:
            for x in range(len(peats4)):
                if ras[i] in peats4:
                    o= peats3.index(ras[i])
                    p = peats2.index(ras[i])
                    q = repeats.index(ras[i])
                    r = uniq.index(ras[i])
                    avg = (3*old2avgpraervs[j]+finsihed[x]+finished[index3[o]]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/8
                    praedata.append(avg)
                    sevavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                if ras[i] not in peats4:
                    if ras[i] in peats3:
                        o= peats3.index(ras[i])
                        p = peats2.index(ras[i])
                        q = repeats.index(ras[i])
                        r = uniq.index(ras[i])
                        avg = (3*old2avgpraervs[j]+finished[index3[o]]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/7
                        praedata.append(avg)
                        sixavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                    if ras[i] not in peats3:
                        if ras[i] in peats2:
                            p = peats2.index(ras[i])
                            q = repeats.index(ras[i])
                            r = uniq.index(ras[i])
                            avg = (3*old2avgpraervs[j]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/6
                            praedata.append(avg)
                            fifavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                        if ras[i] not in peats2:
                            if ras[i] in repeats:
                                q = repeats.index(ras[i])
                                r = uniq.index(ras[i])
                                avg = (3*old2avgpraervs[j]+finished[rindex[q]]+finished[uindex[r]])/5
                                praedata.append(avg)
                                quadavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                            if ras[i] not in repeats:
                                if ras[i] in uniq:
                                        r = uniq.index(ras[i])
                                        avg = (3*old2avgpraervs[j]+finished[uindex[r]])/4
                                        praedata.append(avg)
                                        tripavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                                if ras[i] not in uniq:
                                    print(ras[i])
#3 avg no new avg
for j in range(len(old3avgpraeras)):
    if old3avgpraeras[j] not in ras:
        if old3avgpraervs[j] != 'INDEF':
            praedata.append(old3avgpraervs[j])
            tripavgpraervs.write("%10s %4s \n" % (old3avgpraeras[j],old3avgpraervs[j]))

tripavgpraervs.close()
#4 avg happening now - averaging things that were previously 3 averaged
for i in range(len(ras)):
    for j in range(len(old3avgpraeras)):
        if ras[i] == old3avgpraeras[j]:
            for x in range(len(peats4)):
                if ras[i] in peats4:
                    o= peats3.index(ras[i])
                    p = peats2.index(ras[i])
                    q = repeats.index(ras[i])
                    r = uniq.index(ras[i])
                    avg = (4*old3avgpraervs[j]+finsihed[x]+finished[index3[o]]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/9
                    praedata.append(avg)
                    praestd.append(std)
                    eigavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                if ras[i] not in peats4:
                    if ras[i] in peats3:
                        o= peats3.index(ras[i])
                        p = peats2.index(ras[i])
                        q = repeats.index(ras[i])
                        r = uniq.index(ras[i])
                        avg = (4*old3avgpraervs[j]+finished[index3[o]]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/8
                        praedata.append(avg)
                        sevavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                    if ras[i] not in peats3:
                        if ras[i] in peats2:
                            p = peats2.index(ras[i])
                            q = repeats.index(ras[i])
                            r = uniq.index(ras[i])
                            avg = (4*old3avgpraervs[j]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/7
                            praedata.append(avg)
                            sixavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                        if ras[i] not in peats2:
                            if ras[i] in repeats:
                                q = repeats.index(ras[i])
                                r = uniq.index(ras[i])
                                avg = (4*old3avgpraervs[j]+finished[rindex[q]]+finished[uindex[r]])/6
                                praedata.append(avg)
                                fifavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                            if ras[i] not in repeats:
                                if ras[i] in uniq:
                                        r = uniq.index(ras[i])
                                        avg = (4*old3avgpraervs[j]+finished[uindex[r]])/5
                                        praedata.append(avg)
                                        quadavgpraervs.write("%10s %4s \n" % (ras[i],avg))
                                if ras[i] not in uniq:
                                    print(ras[i])
#no avg - averaging things that were not previously in catalog
for q in range(len(ras)):
    if ras[q] not in oldpraeras:
        if ras[q][0] == '8':
            for x in range(len(peats4)):
                if ras[q] in peats4:
                    o= peats3.index(ras[q])
                    p = peats2.index(ras[q])
                    q = repeats.index(ras[q])
                    r = uniq.index(ras[q])
                    avg = (finsihed[x]+finished[index3[o]]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/5
                    praedata.append(avg)
                    quadavgpraervs.write("%10s %4s \n" % (ras[q],avg))
                if ras[q] not in peats4:
                    if ras[q] in peats3:
                        o= peats3.index(ras[q])
                        p = peats2.index(ras[q])
                        q = repeats.index(ras[q])
                        r = uniq.index(ras[q])
                        avg = (finished[index3[o]]+finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/4
                        praedata.append(avg)
                        tripavgpraervs.write("%10s %4s \n" % (ras[q],avg))
                    if ras[q] not in peats3:
                        if ras[q] in peats2:
                            p = peats2.index(ras[q])
                            q = repeats.index(ras[q])
                            r = uniq.index(ras[q])
                            avg = (finished[index2[p]]+finished[rindex[q]]+finished[uindex[r]])/3
                            praedata.append(avg)
                            doubavgpraervs.write("%10s %4s \n" % (ras[q],avg))
                        if ras[q] not in peats2:
                            if ras[q] in repeats:
                                q = repeats.index(ras[q])
                                r = uniq.index(ras[q])
                                avg = (finished[rindex[q]]+finished[uindex[r]])/2
                                praedata.append(avg)
                                avgpraervs.write("%10s %4s \n" % (ras[q],avg))
                            if ras[q] not in repeats:
                                if ras[q] in uniq:
                                        r = uniq.index(ras[q])
                                        praedata.append(finsihed[uindex[r]])
                                        quadavgpraervs.write("%10s %4s \n" % (ras[q],finsihed[uindex[r]]))
                                if ras[q] not in uniq:
                                    print(ras[q])

for p in range(len(oldpraeras)):
    if oldpraeras[p] not in ras:
        if oldpraervs[p] != 'INDEF':
            praedata.append(oldpraervs[p])
            newnoavgpraervs.write("%10s %4s \n" % (oldpraeras[p],oldpraervs[p]))
newnoavgpraervs.close()
praepoints = []
for l in range(len(praedata)):
    praepoints.append(float(praedata[l]))

plt.hist(praepoints, bins= range(-50,75))
cluster = []
for i in range(len(praepoints)):
    if 20 < praepoints[i] < 50:
        cluster.append(praepoints[i])
print(np.mean(cluster))
plt.show()

"""
#PLEIADES
"""
"""
"""
boldpleiavgras = np.loadtxt('/Volumes/MADDIE/feb2016avgpleirvs.txt',usecols = (0,),dtype = str)
boldpleiavgrvs = np.loadtxt('/Volumes/MADDIE/feb2016avgpleirvs.txt',usecols =(1,),dtype = str)
"""
"""
boldpleiras = np.loadtxt('/Volumes/MADDIE/dec2016noavgpleirvs.txt',usecols = (0,),dtype = str)
boldpleirvs = np.loadtxt('/Volumes/MADDIE/dec2016noavgpleirvs.txt',usecols =(1,),dtype = str)

newnoavgpleirvs = open('/Volumes/MADDIE/newnoavgpleirvs.txt','w')
avgpleirvs = open('/Volumes/MADDIE/new1avgpleirvs.txt','w')
"""
"""
doubavgpleirvs = open('/Volumes/MADDIE/2avgpleirvs.txt','w')
tripavgpleirvs = open('Volumes/MADDIE/3avgpleirvs.txt','w')
"""
"""

oldpleirvs = []
oldpleiras = []
for i in range(len(boldpleiras)):
    if boldpleiras[i][0] == 'b':
        poldpleiras = boldpleiras[i].split("'")
        oldpleiras.append((poldpleiras[1]))
    else:
        oldpleiras.append(boldpleiras[i])
    if boldpleirvs[i][0] == 'b':
        poldpleirvs = boldpleirvs[i].split("'")
        oldpleirvs.append((poldpleirvs[1]))
    else:
        oldpleirvs.append(boldpleirvs[i])


pleidata = []

#noavg
for p in range(len(oldpleiras)):
    if oldpleiras[p] not in ras:
        if oldpleirvs[p] != 'INDEF':
            pleidata.append(oldpleirvs[p])
            newnoavgpleirvs.write("%10s %4s \n" % (oldpleiras[p],oldpleirvs[p]))
newnoavgpleirvs.close()

#1avg - happening now
for i in range(len(ras)):
    for j in range(len(oldpleiras)):
        if ras[i] == oldpleiras[j]:
            if finished[i] != 'INDEF':
                avg = (float(oldpleirvs[j])+float(finished[i]))/2
                pleidata.append(avg)
                avgpleirvs.write("%10s %4s \n" % (ras[i],avg))
avgpleirvs.close()

pleipoints = []
for l in range(len(pleidata)):
    pleipoints.append(float(pleidata[l]))

plt.hist(pleipoints, bins= range(-50,75))
cluster = []
for i in range(len(pleipoints)):
    if -10 < pleipoints[i] < 20:
        cluster.append(pleipoints[i])
print(np.mean(cluster))
plt.show()

"""
