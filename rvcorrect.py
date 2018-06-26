from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join

# rvcorrect.txt with ap and error
#CHANGE NUM. OF LOOPS
#origfilestart =

"""
#Apr2018
images = ['PB1','PB2','PF2','PF3']
obsrun = 'Apr2018'
"""
inpath = '/Volumes/MADDIE/'
"""
#Feb2016
images = ['PF1_1','PF1_2','PF3_1','PF3_2','PF4','PF5','PB3','PF2_1','PF2_2']
obsrun = 'Feb2016'
"""
"""
#Feb2017
images = ['PF1','PF2','PB1']
obsrun = 'Feb2017'
"""
"""
#Mar2018
images = ['PB1','PB2','PB3','PF1']
obsrun = 'Mar2018'
"""
"""
#Dec2016
images = ['PF1','PF2','PB1']
obsrun = 'Dec2016'

"""
#Mar2016
images = ['PB1','PB2','PB3','PF1','PF2','PF4','PF5']
obsrun = 'Mar2016'

def removeb(input):
    if input[0] == 'b':
        split = input.split("'")
        return split[1]
    else:
        return input


fitsfiles = [f for f in listdir(inpath+obsrun+'/final_stuff/donefxcor') if isfile(join(inpath+obsrun+'/final_stuff/donefxcor', f))]
rvtemp = np.loadtxt(inpath+obsrun+'/final_stuff/legend.txt', usecols = (3,))
bids = np.loadtxt(inpath+obsrun+'/final_stuff/legend.txt' , usecols = (1,), dtype = str)

ids = []
for i in bids:
    ids.append(removeb(i))

for i in range(len(ids)):
    splits = ids[i].split(':')
    ids[i] = splits[0]+splits[1]+splits[2]


text = open(inpath+obsrun+'/final_stuff/rvcorrect.txt','w')
files = []
for l in range(len(fitsfiles)):
    splits = fitsfiles[l].split('_')
    if splits[0][0] == 'P':
        files.append(fitsfiles[l])
dates = []
ras =[]
decs = []
uts = []

for i in range(len(images)):
    file = fits.open('/Volumes/MADDIE/'+obsrun+'/images/'+str(images[i])+'_nosky.fits')
    head = file[0].header
    date = head['DATE-OBS'][0:10]
    pieces = date.split('-')
    dates.append(pieces)
    RA = head['RA']
    ras.append(RA)
    DEC = head['DEC']
    decs.append(DEC)
    ut= head['UT']
    uts.append(ut)
for m in range(len(files)):
    parts = files[m].split("_")
    if len(parts) == 3:
        image = parts[0]+'_'+parts[1]
        ra = parts[2][0:8]
    else:
        image = parts[0]
        ra = parts[1][0:8]
    brvls = np.loadtxt(inpath+obsrun+'/final_stuff/donefxcor/'+str(files[m]),usecols = (11,),dtype=str)
    aps = np.loadtxt(inpath+obsrun+'/final_stuff/donefxcor/'+str(files[m]), usecols = (4,))
    berr = np.loadtxt(inpath+obsrun+'/final_stuff/donefxcor/'+str(files[m]), usecols = (13,), dtype = str)
    bhght = np.loadtxt(inpath+obsrun+'/final_stuff/donefxcor/'+str(files[m]), usecols = (7,), dtype = str)
    rvls = removeb(str(brvls))
    err = removeb(str(berr))
    hght = removeb(str(bhght))
    if rvls != 'INDEF':
        rv = float(rvls)
    if image == images[0]:
        if ra in ids:
            index = ids.index(ra)
            if rvls != 'INDEF':
                text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[0][0],dates[0][1],dates[0][2],uts[0],ras[0],decs[0],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image ))
            else:
                text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[0][0],dates[0][1],dates[0][2],uts[0],ras[0],decs[0],0,rvls,ra,ids[index],aps,err,hght,image))

    if image == images[1]:
        if ra in ids:
            index = ids.index(ra)
            if rvls != 'INDEF':
                text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[1][0],dates[1][1],dates[1][2],uts[1],ras[1],decs[1],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image))
            else:
                text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[1][0],dates[1][1],dates[1][2],uts[1],ras[1],decs[1],0,rvls,ra,ids[index],aps,err,hght,image))

    if image == images[2]:
        if ra in ids:
            index = ids.index(ra)
            if rvls != 'INDEF':
                text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[2][0],dates[2][1],dates[2][2],uts[2],ras[2],decs[2],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image))
            else:
                text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[2][0],dates[2][1],dates[2][2],uts[2],ras[2],decs[2],0,rvls,ra,ids[index],aps,err,hght,image))

    if image == images[3]:
        if ra in ids:
            index = ids.index(ra)
            if rvls != 'INDEF':
                text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[3][0],dates[3][1],dates[3][2],uts[3],ras[3],decs[3],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image))
            else:
                text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[3][0],dates[3][1],dates[3][2],uts[3],ras[3],decs[3],0,rvls,ra,ids[index],aps,err,hght,image))

    if image == images[4]:
        if ra in ids:
            index = ids.index(ra)
            if rvls != 'INDEF':
                text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[4][0],dates[4][1],dates[4][2],uts[4],ras[4],decs[4],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image))
            else:
                text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[4][0],dates[4][1],dates[4][2],uts[4],ras[4],decs[4],0,rvls,ra,ids[index],aps,err,hght,image))

    if image == images[5]:
        if ra in ids:
            index = ids.index(ra)
            if rvls != 'INDEF':
                text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[5][0],dates[5][1],dates[5][2],uts[5],ras[5],decs[5],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image))
            else:
                text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[5][0],dates[5][1],dates[5][2],uts[5],ras[5],decs[5],0,rvls,ra,ids[index],aps,err,hght,image))

    if image == images[6]:
        if ra in ids:
            index = ids.index(ra)
            if rvls != 'INDEF':
                text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[6][0],dates[6][1],dates[6][2],uts[6],ras[6],decs[6],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image))
            else:
                text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[6][0],dates[6][1],dates[6][2],uts[6],ras[6],decs[6],0,rvls,ra,ids[index],aps,err,hght,image))
"""
    if image == images[7]:
        if ra in ids:
            index = ids.index(ra)
            if rvls != 'INDEF':
                text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[7][0],dates[7][1],dates[7][2],uts[7],ras[7],decs[7],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image))
            else:
                text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[7][0],dates[7][1],dates[7][2],uts[7],ras[7],decs[7],0,rvls,ra,ids[index],aps,err,hght,image))

    if image == images[8]:
        if ra in ids:
            index = ids.index(ra)
            if rvls != 'INDEF':
                text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[8][0],dates[8][1],dates[8][2],uts[8],ras[8],decs[8],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image))
            else:
                text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[8][0],dates[8][1],dates[8][2],uts[8],ras[8],decs[8],0,rvls,ra,ids[index],aps,err,hght,image))
"""
text.close()

#now run rvcorrect in pyraf >> rvcorrect f=rvcorrect.txt > rv.dat

"""
#creating catalog --Don't use anymore
bvhelio = np.loadtxt('/Volumes/MADDIE/Feb2017/final_stuff/rv.dat',usecols = (2,), dtype = str)
bra = np.loadtxt('/Volumes/MADDIE/Dec2016/final_stuff/rvcorrect.txt', usecols = (8,), dtype = str)
bids = np.loadtxt('/Volumes/MADDIE/Dec2016/final_stuff/rvcorrect.txt', usecols = (9,) , dtype = str)
fitsfiles = [f for f in listdir('/Volumes/MADDIE/Dec2016/final_stuff/donefxcor') if isfile(join('/Volumes/MADDIE/Dec2016/final_stuff/donefxcor', f))]
files = []
for m in fitsfiles:
    if m[0] =='P':
        files.append(m)

praetxt = open('/Volumes/MADDIE/Dec2016/final_stuff/praervs.txt','w')
pleitxt = open('/Volumes/MADDIE/Dec2016/final_stuff/pleirvs.txt', 'w')
vhelio = []
ra = []
ids = []

for i in range(len(bvhelio)):
    parts = files[i].split("_")
    pv = bvhelio[i].split("'")
    pra = bra[i].split("'")
    pids = bids[i].split("'")
    vhelio.append(pv[1])
    ra.append(pra[1])
    ids.append(pids[1])
    if parts[0][0:2] == 'Pl':
        if ra[i] == parts[2]:
            pleitxt.write("%6s %6s %8s %4s %5s\n" % (ra[i],ids[i],parts[0]+parts[1]))
        if ra[i] == parts[1]:
            pleitxt.write("%6s %6s %8s %4s %3s\n" % (ra[i],ids[i],parts[0]))
    else:
        if ra[i] == parts[2]:
            praetxt.write("%6s %6s %8s %4s %5s\n" % (ra[i],ids[i],parts[0]+parts[1]))
        if ra[i] == parts[1]:
            praetxt.write("%6s %6s %8s %4s %3s\n" % (ra[i],ids[i],parts[0]))
praetxt.close()
pleitxt.close()

"""
