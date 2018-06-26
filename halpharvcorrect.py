from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join

# creates a file to run in rvcorrect to get heliocentrically corrected rvs for the halpha stars
#After you have run it paste this line into a pyraf window: 
#rvcorrect f=halpharvcorrect.txt > halhparv.dat

obsrun = 'Apr2018'

inpath = '/Volumes/MADDIE/'
#location of files.txt with list of images and correspondng orig images
file = inpath+obsrun+'/files.txt'
#location of fxcor outputs (created in halphafxcor.py)
fxcorpath = inpath+obsrun+'/final_stuff/donefxcor'
#location of halpha legend file (created in halphafxcor.py)
legend = inpath+obsrun+'/Halphalegend.txt'
#where you want to write the halpharvcorrect.txt file
outpath = inpath+obsrun+'/final_stuff/halpharvcorrect.txt'
#where you have the reduced data
imagepath = inpath+obsrun+'/images/'

def removeb(input):
    if input[0] == 'b':
        split = input.split("'")
        return split[1]
    else:
        return input

images = np.loadtxt(file, usecols = (1,), dtype =str)
    
fitsfiles = [f for f in listdir(fxcorpath) if isfile(join(fxcorpath, f))]
rvtemp = np.loadtxt(legend, usecols = (2,))
bids = np.loadtxt(legend , usecols = (0,), dtype = str)

ids = []
for i in bids:
    ids.append(removeb(i))



text = open(outpath,'w')
files = []
for l in range(len(fitsfiles)):
    splits = fitsfiles[l].split('_')
    if splits[0][0] == 'h':
        files.append(fitsfiles[l])
dates = []
ras =[]
decs = []
uts = []

for i in range(len(images)):
    file = fits.open(imagepath+str(images[i])+'_nosky.fits')
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
    if len(parts) == 4:
        image = parts[1]+'_'+parts[2]
        ra = parts[3][0:8]
    else:
        image = parts[1]
        ra = parts[2][0:8]
    brvls = np.loadtxt(fxcorpath+'/'+str(files[m]),usecols = (11,),dtype=str)
    aps = np.loadtxt(fxcorpath+'/'+str(files[m]), usecols = (4,))
    berr = np.loadtxt(fxcorpath+'/'+str(files[m]), usecols = (13,), dtype = str)
    bhght = np.loadtxt(fxcorpath+'/'+str(files[m]), usecols = (7,), dtype = str)
    rvls = removeb(str(brvls))
    err = removeb(str(berr))
    hght = removeb(str(bhght))
    if rvls != 'INDEF':
        rv = float(rvls)
    if image in images:
        iindex = (images.tolist()).index(image)
        if ra in ids:
            index = ids.index(ra)
            if rvls != 'INDEF':
                text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[iindex][0],dates[iindex][1],dates[iindex][2],uts[iindex],ras[iidex],decs[iindex],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image))
            else:
                text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[iindex][0],dates[iindex][1],dates[iindex][2],uts[iindex],ras[iindex],decs[iindex],0,rvls,ra,ids[index],aps,err,hght,image))
text.close()


