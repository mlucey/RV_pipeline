from astropy.io import fits
import numpy as np
#import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join

# creates a file to run in rvcorrect to get heliocentrically corrected rvs
#**need to edit obsrun and the paths**

#After you have run it paste this line into a pyraf window when you are in the directory with rvcorrect.txt:
#rvcorrect f=rvcorrect.txt > rv.dat

obsrun = ['Apr2018','Mar2018','Dec2016','Feb2017','Feb2016','Mar2016']
def removeb(input):
    if input[0] == 'b':
        split = input.split("'")
        return split[1]
    else:
        return input

for t in obsrun:
    inpath = '/Volumes/MADDIE/'+t
    #location of files.txt with list of images and correspondng orig images
    file = inpath+'/files.txt'
    #location of fxcor outputs (created in runfxcor_new.py)
    fxcorpath = inpath+'/results/donefxcor'
    #location of halpha legend file (created in run_fxcor_new.py)
    legend = inpath+'/results/legend.txt'
    #where you want to write the rvcorrect.txt file
    outpath = inpath+'/results/rvcorrect.txt'
    #where you have the reduced data
    imagepath = inpath+'/images/'
    #if running in python 3, this will be necessary, removes a weird b that shows up in front of strings

    #getting list of images in obsrun
    images = np.loadtxt(file, usecols = (1), dtype =str)

    #getting all fxcor output files
    fitsfiles = [f for f in listdir(fxcorpath) if isfile(join(fxcorpath, f))]
    #pulling out of the legend the rv of the template used for each spectra
    rvtemp = np.loadtxt(legend, usecols = (3))
    #pulling ids out of the legend
    bids = np.loadtxt(legend , usecols = (0), dtype = str)
    #pulling image_ap out of legend
    bimaps = np.loadtxt(legend, usecols =(4), dtype = str)
    #removes bs if they are there (bs appear from a weird formatting issue when running in python 3)
    ids = []
    imaps = []
    for i in range(len(bids)):
        ids.append(removeb(bids[i]))
        imaps.append(removeb(bimaps[i]))


    #create rvcorrect file
    text = open(outpath,'w')
    #grabbing only object fxcor output files
    files = []
    for l in range(len(fitsfiles)):
        splits = fitsfiles[l].split('_')
        if splits[0][0] == 'P':
            files.append(fitsfiles[l])

    dates = []
    ras =[]
    decs = []
    uts = []
    #getting info for rvcorrect out of headers
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

    #getting info out of fxcor output files
    for m in range(len(files)):
        parts = files[m].split("_")
        if len(parts) == 3:
            image = parts[0]+'_'+parts[1]
            ra = parts[2][0:8]
        else:
            image = parts[0]
            ra = parts[1][0:8]
        brvls = np.loadtxt(fxcorpath+'/'+str(files[m]),usecols = (11,),dtype=str)
        aps = np.loadtxt(fxcorpath+'/'+str(files[m]), usecols = (4,))
        berr = np.loadtxt(fxcorpath+'/'+str(files[m]), usecols = (13,), dtype = str)
        bhght = np.loadtxt(fxcorpath+'/'+str(files[m]), usecols = (7,), dtype = str)
        rvls = removeb(str(brvls))
        err = removeb(str(berr))
        hght = removeb(str(bhght))
        imap = str(image)+'_'+str(int(aps))
        if rvls != 'INDEF':
            rv = float(rvls)
        #matching header and fxcor info
        if image in images:
            iindex = (images.tolist()).index(image)
            #matching fxcor and legend info
            if imap in imaps:
                index = imaps.index(imap)
                #writing to rvcorrect file and subtracting out template rv
                if rvls != 'INDEF':
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %.3f %4s\n" % (dates[iindex][0],dates[iindex][1],dates[iindex][2],uts[iindex],ras[iindex],decs[iindex],0,(rv+rvtemp[index]),ra,ids[index],aps,float(err),float(hght),image))
                else:
                    text.write("%4s %2s %2s %8s %8s %8s %1s %5s %8s %4s %3s %6s %5s %4s\n" % (dates[iindex][0],dates[iindex][1],dates[iindex][2],uts[iindex],ras[iindex],decs[iindex],0,rvls,ra,ids[index],aps,err,hght,image))
    text.close()
