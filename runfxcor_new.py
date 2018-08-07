from astropy.io import fits
import numpy as np
#import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join, isdir
from pyraf import iraf
from astropy.time import Time
iraf.rv()

obsrun = 'Mar2018'
#where the images and raw data are
inpath = '/Volumes/MADDIE/'+obsrun+'/images/'
#where you want the legend (which spectrum is run with which template) saved
legendpath = '/Volumes/MADDIE/'+obsrun+'/results/'
#where you have the files.txt files that have which image corresponds with which raw image
filepath = '/Volumes/MADDIE/'+obsrun+'/'
#where you have the templates(the have the same start so I included the start here too)
temppath = '/Volumes/MADDIE/newtemplates/'
#where you want the chosen template fxcor outputs saved(all fxcor outputs for each data run should be in the same place)
outpath = '/Volumes/MADDIE/'+obsrun+'/results/donefxcor/'
#where you want all fxcor outputs saved
alloutpath = '/Volumes/MADDIE/'+obsrun+'/results/alldonefxcor/'
#where you want the halpha fxcor outputs saved
halphaoutpath = '/Volumes/MADDIE/'+obsrun+'/results/halphafxcor/'

templatelist = np.loadtxt(temppath+'RV_singlesource.txt', dtype='string', skiprows=1, usecols=0, unpack=True)
halphatemplatelist = np.loadtxt(temppath+'RV_halpha.txt', dtype='string', skiprows=1, usecols=0, unpack=True)
rvs, rverrs = np.loadtxt(temppath+'RV_singlesource.txt', dtype='float', skiprows=1, usecols=(10,12), unpack=True)
hrvs, hrverrs = np.loadtxt(temppath+'RV_halpha.txt', dtype='float', skiprows=1, usecols=(10,12), unpack=True)


#feb2017 region of spectrum to correlate if you want to avoid halpha
#region = "6425-6550,6570-6800"
#normal region of spectrum to correlate if you want to avoid halpha
region = "6300-6550,6570-6800"

c=2.99792*10**5
halpha =6562.801

def wav(image):
    return np.linspace(image[0].header['CRVAL1'],image[0].header['CRVAL1']+(image[0].header['NAXIS1']-1)*image[0].header['CDELT1'],image[0].header['NAXIS1'])
#keeps track of which spectrum was run with which template
text = open(legendpath+'legend.txt','w')
#keeps track of which spectrum did not correlate the halpha region
halphatext = open(legendpath+'Halphalegend.txt','w')

images = np.loadtxt(filepath+'files.txt', usecols = (1,),dtype = str)
origimages = np.loadtxt(filepath+'files.txt', usecols = (0,),dtype = str)

for j in range(len(images)):
    file = fits.open(inpath+images[j]+'_nosky.fits')
    origfile = fits.open(inpath+origimages[j]+'.fits')
    numstars = file[0].header['NAXIS2']
    ras = []
    writeras = []
    ids = []
    aps = []
    rs =[]
    rows = []
    # pulling info out of headers
    for i in range(1,101):
        splits = (origfile[0].header['SLFIB'+str(i)]).split(" ")
        if splits[1] == "1":
            ra = splits[2]
            id = splits[4]
            r = splits[5]
            aps.append(i)
            ras.append(ra)
            ids.append(id)
            rs.append(r)
            splits = ra.split(':')
            writeras.append(splits[0]+splits[1]+splits[2])
    val= file[0]
    for i in range(numstars):
        #matching info from reduced data to raw data headers
        splits = (file[0].header['APID'+str(i+1)]).split(" ")
        if splits[0] in ids:
            index = ids.index(splits[0])
            try:
                #running the spectrum against every template in fxcor
                    for k in range(len(templatelist)):
                        try:
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+templatelist[k]+'.fits',apertures = str(int(aps[index])), function = "gaussian", background="INDEF", osample = region, rsample = region, output = alloutpath+images[j]+"_"+str(writeras[index])+'_all', verbose='txtonly', interactive='no')
                        except:
                            print(images[j]+' and '+templatelist[k]+' did not work')
                    #pulling the ccf heights and fwhm for each template out of fxcor output
                    height, fwhm = np.loadtxt(alloutpath+images[j]+"_"+str(writeras[index])+'_all.txt', dtype='string', usecols=(7,8), unpack=True)
                    height[height == 'INDEF'] = '0.0'
                    height = np.asfarray(height,float)
                    #finding the highest ccf peak (gives the best match)
                    best = np.argmax(height)
                    #if ccf peak is less than .5 it will not pass the quality check so we try halpha analysis
                    if height[best] < .5:
                        #running the spectrum against all the templates with halpha emission in fxcor
                        for k in range(len(halphatemplatelist)):
                            try:
                                iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+halphatemplatelist[k]+'.fits',apertures = str(int(aps[index])), function = 'gaussian', background='INDEF', osample='*', rsample = '*', output = halphaoutpath+images[j]+'_'+str(writeras[index])+'_halpha', verbose='txtonly', interactive='no')
                            except:
                                print('Halpha '+images[j]+' and '+halphatemplatelist[k]+' did not work')
                        # grabbing the ccf peak heights from the halpha analysis
                        hheight, hfwhm = np.loadtxt(halphaoutpath+images[j]+"_"+str(writeras[index])+'_halpha.txt', dtype='string', usecols=(7,8), unpack=True)
                        hheight[hheight == 'INDEF'] = '0.0'
                        hheight = np.asfarray(hheight,float)
                        #finding the highest ccf peak
                        hbest = np.argmax(hheight)
                        # if the ccf peak from halpha analysis is higher than the previous ccf peak we keep the halpha result
                        if height[best] < hheight[hbest]:
                            #final run with best template to be used in final analysis later
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+halphatemplatelist[hbest]+'.fits',apertures = str(int(aps[index])), function = 'gaussian', background='INDEF', osample='*', rsample = '*', output = outpath+images[j]+'_'+str(writeras[index])+'', verbose='txtonly', interactive='no')
                            #writing a legend to keep track of what the best template spectrum was
                            text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],halphatemplatelist[hbest],hrvs[hbest]))
                           #writing a legend to keep track of which results come from halpha analysis
                            halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],templatelist[best],rvs[best]))
                        else:
                          # if the ccf peak did not increase with halpha analysis we go back and use the best template match from normal analysis
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+templatelist[best]+'.fits',apertures = str(int(aps[index])), function = "gaussian", background="INDEF", osample = region, rsample =region, output = outpath+images[j]+"_"+str(writeras[index]), verbose='txtonly', interactive='no')
                            text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],templatelist[best],rvs[best]))
                    else:
                        #if the ccf peak is greater than 0.5 it passes that quality check and we use it in our final analysis
                        iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+templatelist[best]+'.fits',apertures = str(int(aps[index])), function = "gaussian", background="INDEF", osample = region, rsample =region, output = outpath+images[j]+"_"+str(writeras[index]), verbose='txtonly', interactive='no')
                        text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],templatelist[best],rvs[best]))
            except:
                #prints fails if it fails to run
                text.write("%8s %4s %8s %4s \n" % (ids[index],ras[index],'fail','fail'))
        else:
            #prints ids of spectrum whose header info did not make it into the reduced header info (from mixup in dot iraf files)
            print(splits[0])
text.close()
halphatext.close()
