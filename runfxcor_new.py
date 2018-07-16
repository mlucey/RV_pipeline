from astropy.io import fits
import numpy as np
#import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join, isdir
from pyraf import iraf
from astropy.time import Time
iraf.rv()

obsrun = 'Mar2016'
#where the images and raw data are
inpath = '/Users/Natalie/mann/fxcor_rv/'+obsrun+'/'
#where you want the legend (which spectrum is run with which template) saved
legendpath = '/Users/Natalie/mann/fxcor_rv/'+obsrun+'/results/'
#where you have the files.txt files that have which image corresponds with which raw image
filepath = '/Users/Natalie/mann/fxcor_rv/'+obsrun+'/'
#where you have the templates(the have the same start so I included the start here too)
temppath = '/Users/Natalie/mann/templates/'
#where you want the fxcor outputs saved(all fxcor outputs for each data run should be in the same place)
outpath = '/Users/Natalie/mann/fxcor_rv/'+obsrun+'/results/'
#the templates to correlate against

templatelist = np.loadtxt(temppath+'RV_singlesource.txt', dtype='string', skiprows=1, usecols=0, unpack=True)
rvs, rverrs = np.loadtxt(temppath+'RV_singlesource.txt', dtype='float', skiprows=1, usecols=(10,12), unpack=True)


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
halphatext = open(legendpath+'noHalphalgend.txt','w')

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
            flux = file[0].data[i][200:1800]
            errs = np.sqrt(abs(flux))
            imagewav = wav(file)
            objwav = imagewav[200:1800]
            listflux = flux.tolist()
            maxindex = listflux.index(max(flux))
            try:
                    for k in range(len(templatelist)):
                        try:
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+templatelist[k]+'.fits',apertures = str(int(aps[index])), function = "gaussian", background="INDEF", osample = region, rsample = region, output = outpath+images[j]+"_"+str(writeras[index])+'_all', verbose='txtonly', interactive='no')
                        except:
                            print(images[j]+' and '+templatelist[k]+' did not work')
                    height, fwhm = np.loadtxt(outpath+images[j]+"_"+str(writeras[index])+'_all.txt', dtype='string', usecols=(7,8), unpack=True)
                    height[height == 'INDEF'] = '0.0'
                    height = np.asfarray(height,float)
                    best = np.argmax(height)
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+templatelist[best]+'.fits',apertures = str(int(aps[index])), function = "gaussian", background="INDEF", osample = region, rsample =region, output = outpath+images[j]+"_"+str(writeras[index]), verbose='txtonly', interactive='no')
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],templatelist[best],rvs[best]))
               except:
                    text.write("%8s %4s %8s %4s \n" % (ids[index],ras[index],'fail','fail')
        else:
            print(splits[0])
text.close()
