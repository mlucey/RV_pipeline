from astropy.io import fits
import numpy as np
#import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join, isdir
from pyraf import iraf
from astropy.time import Time
iraf.rv()

#**need to edit obsrun and the paths**

obsrun = ['Apr2018','Mar2018','Dec2016','Feb2017','Feb2016','Mar2016']

def wav(image):
    return np.linspace(image[0].header['CRVAL1'],image[0].header['CRVAL1']+(image[0].header['NAXIS1']-1)*image[0].header['CDELT1'],image[0].header['NAXIS1'])

for t in obsrun:
    #where the images and raw data are
    inpath = '/Volumes/MADDIE/'+t+'/images/'
    #where you want the legend (which spectrum is run with which template) saved
    legendpath = '/Volumes/MADDIE/'+t+'/results/'
    #where you have the files.txt files that have which image corresponds with which raw image
    filepath = '/Volumes/MADDIE/'+t+'/'
    #where you have the templates(the have the same start so I included the start here too)
    temppath = '/Volumes/MADDIE/newtemplates/'
    #where you want the chosen template fxcor outputs saved(all fxcor outputs for each data run should be in the same place)
    outpath = '/Volumes/MADDIE/'+t+'/results/donefxcor/'
    #where you want all fxcor outputs saved
    alloutpath = '/Volumes/MADDIE/'+t+'/results/alldonefxcor/'
    #where you want the halpha fxcor outputs saved
    halphaoutpath = '/Volumes/MADDIE/'+t+'/results/halphafxcor/'

    templatelist = np.loadtxt(temppath+'RV_singlesource.txt', dtype='string', skiprows=1, usecols=0, unpack=True)
    halphatemplatelist = np.loadtxt(temppath+'RV_halpha.txt', dtype='string', skiprows=1, usecols=0, unpack=True)
    rvs, rverrs = np.loadtxt(temppath+'RV_singlesource.txt', dtype='float', skiprows=1, usecols=(10,12), unpack=True)
    hrvs, hrverrs = np.loadtxt(temppath+'RV_halpha.txt', dtype='float', skiprows=1, usecols=(10,12), unpack=True)

    #feb2017 region of spectrum to correlate if you want to avoid halpha
    if t == 'Feb2017':
        region = "6425-6550,6570-6800"
    #normal region of spectrum to correlate if you want to avoid halpha
    else:
        region = "6300-6550,6570-6800"

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
                flux = file[0].data[i][200:1800]
                errs = np.sqrt(abs(flux))
                imagewav = wav(file)
                objwav = imagewav[200:1800]
                listflux = flux.tolist()
                maxindex = listflux.index(max(flux))
                try:
                    for k in range(len(templatelist)):
                        try:
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+templatelist[k]+'.fits',apertures = str(int(aps[index])), function = "gaussian", background="INDEF", osample = region, rsample = region, output = alloutpath+images[j]+"_"+str(writeras[index])+'_all', verbose='txtonly', interactive='no')
                        except:
                            print(images[j]+' and '+templatelist[k]+' did not work')
                    height, fwhm = np.loadtxt(alloutpath+images[j]+"_"+str(writeras[index])+'_all.txt', dtype='string', usecols=(7,8), unpack=True)
                    height[height == 'INDEF'] = '0.0'
                    height = np.asfarray(height,float)
                    best = np.argmax(height)
                    if height[best] < .42:
                        for k in range(len(halphatemplatelist)):
                            try:
                                iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+halphatemplatelist[k]+'.fits',apertures = str(int(aps[index])), function = 'gaussian', background='INDEF', osample='*', rsample = '*', output = halphaoutpath+images[j]+'_'+str(writeras[index])+'_halpha', verbose='txtonly', interactive='no')
                            except:
                                print('Halpha '+images[j]+' and '+halphatemplatelist[k]+' did not work')
                        hheight, hfwhm = np.loadtxt(halphaoutpath+images[j]+"_"+str(writeras[index])+'_halpha.txt', dtype='string', usecols=(7,8), unpack=True)
                        hheight[hheight == 'INDEF'] = '0.0'
                        hheight = np.asfarray(hheight,float)
                        hbest = np.argmax(hheight)
                        if height[best] < hheight[hbest]:
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+halphatemplatelist[hbest]+'.fits',apertures = str(int(aps[index])), function = 'gaussian', background='INDEF', osample='*', rsample = '*', output = outpath+images[j]+'_'+str(writeras[index])+'', verbose='txtonly', interactive='no')
                            text.write("%8s %4s %8s %8.3f %8s \n" % (ids[index],ras[index],halphatemplatelist[hbest],hrvs[hbest],images[j]+'_'+aps[index]))
                            print(images[j],aps[index])
                            halphatext.write("%8s %4s %8s %8.3f %8s \n" % (ids[index],ras[index],templatelist[best],rvs[best],images[j]+'_'+aps[index]))
                        else:
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+templatelist[best]+'.fits',apertures = str(int(aps[index])), function = "gaussian", background="INDEF", osample = region, rsample =region, output = outpath+images[j]+"_"+str(writeras[index]), verbose='txtonly', interactive='no')
                            text.write("%8s %4s %8s %8.3f %8s \n" % (ids[index],ras[index],templatelist[best],rvs[best],images[j]+'_'+aps[index]))
                    else:
                        iraf.rv.fxcor(inpath+images[j]+'_nosky.fits',temppath+'new_'+templatelist[best]+'.fits',apertures = str(int(aps[index])), function = "gaussian", background="INDEF", osample = region, rsample =region, output = outpath+images[j]+"_"+str(writeras[index]), verbose='txtonly', interactive='no')
                        text.write("%8s %4s %8s %8.3f %8s \n" % (ids[index],ras[index],templatelist[best],rvs[best],images[j]+'_'+aps[index]))
                except:
                    text.write("%8s %4s %8s %4s \n" % (ids[index],ras[index],'fail','fail'))
            else:
                print(splits[0])
    text.close()
    halphatext.close()
