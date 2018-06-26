from astropy.io import fits
import os
from os import listdir
from os.path import isfile, join, isdir
from array import array
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import scipy
import scipy.optimize as op
import scipy.ndimage
import find_rv
from astropy.time import Time

origimages= ['23','35','53']
origfilestart =
images = ['PF1','PF2','PB1']
obsrun = 'Feb2017'
inpath = '/Volumes/MADDIE/'+obsrun+'/images/'
outpath = '/Volumes/MADDIE/'+obsrun+'/final_stuff/'
temppath = '/Volumes/MADDIE/Forfxcor/'
c=2.99792*10**5
halpha =6562.801
crosscorr_width = 200
tempimages = ['1514217i','1629119i','1647952i','1649165i','1650473i','1688171i','1741422i','1756192p','1787492i','1789183i']
rvs = [22.190,14.531,-71.084,-9.499,-26.417,-22.016,22.942,18.210,-23.118,-45.630]
rverrs = [.287,0.050,.125,.032,.086,.108,.055,.030,.162,.056]

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def dopplershift(shift):
    return c*((shift/halpha)**2-1)/(1+(shift/halpha)**2)

def dopplershifterror(shift,error):
    return c*(4*shift*halpha**2)/((halpha**2+shift**2)**2)*error

def wav(image):
    return np.linspace(image[0].header['CRVAL1'],image[0].header['CRVAL1']+(image[0].header['NAXIS1']-1)*image[0].header['CDELT1'],image[0].header['NAXIS1'])


def halphatempinfo(temp):
    path = temppath+temp+'plushalpha.fits'
    tempfile = fits.open(path)
    flux = tempfile[0].data
    oldpath = temppath+temp+'.fits'
    oldtempfile = fits.open(oldpath)
    oldflux = oldtempfile[0].data[1]
    noise = oldtempfile[0].data[2]
    snr = np.mean(oldflux/noise)
    err = flux/snr
    return wav(tempfile)[200:1800],flux[200:1800],err[200:1800]

def tempinfo(temp):
    path = temppath+'new1_'+temp+'.fits'
    tempfile = fits.open(path)
    flux = tempfile[0].data
    oldpath = temppath+temp+'.fits'
    oldtempfile = fits.open(oldpath)
    oldflux = oldtempfile[0].data[1]
    noise = oldtempfile[0].data[2]
    snr = np.mean(oldflux/noise)
    err = flux/snr
    return wav(tempfile)[200:1800],flux[200:1800],err[200:1800]

#halpha = open(outpath+'halpha.txt','w+')

text = open(outpath+'rvcorrect.txt','w+')
"""
folders = [f for f in listdir(inpath) if isdir(join(inpath, f))]
datafolders  = []
for i in folders:
    if i[0]=='P':
        if i[1] != '1':
            datafolders.append(i)
"""
#Identifying H-alpha spectra
for j in range(len(images)):
    #feb132016 origfile = fits.open('/Volumes/MADDIE/'+obsrun+'/WIYN/'+origfilestart+str(origimages[j])+'.fits')
    file = fits.open(inpath+images[j]+'_nosky.fits')
    origfile = fits.open(inpath+origfilestart+origimages[j]+'.fits')
    numstars = file[0].header['NAXIS2']
    ras = []
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
    val= file[0]
    for i in range(len(numstars)):
        #pulling infor for rvcorrect
        jddate = file[0].header['JD']
        t = Time(float(jddate), format='jd', scale='utc')
        date = t.iso
        date1 = date.split(" ")
        dates = date1[0].split("-")
        ut = file[0].header['UT']
        ra = file[0].header['RA']
        dec = file[0].header['DEC']
        splits = (file[0].header['APID'+str(i+1)]).split(",")
        if splits[0] in ids:
            index = ids.index(splits[0])
            flux = file[0].data[i][200:1800]
            errs = np.sqrt(abs(flux))
            imagewav = wav(file)
            objwav = imagewav[200:1800]
            listflux = flux.tolist()
            maxindex = listflux.index(max(flux))
            if 6560 < objwav[maxindex] <6570:
                if max(flux) >3.25*np.std(flux)+np.mean(flux):
                    hflux = []
                    herrs = []
                    hobjwav = []
                    start = (np.abs(objwav-halpha-5)).argmin()
                    stop = (np.abs(objwav-halpha+5)).argmin()
                    for l in range(start):
                        hflux.append(flux[l])
                        herrs.append(errs[l])
                        hobjwav.append(objwav[l])
                    for k in range(stop,len(objwav)):
                        hflux.append(flux[k])
                        herrs.append(errs[k])
                        hobjwav.append(objwav[k])
                    #K0
                    if 11 <= float(rs[i]) < 12:
                        rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[9])[0],tempinfo(tempimages[9])[1],tempinfo(tempimages[9])[2],rvs[9],rverrs[9],crosscorr_width,outpath+name+'_'+str(aps[i]))
                        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'K0'))
                    #K3
                    if 12 <= float(rs[i]) < 13:
                        rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[4])[0],tempinfo(tempimages[4])[1],tempinfo(tempimages[4])[2],rvs[4],rverrs[4],crosscorr_width,outpath+name+'_'+str(aps[i]))
                        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'K3'))
                    #K5
                    if 13 <= float(rs[i]) < 13.5:
                        rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[8])[0],tempinfo(tempimages[8])[1],tempinfo(tempimages[8])[2],rvs[8],rverrs[8],crosscorr_width,outpath+name+'_'+str(aps[i]))
                        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'K5'))
                    #K7
                    if 13.5 <= float(rs[i]) < 14.25:
                        rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[0])[0],tempinfo(tempimages[0])[1],tempinfo(tempimages[0])[2],rvs[0],rverrs[0],crosscorr_width,outpath+name+'_'+str(aps[i]))
                        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'K7'))
                    #M0
                    if 14.25 <= float(rs[i]) <15:
                        rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[1])[0],tempinfo(tempimages[1])[1],tempinfo(tempimages[1])[2],rvs[1],rverrs[1],crosscorr_width,outpath+name+'_'+str(aps[i]))
                        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'M0'))
                    #M1
                    if 15 <= float(rs[i]) < 15.5:
                        rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[3])[0],tempinfo(tempimages[3])[1],tempinfo(tempimages[3])[2],rvs[3],rverrs[3],crosscorr_width,outpath+name+'_'+str(aps[i]))
                        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'M1'))
                    #M2
                    if 15.5 <= float(rs[i]) <16.25:
                        rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[2])[0],tempinfo(tempimages[2])[1],tempinfo(tempimages[2])[2],rvs[2],rverrs[2],crosscorr_width,outpath+name+'_'+str(aps[i]))
                        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'M2'))
                    #M3
                    if 16.25 <= float(rs[i]) < 17:
                        rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[6])[0],tempinfo(tempimages[6])[1],tempinfo(tempimages[6])[2],rvs[6],rverrs[6],crosscorr_width,outpath+name+'_'+str(aps[i]))
                        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'M3'))
                    #M3.5
                    if 17<= float(rs[i]) <18:
                        rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[7])[0],tempinfo(tempimages[7])[1],tempinfo(tempimages[7])[2],rvs[7],rverrs[7],crosscorr_width,outpath+name+'_'+str(aps[i]))
                        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'M3.5'))
            else:
                #K0
                if 11 <= float(rs[i]) < 12:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[9])[0],tempinfo(tempimages[9])[1],tempinfo(tempimages[9])[2],rvs[9],rverrs[9],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'K0'))
                #K3
                if 12 <= float(rs[i]) < 13:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[4])[0],tempinfo(tempimages[4])[1],tempinfo(tempimages[4])[2],rvs[4],rverrs[4],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'K3'))
                #K5
                if 13 <= float(rs[i]) < 13.5:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[8])[0],tempinfo(tempimages[8])[1],tempinfo(tempimages[8])[2],rvs[8],rverrs[8],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'K5'))
                #K7
                if 13.5 <= float(rs[i]) < 14.25:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[0])[0],tempinfo(tempimages[0])[1],tempinfo(tempimages[0])[2],rvs[0],rverrs[0],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'K7'))
                #M0
                if 14.25 <= float(rs[i]) <15:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[1])[0],tempinfo(tempimages[1])[1],tempinfo(tempimages[1])[2],rvs[1],rverrs[1],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'M0'))
                #M1
                if 15 <= float(rs[i]) < 15.5:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[3])[0],tempinfo(tempimages[3])[1],tempinfo(tempimages[3])[2],rvs[3],rverrs[3],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'M1'))
                #M2
                if 15.5 <= float(rs[i]) <16.25:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[2])[0],tempinfo(tempimages[2])[1],tempinfo(tempimages[2])[2],rvs[2],rverrs[2],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'M2'))
                #M3
                if 16.25 <= float(rs[i]) < 17:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[6])[0],tempinfo(tempimages[6])[1],tempinfo(tempimages[6])[2],rvs[6],rverrs[6],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'M3'))
                #M3.5
                if 17<= float(rs[i]) <18:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[7])[0],tempinfo(tempimages[7])[1],tempinfo(tempimages[7])[2],rvs[7],rverrs[7],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),images[j],'M3.5'))
#halpha.close()
text.close()
"""
        if 6560 < objwav[maxindex] <6570:
            if max(flux) >3.25*np.std(flux)+np.mean(flux):
                #K0
                if 11 <= float(rs[i]) < 12:
                    rv =find_rv.radial_velocity(objwav,flux,errs,halphatempinfo(tempimages[9])[0],halphatempinfo(tempimages[9])[1],halphatempinfo(tempimages[9])[2],rvs[9],rverrs[9],crosscorr_width,outpath+'Halpha_'+name+'_'+str(aps[i]))
                    halpha.write("%8s %4s %3s %4s %3s\n" % (ras[i],ids[i],aps[i],name,'K0'))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'K0'))
                #K3
                if 12 <= float(rs[i]) < 13:
                    rv =find_rv.radial_velocity(objwav,flux,errs,halphatempinfo(tempimages[4])[0],halphatempinfo(tempimages[4])[1],halphatempinfo(tempimages[4])[2],rvs[4],rverrs[4],crosscorr_width,outpath+'Halpha_'+name+'_'+str(aps[i]))
                    halpha.write("%8s %4s %3s %4s %3s\n" % (ras[i],ids[i],aps[i],name,'K3'))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'K3'))
                #K5
                if 13 <= float(rs[i]) < 13.5:
                    rv =find_rv.radial_velocity(objwav,flux,errs,halphatempinfo(tempimages[8])[0],halphatempinfo(tempimages[8])[1],halphatempinfo(tempimages[8])[2],rvs[8],rverrs[8],crosscorr_width,outpath+'Halpha_'+name+'_'+str(aps[i]))
                    halpha.write("%8s %4s %3s %4s %3s\n" % (ras[i],ids[i],aps[i],name,'K5'))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'K5'))
                #K7
                if 13.5 <= float(rs[i]) < 14.25:
                    rv =find_rv.radial_velocity(objwav,flux,errs,halphatempinfo(tempimages[0])[0],halphatempinfo(tempimages[0])[1],halphatempinfo(tempimages[0])[2],rvs[0],rverrs[0],crosscorr_width,outpath+'Halpha_'+name+'_'+str(aps[i]))
                    halpha.write("%8s %4s %3s %4s %3s\n" % (ras[i],ids[i],aps[i],name,'K7'))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'K7'))
                #M0
                if 14.25 <= float(rs[i]) <15:
                    rv =find_rv.radial_velocity(objwav,flux,errs,halphatempinfo(tempimages[1])[0],halphatempinfo(tempimages[1])[1],halphatempinfo(tempimages[1])[2],rvs[1],rverrs[1],crosscorr_width,outpath+'Halpha_'+name+'_'+str(aps[i]))
                    halpha.write("%8s %4s %3s %4s %3s\n" % (ras[i],ids[i],aps[i],name,'M0'))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'M0'))
                #M1
                if 15 <= float(rs[i]) < 15.5:
                    rv =find_rv.radial_velocity(objwav,flux,errs,halphatempinfo(tempimages[3])[0],halphatempinfo(tempimages[3])[1],halphatempinfo(tempimages[3])[2],rvs[3],rverrs[3],crosscorr_width,outpath+'Halpha_'+name+'_'+str(aps[i]))
                    halpha.write("%8s %4s %3s %4s %3s\n" % (ras[i],ids[i],aps[i],name,'M1'))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'M1'))
                #M2
                if 15.5 <= float(rs[i]) <16.25:
                    rv =find_rv.radial_velocity(objwav,flux,errs,halphatempinfo(tempimages[2])[0],halphatempinfo(tempimages[2])[1],halphatempinfo(tempimages[2])[2],rvs[2],rverrs[2],crosscorr_width,outpath+'Halpha_'+name+'_'+str(aps[i]))
                    halpha.write("%8s %4s %3s %4s %3s\n" % (ras[i],ids[i],aps[i],name,'M2'))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'M2'))
                #M3
                if 16.25 <= float(rs[i]) < 17:
                    rv =find_rv.radial_velocity(objwav,flux,errs,halphatempinfo(tempimages[6])[0],halphatempinfo(tempimages[6])[1],halphatempinfo(tempimages[6])[2],rvs[6],rverrs[6],crosscorr_width,outpath+'Halpha_'+name+'_'+str(aps[i]))
                    halpha.write("%8s %4s %3s %4s %3s\n" % (ras[i],ids[i],aps[i],name,'M3'))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'M3'))
                #M3.5
                if 17<= float(rs[i]) <18:
                    rv =find_rv.radial_velocity(objwav,flux,errs,halphatempinfo(tempimages[7])[0],halphatempinfo(tempimages[7])[1],halphatempinfo(tempimages[7])[2],rvs[7],rverrs[7],crosscorr_width,outpath+'Halpha_'+name+'_'+str(aps[i]))
                    halpha.write("%8s %4s %3s %4s %3s\n" % (ras[i],ids[i],aps[i],name,'M3.5'))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'M3.5'))

            else:
                #K0
                if 11 <= float(rs[i]) < 12:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[9])[0],tempinfo(tempimages[9])[1],tempinfo(tempimages[9])[2],rvs[9],rverrs[9],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'K0'))
                #K3
                if 12 <= float(rs[i]) < 13:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[4])[0],tempinfo(tempimages[4])[1],tempinfo(tempimages[4])[2],rvs[4],rverrs[4],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'K3'))
                #K5
                if 13 <= float(rs[i]) < 13.5:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[8])[0],tempinfo(tempimages[8])[1],tempinfo(tempimages[8])[2],rvs[8],rverrs[8],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'K5'))
                #K7
                if 13.5 <= float(rs[i]) < 14.25:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[0])[0],tempinfo(tempimages[0])[1],tempinfo(tempimages[0])[2],rvs[0],rverrs[0],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'K7'))
                #M0
                if 14.25 <= float(rs[i]) <15:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[1])[0],tempinfo(tempimages[1])[1],tempinfo(tempimages[1])[2],rvs[1],rverrs[1],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'M0'))
                #M1
                if 15 <= float(rs[i]) < 15.5:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[3])[0],tempinfo(tempimages[3])[1],tempinfo(tempimages[3])[2],rvs[3],rverrs[3],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'M1'))
                #M2
                if 15.5 <= float(rs[i]) <16.25:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[2])[0],tempinfo(tempimages[2])[1],tempinfo(tempimages[2])[2],rvs[2],rverrs[2],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'M2'))
                #M3
                if 16.25 <= float(rs[i]) < 17:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[6])[0],tempinfo(tempimages[6])[1],tempinfo(tempimages[6])[2],rvs[6],rverrs[6],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'M3'))
                #M3.5
                if 17<= float(rs[i]) <18:
                    rv =find_rv.radial_velocity(objwav,flux,errs,tempinfo(tempimages[7])[0],tempinfo(tempimages[7])[1],tempinfo(tempimages[7])[2],rvs[7],rverrs[7],crosscorr_width,outpath+name+'_'+str(aps[i]))
                    text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %2f %4s %3s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv[0],ras[i],ids[i],aps[i],float(rv[1]),float(rv[2]),name,'M3.5'))
"""
