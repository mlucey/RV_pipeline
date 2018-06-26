from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join
import itertools
import pandas
from scipy import optimize as op
import scipy.integrate as integrate
import math

from pyraf import iraf
iraf.rv()



obsrun = 'Apr2018'

c=2.99792*10**5
halpha =6564.614
height = 200
inpath = '/Volumes/MADDIE/'+obsrun+'/images/'

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def dopplershift(shift):
    return c*((shift/halpha)**2-1)/(1+(shift/halpha)**2)

def dopplershifterror(shift,error):
    return c*(4*shift*halpha**2)/((halpha**2+shift**2)**2)*error

#change wav for Feb 2017 data
#feb2017 wav = np.linspace(6402.666,6402.666+2031*.2072,2031)
# feb132016wav = np.linspace(6339.684,6339.684+2032*.2108,2032)
def wav(image):
    return np.linspace(image[0].header['CRVAL1'],image[0].header['CRVAL1']+(image[0].header['NAXIS1']-1)*image[0].header['CDELT1'],image[0].header['NAXIS1'])
# March2016 wav=np.linspace(6339.670,6339.670+2029*.2109,2030)

halphafiles = []
halphaaps = []
halpharas = []


hght = np.loadtxt('/Volumes/MADDIE/'+obsrun+'/final_stuff/rvcorrect.txt',usecols = (12,),dtype = str)
files = np.loadtxt('/Volumes/MADDIE/'+obsrun+'/final_stuff/rvcorrect.txt',usecols = (13,),dtype = str)
aps = np.loadtxt('/Volumes/MADDIE/'+obsrun+'/final_stuff/rvcorrect.txt',usecols = (10,))
ra = np.loadtxt('/Volumes/MADDIE/'+obsrun+'/final_stuff/rvcorrect.txt',usecols = (9,),dtype = str)

for i in range(len(hght)):
    if hght[i] == 'INDEF' or float(hght[i]) < .5:
        halphafiles.append(files[i])
        halphaaps.append(aps[i])
        halpharas.append(ra[i])

file = open('/VOlumes/MADDIE/'+obsrun+'/final_stuff/Halphalegend.txt','w')
for k in range(len(halphafiles)):
    iraf.rv.fxcor(inpath+halphafiles[k]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/1647952iplushalpha.fits',apertures = str(int(halphaaps[k])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/halpha_"+halphafiles[k]+"_"+str(halpharas[k]),verbose='txtonly')
    file.write("%8s %8s %8.3f \n" % (halpharas[k],'1647952i',-71.084))
file.close()
"""
sigmas = []
cs= []
ras = []


for i in halpharasi:
    ras1 = i.split(":")
    ra = ras1[0]+ras1[1]+ras1[2]
    ras.append(ra)

"""
"""
fitfiles = ['PF1','PlF1','PlF1','PlF2','PlF2']
fitnums = [5,16,9,7,10]

for i in range(len(fitfiles)):
    pic = fits.open('/Volumes/Maddie/'+obsrun+'/images/'+fitfiles[i]+'_nosky.fits')
    val=pic[0]
    flux=val.data[fitnums[i]][0:2030]
    plt.plot(wav,flux)
    plt.show()
    height = max(flux)
    opt,cov = op.curve_fit(gauss_function,wav,flux,p0=(height,halpha,1,500))
    sigmas.append(opt[2])
    cs.append(opt[3])



tempflux = gauss_function(oldtempwav,height,halpha,np.mean(sigmas),np.mean(cs))

shiftedwav = dopplershift(v,oldtempwav)


image = fits.open('/Volumes/MADDIE/Forfxcor/new1_1647952i.fits')
values = image[0]
flux = values.data



file = open(os.path.join('/Volumes/MADDIE/Forfxcor/Halphaplustemp.ascii'),'w')
for i in range(len(tempflux)):
    file.write('%10.5f %10.5f \n' % (oldtempwav[i],(tempflux[i]/height)+flux[i]))
file.close()

iraf.rspectext('/Volumes/MADDIE/Forfxcor/halphaplustemp.ascii','/Volumes/MADDIE/Forfxcor/halphaplustemp.fits', dtype = 'interp')
"""
"""

text = open('/Volumes/MADDIE/Halpha/'+obsrun+'/final_stuff/gaussrvcorrect.txt','w+')

for i in range(len(Halphafiles)):
    #iraf.scopy('/Volumes/MADDIE/'+obsrun+'/images/'+Halphafiles[i]+'.ms.fits',
    #'/Volumes/MADDIE/'+obsrun+'/images/'+Halphafiles[i], apertures = halphaaps[i],
    #format = 'onedspec')
    if len(str(halphaaps[i])) < 2:
        filepath = '/Volumes/MADDIE/'+obsrun+'/images/'+Halphafiles[i]+'.000'+str(halphaaps[i])+'.fits'
        file = fits.open(filepath)
        data = file[0].data
        date = file[0].header['DATE-OBS']
        date1 = date.split("T")
        dates = date1[0].split("-")
        ut = file[0].header['UT']
        ra = file[0].header['RA']
        dec = file[0].header['DEC']

        fits.setval(filepath, "WAT1_001", value = 'wtype=linear label=Wavelength units=Angstrom')
        sp = pyspeckit.Spectrum('/Volumes/MADDIE/'+obsrun+'/images/'+Halphafiles[i]+".000"+str(halphaaps[i])+'.fits')
        start = (np.abs(wav(file)-6555)).argmin()
        stop = (np.abs(wav(file)-6561)).argmin()
        start1 = (np.abs(wav(file)-6565)).argmin()
        stop1 = (np.abs(wav(file)-6570)).argmin()
        foravg = data[start:stop].tolist()+data[start1:stop1].tolist()
        avg = np.mean(foravg)
        avglist = np.linspace(avg,avg,len(wav(file)))
        sp_avg = pyspeckit.Spectrum(data=avglist, xarr= wav(file))
        sp_contsub = sp.copy()
        sp_contsub.data -= sp_avg.data
        sp_contsub.specfit(fittype = 'gaussian',guesses = [height, halpha, 1],limits=[(0,0), (6555,6570), (0,3)],
        limited=[(True,False), (True,True), (True,True)])
        plt.plot(wav(file),sp_contsub.data)
        plt.plot(wav(file),gauss_function(wav(file),*(sp_contsub.specfit.parinfo)))
        plt.xlim(6555,6570)
        #plt.show()
        #print(sp_contsub.specfit.parinfo)
        #print(dopplershift(sp_contsub.specfit.parinfo.values[1]))
        rv= dopplershift(sp_contsub.specfit.parinfo.values[1])
        err = dopplershifterror(sp_contsub.specfit.parinfo.values[1],sp_contsub.specfit.parinfo[1].error)
        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %4s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv,ras[i],halphaids[i],halphaaps[i],float(err),Halphafiles[i]))
    else:
        filepath = '/Volumes/MADDIE/'+obsrun+'/images/'+Halphafiles[i]+'.00'+str(halphaaps[i])+'.fits'
        file = fits.open(filepath)
        data = file[0].data
        date = file[0].header['DATE-OBS']
        date1 = date.split("T")
        dates = date1[0].split("-")
        ut = file[0].header['UT']
        ra = file[0].header['RA']
        dec = file[0].header['DEC']

        fits.setval(filepath, "WAT1_001", value = 'wtype=linear label=Wavelength units=Angstrom')
        sp = pyspeckit.Spectrum('/Volumes/MADDIE/'+obsrun+'/images/'+Halphafiles[i]+".00"+str(halphaaps[i])+'.fits')
        start = (np.abs(wav(file)-6555)).argmin()
        stop = (np.abs(wav(file)-6561)).argmin()
        start1 = (np.abs(wav(file)-6565)).argmin()
        stop1 = (np.abs(wav(file)-6570)).argmin()
        foravg = data[start:stop].tolist()+data[start1:stop1].tolist()
        avg = np.mean(foravg)
        avglist = np.linspace(avg,avg,len(wav(file)))
        sp_avg = pyspeckit.Spectrum(data=avglist, xarr= wav(file))
        sp_contsub = sp.copy()
        sp_contsub.data -= sp_avg.data
        sp_contsub.specfit(fittype = 'gaussian',guesses = [height, halpha, 1],limits=[(0,0), (6555,6570), (0,3)],
        limited=[(True,False), (True,True), (True,True)])
        #plt.plot(wav(file),sp_contsub.data)
        #plt.plot(wav(file),gauss_function(wav(file),*(sp_contsub.specfit.parinfo)))
        #plt.xlim(6555,6570)
        #plt.show()
        #print(dopplershift(sp_contsub.specfit.parinfo.values[1]))
        rv= dopplershift(sp_contsub.specfit.parinfo.values[1])
        err = dopplershifterror(sp_contsub.specfit.parinfo.values[1],sp_contsub.specfit.parinfo[1].error)
        text.write("%4s %2s %2s %8s %8s %8s %1s %6.3f %8s %4s %3s %6.3f %4s\n" % (dates[0],dates[1],dates[2],ut,ra,dec,0,rv,ras[i],halphaids[i],halphaaps[i],float(err),Halphafiles[i]))
text.close()
    #feb132016 iraf.rv.fxcor('/Volumes/MADDIE/'+obsrun+'/'+Halphafiles[i]+'_nosky.fits',
    #'/Volumes/MADDIE/Forfxcor/halphatemp.fits',apertures = str(int(halphaaps[i])),
    #output = "/Volumes/MADDIE/Halpha/donefxcor/"+obsrun+'/'+str(Halphafiles[i])+"_"+str(ras[i])+'_'+str(halphaids[i]),
    #verbose='txtonly')

text1 = open('/Volumes/MADDIE/Forfxcor/halphafile.txt','w+')

M2s = []
for k in range(len(Halphafiles)):
    if 15.5 <= float(halphars[k]) <16.25:
        M2s.append(k)
print(ras[M2s[0]])
"""
"""
filepath = '/Volumes/MADDIE/'+obsrun+'/images/'+Halphafiles[M2s[0]]+'.000'+str(halphaaps[M2s[0]])+'.fits'
file = fits.open(filepath)
data = file[0].data
for i in range(len(data)):
    text1.write('%10.5f %10.5f \n' % (wav(file)[i],data[i]))


text1.close()
"""
