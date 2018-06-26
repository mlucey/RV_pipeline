from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io.fits import getdata
from astropy.io.fits import getheader
from astropy.wcs import WCS
import scipy.ndimage
from pysynphot import observation
from pysynphot import spectrum


"""
'1650473i'
'1288857i','1289045i'
"""
list = ['1650473i','1514217i','1629119i','1647952i','1649165i','1741422i','1787492i','1756192p','1688171i','1789183i']
list1 = ['sun']
"""
#BAD
for k in list:
    temp = fits.open('/Volumes/MADDIE/RV_std/'+k+'.fits')
    wavestart = (temp[0].data[0][0])*10
    dwave = (temp[0].data[0][1]-temp[0].data[0][0])*10
    outfile = '/Volumes/MADDIE/RV_std/new1_'+str(k)+'.fits'
    head= temp[0].header
    flux = temp[0].data[1]
    fits.writeto(outfile,flux,header=head)

    fits.setval(outfile, 'CDELT1', value = dwave)
    fits.setval(outfile, 'CRVAL1', value= wavestart)
    fits.setval(outfile, 'CTYPE1', value = 'LINEAR')
    fits.setval(outfile, 'CRPIX1' , value = 1.)
    fits.setval(outfile, 'CD1_1', value = dwave)

    w=WCS('/Volumes/MADDIE/RV_std/'+k+'.fits')

    w.all_pix2world(temp[0].data[0], 2)

    temp.close()
"""
"""
#ALSO BAD
for j in list:
    temp =fits.open('/Volumes/MADDIE/RV_std/'+j+'.fits')
    val=temp[0]
    wav=val.data[0]
    flux=val.data[1]
    file= open(os.path.join('/Volumes/MADDIE/RV_std',j+".ascii"),"w")
    newwavj =[]
    newfluxj = []

    for k in range(len(wav)):
        if wav[k]> 633.8:
            if wav[k]<676.8:
                newwavj.append(wav[k])
                newfluxj.append(flux[k])

    newerwav = scipy.ndimage.interpolation.zoom(newwavj,.13)
    newerflux = scipy.ndimage.interpolation.zoom(newfluxj,.13)
    for h in range(len(newerwav)):
        file.write('%10.5f %10.5f\n' % ( 10*newerwav[h],newerflux[h]))
    file.close()
    temp.close()

"""

#create good ascii with good rebinned data

def rebin_spec(wave, specin, wavnew):
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew, force='taper')

    return obs.binflux



"""
for j in list1:

    repeats = []
    fluxj2 = []
    wavj2 = []
    newwavj =[]
    newfluxj = []
    temp =fits.open('/Volumes/MADDIE/RV_std/'+j+'.fits')
    val=temp[0]
    wav=val.data[0]
    flux=val.data[1]
    file= open(os.path.join('/Volumes/MADDIE/RV_std',j+".ascii"),"w")
    for k in range(len(wav)):
        if wav[k]> 633.8:
            if wav[k]<676.8:
                newwavj.append(10*wav[k])
                newfluxj.append(flux[k])
    temptuple=zip(newwavj,newfluxj)
    sortbywave = sorted(temptuple,key = lambda tup:tup[0])
    wav1,flux1 =zip(*sortbywave)
    for l in range(len(wav1)-1):
        if wav1[l]== wav1[l+1]:
            repeats.append(wav1[l])
        else:
            fluxj2.append(flux1[l])
            wavj2.append(wav1[l])

    temptuple=zip(wavj2,fluxj2)
    sortbywave = sorted(temptuple,key = lambda tup:tup[0])
    wav3,flux3 =zip(*sortbywave)
"""
objwav = np.linspace(6339.68457,6339.68457+.21082*2031,2032)


wav4 = np.loadtxt('/Volumes/MADDIE/RV_std/ogsun.ascii', usecols= (0,))
flux4 = np.loadtxt('/Volumes/MADDIE/RV_std/ogsun.ascii', usecols = (1,))

newdataj = rebin_spec(wav4*10,flux4,objwav)
file= open(os.path.join('/Volumes/MADDIE/RV_std',"sun.ascii"),"w")
for h in range(len(newdataj)):
    file.write('%10.5f %10.5f\n' % (objwav[h],newdataj[h]))
file.close()

"""
#turn ascii to fits
import sys
from pyraf import iraf
iraf.noao.onedspec()

for h in list1:
    iraf.rspectext('/Volumes/MADDIE/RV_std/'+h+'.ascii', '/Volumes/MADDIE/Forfxcor/new1_'+h+'.fits', dtype= 'interp')
"""
