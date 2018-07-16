from astropy.io import fits
import numpy as np
#import matplotlib.pyplot as plt
import os
from astropy.io.fits import getdata
from astropy.io.fits import getheader
from astropy.wcs import WCS
import scipy.ndimage
from pysynphot import observation
from pysynphot import spectrum
from pyraf import iraf

#run in pyraf window
print('loading RV info')
#list = np.loadtxt('RV_full_info.txt', dtype='string', skiprows=1, usecols=0, unpack=True)
list = np.loadtxt('RV_full_ESO.txt', dtype='string', skiprows=1, usecols=0, unpack=True)


#list of template files you want to edit, assumes these files have wav column first (in angstroms) and flux column second. 
#list = ['1650473i','1514217i','1629119i','1647952i','1649165i','1741422i','1787492i','1756192p','1688171i','1789183i']
#where you have the orig ascii filaes saved
inpath = '/Users/Natalie/mann/templates/'
#what you want to new begining of the new ascii and fits files to be (or you can go down and totally change the file name)
outname = 'new_'
#where you want the new fits file saved
outpath = inpath


def rebin_spec(wave, specin, wavnew):
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew, force='taper')
    return obs.binflux

#final wavelength dispersion you want, the wavelgenth of one of the data files
objwav = np.linspace(6339.68457,6339.68457+.21082*2031,2032)

for i in list:
    goodwav = []
    goodflux = []
    wav4 = np.loadtxt(inpath+i+'.ascii', usecols= (0), skiprows=1)
    flux4 = np.loadtxt(inpath+i+'.ascii', usecols = (1), skiprows=1)
    for j in range(len(wav4)):
        if 6339.68457  <= float(wav4[j]) <= 6767.86:
            goodwav.append(wav4[j])
            goodflux.append(flux4[j])
    try:
        print('running '+inpath+i+'.ascii...')
        goodwav = np.array(goodwav)
        goodflux = np.array(goodflux)
        newdata = rebin_spec(goodwav,goodflux,objwav)
        file= open(inpath+outname+i+'.ascii',"w")
        for h in range(len(newdata)):
            file.write('%10.5f %10.5f\n' % (objwav[h],newdata[h]))
        file.close()
          #turn ascii to fits
        iraf.rspectext(inpath+outname+i+'.ascii', outpath+outname+i+'.fits', dtype= 'interp')
    except:
        print(inpath+i+'.ascii'+' failed')
    #newdata = rebin_spec(goodwav,goodflux,objwav)
    #file= open(inpath+outname+i+'.ascii',"w")
    
