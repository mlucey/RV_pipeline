#!/usr/bin/env python
# encoding: utf-8
"""
find_rv.py
Created by Vivienne Baldassare on 2012-03-30.
Adapated by Natalie Gosnell on 2018-05-22
and Maddie Lucey on 2018-06-04

Description:
    This code finds the radial velocity of a target when supplied with data for the target and data for a standard object
    whose radial velocity is known.
Usage:
    Note: Data is not corrected for heliocentric velocities!

    Inputs:
        wv_obj, fx_obj, and sig_obj are arrays containing data for the the wavelength, flux, and flux uncertainty of the target.
        wv_std, fx_std, and sig_std are arrays containing data for the the wavelength, flux, and flux uncertainty of the standard.
        rv_std is the radial velocity of the standard star
        rv_std_err is the error of the standard star radial velocity
        crosscorr_width is the pixel shift range used to calculate cross correlations
        figurename is a string for the final figure, 'figurename'.pdf
    Example:
        >>> import find_rv
        >>> find_rv.radial_velocity(wv_obj,fx_obj,sig_obj,wv_std,fx_std,sig_std,rv_std,rv_std_err,crosscorr_width,figurename)
"""
from astropy.io import fits
import os
from os import listdir
from os.path import isfile, join
from array import array
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import scipy
import scipy.optimize as op
import scipy.ndimage

def quadratic(x,a,b,c):
    return a + b*x + c*x**2

def cubic(x,a,b,c,d):
    return a + b*x + c*x**2 + d*x**3

def quartic(x,a,b,c,d,e):
    return a + b*x + c*x**2 + d*x**3 + e*x**4


def lsf_rotate(deltav,vsini,epsilon=None,velgrid=None):
    # based on the IDL routine LSF_ROTATE.PRO
    if epsilon == None:
        epsilon = 0.6

    e1 = 2.0*(1.0 - epsilon)
    e2 = np.pi*epsilon/2.0
    e3 = np.pi*(1.0 - epsilon/3.0)

    npts = np.ceil(2*vsini/deltav)
    if npts % 2 == 0:
        npts += 1
    nwid = np.floor(npts/2)
    x = np.arange(npts) - nwid
    x = x*deltav/vsini
    x1 = np.abs(1.0 - x**2)
    if velgrid != None:
        velgrid = x*vsini
        return (e1*np.sqrt(x1) + e2*x1)/e3,velgrid

    return (e1*np.sqrt(x1) + e2*x1)/e3

def radial_velocity(wv_obj,fx_obj,sig_obj,wv_std,fx_std,sig_std,rv_std,rv_std_err,crosscorr_width,figurename):


    # The more random iterations, the better... but it takes longer
    n_iter = 1000

    # Step 1: Fix the spectra:
    # * Select only the region in which they overlap
    # * Make a new stretched wavelength array (for sub-pixel precision work)
    # * Interpolate the data onto the new wavelength array
    # * Remove large scale slopes so we only compare line and band features

    # Find where standard and object overlap ---------------
    wv_min = max([min(wv_std),min(wv_obj)])
    wv_max = min([max(wv_std),max(wv_obj)])

    n_pix_std = len(wv_std)

    # Creates ln standard wavelength array ---------------------------------
    # AR 2013.0423 The wavelength array only covers the overlap region.  Also, I'm folding the rebinning by 10 into this statement.
    acoef_std = (n_pix_std*10 -1)/(math.log(wv_max) - math.log(wv_min))
    bcoef_std = (n_pix_std*10) - (acoef_std * math.log(wv_max))

    arr = np.arange(n_pix_std*10)+1
    wv_ln_std = np.exp((arr - bcoef_std)/acoef_std)

    # AR 2012.1018: Find the conversion between pixels and velocity.  This will vary from instrument
    #  to instrument and spectral order to spectral order, so we should preferentially calculate this
    #  based on the actual input spectrum.
    # AR 2013.0422: Change the calculation to happen AFTER the corrected wavelength scale has been made
    # Find the average pixel/spectrum offset
    #  Note: even though it's called micron_per_pix, it will still work if the wavelengths are
    #  angstroms instead (it really converts <wavelength unit> to km/s)

    # Interpolate data onto same ln wavelength scale -------------------------------

    fx_interp_std = np.interp(wv_ln_std, wv_std, fx_std)
    fx_interp_obj = np.interp(wv_ln_std, wv_obj, fx_obj)
    sig_interp_std = np.interp(wv_ln_std, wv_std, sig_std) # AR 2012.1018 Also need to rebin sig
    sig_interp_obj = np.interp(wv_ln_std, wv_obj, sig_obj) # AR 2012.1018 Also need to rebin sig

    # Rebin Data ----------------------------

    wv_arr_std=np.asarray(wv_ln_std,dtype=float)
    fx_arr_obj=np.asarray(fx_interp_obj,dtype=float)
    fx_arr_std=np.asarray(fx_interp_std,dtype=float)
    sig_arr_obj=np.asarray(sig_interp_obj,dtype=float)
    sig_arr_std=np.asarray(sig_interp_std,dtype=float)

    datalen = len(fx_arr_obj)

    pix_scale = (2.99792458*10**5)/acoef_std

    # 3. Cross-correlate the data, using n_iter trials:
    #  * Generate two random gaussian noises scaled to the uncertainty on the fluxes
    #  * Apply those gaussian noises to the standard and target stars
    #  * Cross-correlate the standard and target stars
    #  * Find and then cut out just the part of the cross-correlation curve near the maximum
    #  * Set up gaussian
    #  * Fit gaussian to that center part
    #  * Save fitted parameters (pixel shift aka mean of gaussian, width aka stddev of gaussian)
    #  * Repeat n_iter times

    # Cross correlation loop --------------------------------
    pix_shift=np.zeros(n_iter)    #initialize array for pixel shift values
    pix_width=np.zeros(n_iter)    #initialize array for pixel width values
    ccf_peak=np.zeros(n_iter)   #initialize array for ccf peak values
    l = 0


    # using the xrange generator rather than making a full list saves memory
    for l in xrange(n_iter):
        # prepare the randomized data
        # GETTING ARRAYS READY FOR CROSS CORRELATION


        # Randomize noise:
        # create gaussian distribution of random numbers b/t 1 and -1, multiply err by numbers, add numbers to flux
        # I have drastically simplified the arrays here AR 2013.0319
        # AR 2013.0318: There was a problem, previously: noise was a fixed value, not linked to the known error values

        # AR 2013.0321: Speed fix.  Rather than step through the array and generate one
        #  normally-distributed error value scaled to the SNR at that point, I will generate an
        #  array of normally-distributed error values scaled to 1, and then multiply by the SNR:
        #  One array generation, one array multiplication.

        rand_dist = np.random.normal(loc=0.0,scale=1.0,size=datalen)
        rand_dist2 = np.random.normal(loc=0.0,scale=1.0,size=datalen)

        fx_temp_obj = np.asarray(fx_arr_obj + rand_dist * sig_arr_obj)
        fx_temp_std = np.asarray(fx_arr_std + rand_dist2 * sig_arr_std)
        mean_obj=np.mean(fx_temp_obj)
        mean_std=np.mean(fx_temp_std)
        stddev_obj=np.std(fx_temp_obj,ddof=1)
        stddev_std=np.std(fx_temp_std,ddof=1)

        # Regularize data (subtract mean, divide by std dev) (Should definitely be done AFTER noise was added)
        fx_reg_temp_obj = fx_temp_obj-mean_obj
        fx_reg_temp_obj = fx_reg_temp_obj/stddev_obj
        fx_reg_temp_std = fx_temp_std-mean_std
        fx_reg_temp_std = fx_reg_temp_std/stddev_std

        # curve fit - remove a cubic AR 2012.1113
        coeff,pcov = op.curve_fit(cubic,wv_arr_std,fx_reg_temp_obj)
        fx_reg_temp_obj = fx_reg_temp_obj - (coeff[0] + coeff[1]*wv_arr_std + coeff[2]*wv_arr_std**2 + coeff[3]*wv_arr_std**3)
        coeff,pcov = op.curve_fit(cubic,wv_arr_std,fx_reg_temp_std)
        fx_reg_temp_std = fx_reg_temp_std - (coeff[0] + coeff[1]*wv_arr_std + coeff[2]*wv_arr_std**2 + coeff[3]*wv_arr_std**3)

        # CROSS CORRELATION

        #cross-correlation of template with itself
        ystd = np.correlate(fx_reg_temp_std,fx_reg_temp_std)
        # compute the cross-correlation between the two spectra
        ycorr = np.correlate(fx_reg_temp_obj, fx_reg_temp_std, mode='full')
        #Cross-correlation spectra with itself
        ynormal = np.correlate(fx_reg_temp_obj,fx_reg_temp_obj)
        # time required: 0.045 seconds average

        # create the x offset axis (same length as ycorr, with 0 in the MIDDLE)
        length = len(ycorr)
        xcorr = np.arange(length) - length//2
        # AR 2012.1126 Select a tiny piece around the maximum to fit with a gaussian.
        xmid = np.argmax(ycorr)
        ymax = np.max(ycorr)
        # now take just the portion of the array that matters
        xcorr_min=int(xmid-crosscorr_width)
        xcorr_max=int(xmid+crosscorr_width)
        ycorr1=ycorr[xcorr_min:xcorr_max]    #isolate section of array with gaussian
        xcorr1=xcorr[xcorr_min:xcorr_max]       #isolate the same section of the pixel range
        ycorr2=ycorr[xcorr_min-50:xcorr_max+50]
        xcorr2=xcorr[xcorr_min-50:xcorr_max+50]

        # suggestion from D. Hogg 12/15/12: Add extra linear feature to fit.
        # suggestion from D. Hogg 12/15/12: operate on ln(amp) so that the amplitude CANNOT be negative.
        def chi2(p):    #define gaussian function for fitting
            sig2=p[2] ** 2
            m = (np.exp(p[0]) * np.exp(-0.5 * (xcorr1 - p[1]) ** 2 / sig2)) + p[3] + p[4]*xcorr1
            return (ycorr1 - m)

        # set up initial values for chi2
        sig = 10
        sky = np.min(ycorr1)/1.2
        #                print ycorr1[-1],ycorr1[0],xcorr1[-1],xcorr1[0]
        sky2 = (ycorr1[-1]-ycorr1[0])/(xcorr1[-1]-xcorr1[0])
        lnamp = np.log(ymax/1.2-sky)    # guess some values
        mean = xcorr[xmid]

        amp = np.exp(lnamp)
        sig2 = sig**2

        popt, ier = op.leastsq(chi2, [lnamp, mean, sig, sky, sky2])
        lnamp, mean, sig, sky, sky2 = popt

        amp = np.exp(lnamp)

        print_num=l%np.floor(n_iter/2)        #prints data 10 times during calculation. Sort of a progress meter.
        if print_num == 0:
            ## Uncomment the following to make a plot every 500 fits.
            #fig = plt.figure(l)
            #ax = fig.add_subplot(111)
            #my_gauss = (amp * (np.exp(-0.5 * ((xcorr1 - mean) ** 2) / sig**2))) + sky + sky2 * xcorr1
            #ax.plot(xcorr1,my_gauss,'r--')
            #ax.plot(xcorr2,ycorr2,'#000000')
            #ax.plot(xcorr1,ycorr1-my_gauss,'#00CC00')
            ##if abs(mean - xcorr[xmid]) > 5:
            ##    print "Mean is off",mean,xcorr[xmid]
            #figname='rv_{0:}_{1:}_{2:}_{3:}.png'.format(std_name,obj_name,order,l)
            #ax.set_xlim(xcorr[xcorr_min-50],xcorr[xcorr_max+50])
                #fig.savefig(figname)
            #fig.clf()
            #plt.close()
            #print "amp={0: 12.4f}  mu={1: 10.4f}  sig={2: 9.4f}  sky={3: 11.4f}  sky2={4: 8.4f}".format(amp,mean,sig,sky,sky2)
            print("amp={0: 12.4f}  mu={1: 10.4f}  sig={2: 9.4f}".format(amp,mean,sig))


        # if ier < 5:
        pix_shift[l] = mean
        ccf_index = (np.abs(mean-xcorr)).argmin()
        ccf_peak[l] = ycorr[ccf_index]/ynormal
        # I'm calculating the vsini now because I need errors, and the vsini calculation is not linear.
        #pix_width[l] = vsinicoeff[0] + vsinicoeff[1] * sig + vsinicoeff[2] * sig**2 + vsinicoeff[3] * sig**3 + vsinicoeff[4] * sig**4



    # End cross correlation loop ---------------------------------

    # 4. Find the RV
    # All 5000 rv fits have been calculated and stored in arrays
    # 4a. Compute the mean pixel shift and pixel shift uncertainty.
    # 4b. Convert pixel shift into RV
    # 4c. Shift the wavelength array appropriately - all lines should now line up.

    ## Uncomment this to print out an example cross-correlation diagram
    #fig = plt.figure(2)
    #ax = fig.add_subplot(111)
    #ax.plot(xcorr,ycorr,'k')
    #figname='rv_{0:}_{1:}_{2:}_xcorr.png'.format(std_name,obj_name,order)
    #fig.savefig(figname)
    #fig.clf()
    #plt.close()



    # Turn the list of pixel shifts into a numpy array
    pix_shift = np.asarray(pix_shift)

    # 4a. Compute the mean pixel shift (rv value) and pixel shift uncertainty (RV uncertainty).

    mu = np.mean(pix_shift)
    sigma = np.std(pix_shift,ddof=1)


    ccfpeak = np.mean(ccf_peak)
    #vsini = np.mean(pix_width)
    #vsini_err = np.std(pix_width,ddof=1)

    #axh = figv.add_subplot(212)
    #n, bins, patches=axh.hist(pix_width,bins=30,normed=1.0,facecolor='green',align='mid')
    #figv.savefig('vsiniplot.png')
    #plt.clf()
    #plt.close()

    # 4b. Transform pixel shift to shift in radial velocity

    # AR 2013.0423: The actually appropriate method requires a speed-of-light correction. This works for both angstroms and microns.
    rv_meas = (2.99792458*10**5 * mu)/acoef_std
    rv_meas_err = (2.99792458*10**5 * sigma)/acoef_std

    rv_corr = rv_std + rv_meas
    rv_corr_err = (rv_std_err**2 + rv_meas_err**2)**(0.5)

    # 4c. Apply shift to arrays
    wv_rvcorr_obj=wv_arr_std * (1 - rv_meas/(2.99792458*10**5))

    ## 5. Create plots ---------------------------------
    # The plots are the only reason find_rv.py needs to know figurename

    # Plot object and standard so you can clearly see that shift exists --------------------------------
    fig = plt.figure(figsize = (10, 10))

    # AR 2013.0703 Regularize the spectra for display purposes in the final graph
    # I'm using the mean and stddev of the last random-added attempt so it won't be perfect...
    fx_reg_obj = fx_arr_obj-mean_obj
    fx_reg_obj = fx_reg_obj/stddev_obj
    fx_reg_std = fx_arr_std-mean_std
    fx_reg_std = fx_arr_std/stddev_std

    #Plots target and standard with shift applied
    ax1 = fig.add_subplot(311)
    ax1.plot(wv_rvcorr_obj, fx_reg_obj, 'red')
    ax1.plot(wv_arr_std, fx_reg_std, 'blue')
    ax1.set_xlabel('Wavelength (A)')
    ax1.set_ylabel('Normalized Flux')
    #to use the lines below, will need to add "target" string back into the radial_velocity function
    #target = 'Target: %s' %(obj_name)
    #ax1.annotate(target,xy=(.7,.9),xycoords='axes fraction',xytext=(.6,.9),textcoords='axes fraction',color='red')

    sig2=sig ** 2
    my_gauss = (amp * (np.exp(-0.5 * ((xcorr1 - mu) ** 2) / sig2))) + sky + sky2 * xcorr1

    #Plots example of gaussian fit to cross correlation function
    ax2 = fig.add_subplot(312)
    ax2.plot(xcorr1,  ycorr1, 'k.')
    ax2.plot(xcorr1, my_gauss, 'r--', linewidth=2)
    ax2.plot(xcorr1,ycorr1-my_gauss,'#00CC00')
    ax2.set_xlabel('Fit to sample cross correlation function')
    ax2.set_xlim(xcorr[xcorr_min-50],xcorr[xcorr_max+50])
    #print pix_shift


    ## Plot histogram of pixel shift values --------------------------------
    ax3 = fig.add_subplot(313)
    n, bins, patches=plt.hist(ccf_peak,bins=30,normed=1.0,facecolor='green',align='mid')
    ##Plot best fit gaussian over histogram
    y=mlab.normpdf(bins,mu,sigma)
    ax3.plot(bins,y,'r--',linewidth=2)
    ax3.set_xlabel('Pixel shift of target')
    ax3.set_ylabel('Frequency (Normalized)')
    rad='RV = %.3f +/- %.3f' %(rv_meas,rv_meas_err)
    corr = 'RV (corr) = %.3f +/- %.3f' %(rv_std + rv_meas, (rv_std_err**2 + rv_meas_err**2)**(0.5))
    pixfit = 'Pix shift = %5.2f +/- %5.2f' %(mu, sigma)
    pixdisp = 'Dispersion = %5.3f km/s/pix' %((2.99792458*10**5)/acoef_std)
    #vsinistr = 'VsinI = %.3f +/- %.3f' % (vsini,vsini_err)
    ax3.annotate(rad,xy=(.7,.9),xycoords='axes fraction',xytext=(.66,.9),textcoords='axes fraction',color='black')
    ax3.annotate(corr,xy=(.6,.8),xycoords='axes fraction',xytext=(.60,.8),textcoords='axes fraction',color='black')
    ax3.annotate(pixfit,xy=(.05,.9),xycoords='axes fraction',xytext=(.05,.9),textcoords='axes fraction',color='black')
    ax3.annotate(pixdisp,xy=(.05,.8),xycoords='axes fraction',xytext=(.05,.8),textcoords='axes fraction',color='black')
    fig.tight_layout()

    figname = '%s.pdf' %(figurename)
    fig.savefig(figname, dpi = 300)
    fig.clf()
    plt.close()

    #plt.figure(l+1)
    #plt.hist(pix_shift)

    #END RADIAL VELOCITY FUNCTION -----------------------------------------
    return rv_corr,rv_corr_err,ccfpeak
