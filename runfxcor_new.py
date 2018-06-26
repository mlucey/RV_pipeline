from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join, isdir
from pyraf import iraf
from astropy.time import Time

iraf.rv()
obsrun = 'Mar2018'
images = np.loadtxt('/Volumes/MADDIE/'+obsrun+'/files.txt', usecols = (1,),dtype = str)
origimages = np.loadtxt('/Volumes/MADDIE/'+obsrun+'/files.txt', usecols = (0,),dtype = str)
#feb2017
#region = "6425-6550,6570-6800"
#normal
region = "6300-6550,6570-6800"

"""
#Mar2016
origimages= ['March2016073','March2016081','March2016047','March2016053','March2016042','march2016025','march2016033']
#origfilestart =
images = ['PF1','PF2','PB1','PB2','PB3','PF4','PF5']
obsrun = 'Mar2016'
"""
"""
#Dec2016
origimages= ['dec16044','dec16080','dec16049']
images = ['PF1','PB1','PF2']
obsrun = 'Dec2016'
"""
"""
#Apr2018
#Apr2018

origimages = ['Apr18_03_023','Apr18_03_028','Apr18_04_005','Apr18_04_010']
images = ['PF2','PB1','PB2','PF3']
obsrun = 'Apr2018'
"""
inpath = '/Volumes/MADDIE/'+obsrun+'/images/'
outpath = '/Volumes/MADDIE/'+obsrun+'/final_stuff/'
temppath = '/Volumes/MADDIE/Forfxcor/'
c=2.99792*10**5
halpha =6562.801
crosscorr_width = 200
tempimages = ['1514217i','1629119i','1647952i','1649165i','1650473i','1688171i','1741422i','1756192p','1787492i','1789183i']
rvs = [22.190,14.531,-71.084,-9.499,-26.417,-22.016,22.942,18.210,-23.118,-45.630]
rverrs = [.287,0.050,.125,.032,.086,.108,.055,.030,.162,.056]

def wav(image):
    return np.linspace(image[0].header['CRVAL1'],image[0].header['CRVAL1']+(image[0].header['NAXIS1']-1)*image[0].header['CDELT1'],image[0].header['NAXIS1'])

text = open('/Volumes/MADDIE/'+obsrun+'/final_stuff/legend.txt','w')
halphatext = open('/Volumes/MADDIE/'+obsrun+'/final_stuff/noHalphalgend.txt','w')

for j in range(len(images)):
    #feb132016 origfile = fits.open('/Volumes/MADDIE/'+obsrun+'/WIYN/'+origfilestart+str(origimages[j])+'.fits')
    file = fits.open(inpath+images[j]+'_nosky.fits')
    origfile = fits.open(inpath+origimages[j]+'.fits')
    #origfile = fits.open(inpath+origimages[j]+'.fits')
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
        #pulling infor for rvcorrect
        jddate = file[0].header['JD']
        t = Time(float(jddate), format='jd', scale='utc')
        date = t.iso
        date1 = date.split(" ")
        dates = date1[0].split("-")
        ut = file[0].header['UT']
        ra = file[0].header['RA']
        dec = file[0].header['DEC']
        splits = (file[0].header['APID'+str(i+1)]).split(" ")
        if splits[0] in ids:
            index = ids.index(splits[0])
            flux = file[0].data[i][200:1800]
            errs = np.sqrt(abs(flux))
            imagewav = wav(file)
            objwav = imagewav[200:1800]
            listflux = flux.tolist()
            maxindex = listflux.index(max(flux))
            if (images[j] == 'PF2' and obsrun == 'Dec2016' and str(int(aps[index]))  == '17') or (images[j] == 'PF1' and obsrun == 'Feb2017' and str(int(aps[index]))  == '9'):
                        #K0
                        if 11 <= float(rs[index]) < 12:
                            #iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1789183i.fits',apertures = str(int(aps[index])), output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]), verbose='txtonly',osample = "6300-6550,6570-6800", rsample = "6300-6550,6570-6800")
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1789183i.fits',apertures = str(int(aps[index])), output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]), verbose='txtonly')
                            text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1789183i',rvs[9]))
                        #K3
                        if 12 <= float(rs[index]) < 13:
                            #iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1650473i.fits',apertures = str(int(aps[index])), output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]), verbose='txtonly',osample = "6300-6550,6570-6800", rsample = "6300-6550,6570-6800")
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1650473i.fits',apertures = str(int(aps[index])), output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]), verbose='txtonly')
                            text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1650473i',rvs[4]))
                        #K5
                        if 13 <= float(rs[index]) < 13.5:
                            #iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1787492i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly',osample = "6300-6550,6570-6800", rsample = "6300-6550,6570-6800")
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1787492i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly')
                            text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1787492i',rvs[8]))
                        #K7
                        if 13.5 <= float(rs[index]) < 14.25:
                            #iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1514217i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly',osample = "6300-6550,6570-6800", rsample = "6300-6550,6570-6800")
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1514217i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly')
                            text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1514217i',rvs[0]))
                        #M0
                        if 14.25 <= float(rs[index]) <15:
                            #iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1629119i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly',osample = "6300-6550,6570-6800", rsample = "6300-6550,6570-6800")
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1629119i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly')
                            text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1629119i',rvs[1]))
                        #M1
                        if 15 <= float(rs[index]) < 15.5:
                            #iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1649165i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly',osample = "6300-6550,6570-6800", rsample = "6300-6550,6570-6800")
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1649165i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly')
                            text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1649165i',rvs[3]))
                        #M2
                        if 15.5 <= float(rs[index]) <16.25:
                            #iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1647952i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly',osample = "6300-6550,6570-6800", rsample = "6300-6550,6570-6800")
                            iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1647952i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly')
                            text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1647952i',rvs[2]))
            else:
                #K0
                if 11 <= float(rs[index]) < 12:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1789183i.fits',apertures = str(int(aps[index])), output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]), verbose='txtonly' , osample = region, rsample = region)
                    halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1789183i',rvs[9]))
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1789183i',rvs[9]))
                #K3
                if 12 <= float(rs[index]) < 13:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1650473i.fits',apertures = str(int(aps[index])), output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]), verbose='txtonly', osample = region, rsample = region)
                    halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1650473i',rvs[4]))
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1650473i',rvs[4]))
                #K5
                if 13 <= float(rs[index]) < 13.5:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1787492i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly', osample = region, rsample = region)
                    halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1787492i',rvs[8]))
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1787492i',rvs[8]))
                #K7
                if 13.5 <= float(rs[index]) < 14.25:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1514217i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly', osample =region, rsample = region)
                    halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1514217i',rvs[0]))
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1514217i',rvs[0]))
                #M0
                if 14.25 <= float(rs[index]) <15:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1629119i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly', osample = region, rsample = region)
                    halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1629119i',rvs[1]))
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1629119i',rvs[1]))
                #M1
                if 15 <= float(rs[index]) < 15.5:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1649165i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly', osample = region, rsample = region)
                    halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1649165i',rvs[3]))
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1649165i',rvs[3]))
                #M2
                if 15.5 <= float(rs[index]) <18:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1647952i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly', osample = region, rsample = region)
                    halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1647952i',rvs[2]))
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1647952i',rvs[2]))
        else:
            print(splits[0])
            """
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
                    if 11 <= float(rs[index]) < 12:
                        iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1789183i.fits',apertures = str(int(aps[index])), output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]), verbose='txtonly' , osample = "6300-6550,6570-6800")
                        halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1789183i',rvs[9]))
                        text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1789183i',rvs[9]))
                    #K3
                    if 12 <= float(rs[index]) < 13:
                        iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1650473i.fits',apertures = str(int(aps[index])), output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(rwriteras[index]), verbose='txtonly', osample = "6300-6550,6570-6800")
                        halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1650473i',rvs[4]))
                        text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1650473i',rvs[4]))
                    #K5
                    if 13 <= float(rs[index]) < 13.5:
                        iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1787492i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly', osample = "6300-6550,6570-6800")
                        halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1787492i',rvs[8]))
                        text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1787492i',rvs[8]))
                    #K7
                    if 13.5 <= float(rs[index]) < 14.25:
                        iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1514217i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly', osample = "6300-6550,6570-6800")
                        halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1514217i',rvs[0]))
                        text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1514217i',rvs[0]))
                    #M0
                    if 14.25 <= float(rs[index]) <15:
                        iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1629119i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly', osample = "6300-6550,6570-6800")
                        halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1629119i',rvs[1]))
                        text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1629119i',rvs[1]))
                    #M1
                    if 15 <= float(rs[index]) < 15.5:
                        iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1649165i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly', osample = "6300-6550,6570-6800")
                        halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1649165i',rvs[3]))
                        text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1649165i',rvs[3]))
                    #M2
                    if 15.5 <= float(rs[index]) <18:
                        iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1647952i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly', osample = "6300-6550,6570-6800")
                        halphatext.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1647952i',rvs[2]))
                        text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1647952i',rvs[2]))
            """
            """
            else:
                #K0
                if 11 <= float(rs[index]) < 12:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1789183i.fits',apertures = str(int(aps[index])), output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]), verbose='txtonly',osample = "*")
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1789183i',rvs[9]))
                #K3
                if 12 <= float(rs[index]) < 13:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1650473i.fits',apertures = str(int(aps[index])), output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]), verbose='txtonly',osample = "*")
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1650473i',rvs[4]))
                #K5
                if 13 <= float(rs[index]) < 13.5:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1787492i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly',osample = "*")
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1787492i',rvs[8]))
                #K7
                if 13.5 <= float(rs[index]) < 14.25:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1514217i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly',osample = "*")
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1514217i',rvs[0]))
                #M0
                if 14.25 <= float(rs[index]) <15:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1629119i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly',osample = "*")
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1629119i',rvs[1]))
                #M1
                if 15 <= float(rs[index]) < 15.5:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1649165i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly',osample = "*")
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1649165i',rvs[3]))
                #M2
                if 15.5 <= float(rs[index]) <16.25:
                    iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1647952i.fits',apertures = str(int(aps[index])),output = "/Volumes/MADDIE/"+obsrun+"/final_stuff/donefxcor/"+images[j]+"_"+str(writeras[index]),verbose='txtonly',osample = "*")
                    text.write("%8s %4s %8s %8.3f \n" % (ids[index],ras[index],'1647952i',rvs[2]))
"""
text.close()
