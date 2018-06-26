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
#where the images and raw data are
inpath = '/Volumes/MADDIE/'+obsrun+'/images/'
#where you want the legend (which spectrum is run with which template) saved
legendpath = '/Volumes/MADDIE/'+obsrun+'/final_stuff/'
#where you have the files.txt files that have which image corresponds with which raw image
filepath = '/Volumes/MADDIE/'+obsrun
#where you have the templates(the have the same start so I included the start here too)
temppath = '/Volumes/MADDIE/Forfxcor/new1_'
#where you want the fxcor outputs saved(all fxcor outputs for each data run should be in the same place)
outpath = '/Volumes/MADDIE/'+obsrun+'/final_stuff/donefxcor/'
#the templates to correlate against
tempimages = ['1514217i','1629119i','1647952i','1649165i','1650473i','1688171i','1741422i','1756192p','1787492i','1789183i']
#the rvs of the templates (make sure this is in the right order)
rvs = [22.190,14.531,-71.084,-9.499,-26.417,-22.016,22.942,18.210,-23.118,-45.630]
#rv errs of the templates
rverrs = [.287,0.050,.125,.032,.086,.108,.055,.030,.162,.056]
#feb2017 region of spectrum to correlate if you want to avoid halpha
#region = "6425-6550,6570-6800"
#normal region of spectrum to correlate if you want to avoid halpha
region = "6300-6550,6570-6800"

c=2.99792*10**5
halpha =6562.801

def wav(image):
    return np.linspace(image[0].header['CRVAL1'],image[0].header['CRVAL1']+(image[0].header['NAXIS1']-1)*image[0].header['CDELT1'],image[0].header['NAXIS1'])
#keeps track of which spectrum was run with which template
text = open(legendpath+'/legend.txt','w')
#keeps track of which spectrum did not correlate the halpha region
halphatext = open(legendpath+'/noHalphalgend.txt','w')

images = np.loadtxt(filepath+'/files.txt', usecols = (1,),dtype = str)
origimages = np.loadtxt(filepath+'/files.txt', usecols = (0,),dtype = str)

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
            # have this if statement because 2 of the spectra would get strange errors if we tried to not correlate the h-alpha region
            #This is where i decide which template it should be ran against based on the r magnitude,
            #right now it is set up to run so that every spectra does not correlate the halpha region except for the 2 that it doesn't work for
            #there is also part of this code that runs everything so that it correlates the whole spectrum unless an halpha emission feature is detected
            #That part is commented out at the bottom of this script, so you can change it back to that if you want
            #im guessing you'll need to make a lot of changes here, for each fxcor command it goes:
            #iraf.rv.fxcor(input image, template image, apertures = str(int(aps[index])), output = outpath+images[j]+'_'+str(writeras[index]), verbose = 'txtonly', osample = '*' or region, rsample = '*' or region)
            #"*' mean correlate the whole spectrum, region is the region on either side of halpha
            # I think it's pretty important that you keep the outputs file names formatted in the way that I have them to be able to use the later scripts I have written
            if (images[j] == 'PF2' and obsrun == 'Dec2016' and str(int(aps[index]))  == '17') or (images[j] == 'PF1' and obsrun == 'Feb2017' and str(int(aps[index]))  == '9'):
                        #K0
                        if 11 <= float(rs[index]) < 12:
                            #iraf.rv.fxcor(inpath+images[j]+'_nosky.fits','/Volumes/MADDIE/Forfxcor/new1_1789183i.fits',apertures = str(int(aps[index])), output = outpath+images[j]+"_"+str(writeras[index]), verbose='txtonly',osample = "6300-6550,6570-6800", rsample = "6300-6550,6570-6800")
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
            #detects if there is halpha emission then only correlates the regions around halpha
            if 6560 < objwav[maxindex] <6570:
                if max(flux) >3.25*np.std(flux)+np.mean(flux):
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
            #correlates the entire spectrum
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
