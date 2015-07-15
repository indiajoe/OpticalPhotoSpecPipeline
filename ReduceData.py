#!/usr/bin/env python
#This script is to semi-automate basic data reduction of Optical Imaging and Spectroscopic data.
#Instruments it supports are
# HFOSC  @HCT,IAO,IIA,India
# IFOSC  @IGO,IUCAA,India

#------------------------------------Take a deep breath
# There was this programmer, the story goes, who was
# deeply mystified by a piece of code.  It had no comments
# at all, and he couldn't for the life of him figure out
# how it did what it did.  For years, he cursed the author
# of that code, but it continued to fascinate and trouble him.
# One day, it came to him in a flash.  He understood it all.
# In fact, it was so obvious that he also understood why
# it didn't need any comments!

# -- _Advanced Perl Programming_, Sriram Srinivasan

#----------------Relax, Don't PANIC, we will be coding in python..

# Advice for brave souls reading this code. Start at the bottom.
# And then slowly read functions one by one from bottom to top.
# I apologize for this upside down code, 
# I wish I could go back in time..
#---------------------------------------------indiajoe@gmail.com

import os
import os.path
import glob
#import pyfits as fits  #Deprecated and merged to astropy
import astropy.io.fits as fits
#import pyfits.convenience
import sys, traceback 
import numpy as np
import warnings
import re
import shlex
import readline
import shutil
import subprocess
import time
import datetime
from astropy.time import Time
from astropy.io import ascii
import astropy.table as table 

# Required for certain functions... Remove it if not needed for your module
from pyraf import iraf
iraf.set(stdimage="imt2048") #Setting the image size to 2048. Not important if people use ds9 7.1

try:
    from scipy.ndimage import filters
except ImportError:
    print('Scipy module missing.')
    print('You will need to install Scipy if and only if you need to do median filtered continuum division')

def SpectralExtraction_subrout(PC):
    """ Extracts spectra from 2d image and also applies wavelength calibration """
    iraf.noao(_doprint=0) 
    iraf.twodspec(_doprint=0) 
    iraf.onedspec(_doprint=0) 
    iraf.apextract(_doprint=0)
    iraf.apextract.unlearn()
    iraf.apall.unlearn()
    iraf.apsum.unlearn()
    iraf.reidentify.unlearn()
    iraf.dispcor.unlearn()
    iraf.scombine.unlearn()

    iraf.apextract.setParam('dispaxis',PC.DISPAXIS) # 2 is vertical and 1 is horizontal
    iraf.cd(os.path.join(PC.MOTHERDIR,PC.OUTDIR))
    directories = LoadDirectories(PC,CONF=False)
    for night in directories:
        PC.currentnight = night # Upgating the night directory for using PC.GetFullPath()
        print('Working on night: '+night)
        iraf.cd(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night))        
        try:
            #First we load the Spectrastoextract_Lamp_BandFilter.txt
            SpecslistFILE=open('Spectrastoextract_Lamp_BandFilter.txt','r')
        except IOError:
            print('Could not fine Spectrastoextract_Lamp_BandFilter.txt in this directory. Hence skipping %s directory'%(night))
            continue
        #Image to Lamp file dictionary
        Img2Lamp = []  #They are first born as lists...
        Img2Filt = []
        AllImglist = []
        for line in SpecslistFILE:
            line = line.rstrip().split()
            Img2Lamp.append((line[0],line[1]))
            Img2Filt.append((line[0],line[2]))
            AllImglist.append(line[0])
        SpecslistFILE.close()
        Img2Lamp = dict(Img2Lamp)
        Img2Filt = dict(Img2Filt)
        #Initialising an empty dictionary for storing spectra of each filter
        Filt2finalspecs = dict()
        for filt in set(Img2Filt.values()) : Filt2finalspecs[filt] = []

        for img in AllImglist:
            print("Working on image "+night+" "+img)
            leftover = glob.glob(os.path.splitext(img)[0]+'_*.fits')  #Clean up of some previous attempts if any..
            leftover += glob.glob(os.path.splitext(img)[0]+'.ms.fits')
            if len(leftover) > 0 :
                if not os.path.isdir('Leftover'): os.makedirs('Leftover') 
                for lft in leftover :
                    shutil.move(lft,'Leftover/')

            #Calculate the effective epadu gain for apall
            try:
                ImagesAveraged = int(fits.convenience.getval(img,'NCOMBINE'))  #Reading from Headers
            except KeyError:  #This image is not combined by anything.
                ImagesAveraged = 1
            ImageScaleFactor = 1  # Scaling of image if any was done.
            EffectiveGain = PC.EPADU*ImagesAveraged*ImageScaleFactor

            # Running apall
            iraf.apall(input=img,nfind=1,lower=-15,upper=15,llimit=PC.SPECAPERTURE_LLIMIT,ulimit=PC.SPECAPERTURE_ULIMIT,b_sample=PC.BACKGROUND,background ='fit',weights ='none',readnoi=PC.READNOISE,gain=EffectiveGain,t_function=PC.TRACEFUNC,t_order=PC.TRACEORDER,t_niterate=1,ylevel=PC.SPECAPERTURE,interactive=PC.VER)  #weights= 'variance' seems to be unstable for our high effective gain
            #Extracting the Lamp arc for this spectra as img_arc.fits
            iraf.apall(input=os.path.join(PC.MOTHERDIR,night,Img2Lamp[img]),reference=img,out=os.path.splitext(img)[0]+'_arc',recenter='no',trace='no',background='none',interactive='no')
            #Now reidentify the lines in this spectra
            RepoLamp = 'RepoLamp_'+Img2Filt[img]+'.fits'
            iraf.reidentify(reference=RepoLamp, images=os.path.splitext(img)[0]+'_arc',verbose='yes',interactive=PC.VER)
            #Edit the header of img to add ref lamp
            iraf.hedit(os.path.splitext(img)[0]+'.ms.fits', "REFSPEC1",os.path.splitext(img)[0]+'_arc.fits', add=1, ver=0)
            # dispcor to apply the calibration
            iraf.dispcor(input=os.path.splitext(img)[0]+'.ms.fits',output=os.path.splitext(img)[0]+'_wc.ms.fits')
            #Saving the output file for future
            Filt2finalspecs[Img2Filt[img]].append(os.path.splitext(img)[0]+'_wc.ms.fits')
        
        #At the end of the night Appending the name to the final spectra list of each band
        for filt in Filt2finalspecs.keys():
            with open('FinalwcSpectralistin_'+filt+'.txt','w') as foo:
                foo.write(' \n'.join(Filt2finalspecs[filt])+' \n')

            print('List of final spectra in FinalwcSpectralistin_'+filt+'.txt')
            if PC.SCOMBINE == 'YES':
                try:
                    iraf.scombine(input='@FinalwcSpectralistin_'+filt+'.txt',output=filt+'_avg_'+Filt2finalspecs[filt][0],combine='average',scale='median')
                    print('Averaging the spectra to final output '+filt+'_avg_'+Filt2finalspecs[filt][0])
                except iraf.IrafError as e :
                    print(e)
                    print('ERROR: Could not scombine images in FinalwcSpectralistin_'+filt+'.txt')

            
    print('All nights over...')            
            
            

def SpectralPairSubtraction_subrout(PC):
    """ This will display all the spectra to be extracted to choose the mutual subtraction pairs """
    directories = LoadDirectories(PC,CONF=False)
    for night in directories:
        PC.currentnight = night # Upgating the night directory for using PC.GetFullPath()
        print('Working on night: '+night)
        try:
            #First we load a dictionary of raw images to their filters
            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects.List'),'r') as FiltrFILE :
                Filtrfiledic = dict([(filtset.split()[0],shlex.split(filtset.rstrip())[1]) for filtset in FiltrFILE])  #Dictionary of filterset for each image.

            #Secondly we load a dictionary of raw images to their Calibration Lamp lamb file
            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-BSFinalLamp.List'),'r') as LampFILE :
                Lampfiledic = dict([(Lampset.split()[0],shlex.split(Lampset.rstrip())[1]) for Lampset in LampFILE])  #Dictionary of Lamp file for each image.

            #Secondly, we load a dictionary of Dither Frame to first Raw images
            ProcessedImagesfilename = 'AllObjects-ProcessedCombinedImg.List' if PC.IMGCOMBINE == 'Y' else 'AllObjects-ProcessedImg.List'
            with open(PC.GetFullPath(ProcessedImagesfilename),'r') as DitherFILE :
                DitherFILElines = list(DitherFILE)

            Ditherfiledic = dict([(ditherset.rstrip().split()[1],ditherset.split()[0]) for ditherset in DitherFILElines if len(ditherset.split()) == 2])  #Dictionary of First image of each Dither set.
            #We shall also keep a list of images in proper order.
            Allimglist = [ditherset.rstrip().split()[1]  for ditherset in DitherFILElines if len(ditherset.split()) == 2]

        except IOError as e:
            print('Cannot open the image file list.')
            print(e)
            print('So skipping this directory.')
            print('-'*60)
            continue
        #Output file to write the table of final image, lamp and filter band
        outlog = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'Spectrastoextract_Lamp_BandFilter.txt'),'w')
        #First Lets create the list of filters to iterate through
        FilterList = list(set(Filtrfiledic.values()))
        print("You have %d orders of spectra for this object on this night %s"%(len(FilterList),night))
        OutFilePrefix = raw_input('Enter the prefix of you want for reduced 1d spectra: ')
        for filt in FilterList:
            print("Working on filter : "+filt)
            #List of images with this filter.
            Imglist = [img for img in Allimglist if Filtrfiledic[Ditherfiledic[img]][0] == filt ]
            if len(Imglist) < 16 :
                for i,img in enumerate(Imglist): iraf.display(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,img),i+1)
            else : print('Number of images more that 16, Hence not displaying in ds9')
            if len(Imglist) <= 26 and len(Imglist) != 0:
                ABCDtoimg = dict()
                Anumb = ord('A')
                with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'ABCDtoImageTable_'+filt+'.txt'),'w') as AlphatoFILE :
                    for i,img in enumerate(Imglist):
                        alpha = chr(Anumb+i)
                        ABCDtoimg[alpha] = img
                        print("%s : %s"%(alpha,img))
                        AlphatoFILE.write("%s  %s"%(alpha,img)+' \n')

                print("Enter the pairs to subtract in space separated form ")
                print("For Example an input: AB BC CB D")
                print("Corresponding images produced by subtraction or not are : A-B, B-C, C-B and D")
                print("Note: the final D is not a subtracted image ")
                subpairs = raw_input('Pairs to process: ')
                subpairs = subpairs.split()
                for instr in subpairs:
                    instr = instr.upper()
                    if len(instr) == 2 :
                        Outimg = OutFilePrefix+'_'+filt+'_'+instr[0]+'-'+instr[1]+'.fits'
                        try:
                            iraf.imarith(operand1=PC.GetFullPath(ABCDtoimg[instr[0]]),op="-",operand2=PC.GetFullPath(ABCDtoimg[instr[1]]),result=PC.GetFullPath(Outimg))
                        except iraf.IrafError as e :
                            print(e)
                            print('Skipping to next instruction')
                            continue
                    elif len(instr) == 1 : 
                        Outimg = OutFilePrefix+'_'+filt+'_'+instr[0]+'.fits'
                        shutil.copy(PC.GetFullPath(ABCDtoimg[instr[0]]),PC.GetFullPath(Outimg))
                    else : 
                        print("Could not understand "+instr)
                        continue
                    
                    outlog.write(Outimg+' '+Lampfiledic[Ditherfiledic[ABCDtoimg[instr[0]]]]+' '+filt+' \n')
                #Now Copy the already identified Repository Lamps to this directory.
                print("Copying already identified lines of this filter %s from Repository.."%(filt))
                if not os.path.isdir(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'database')): os.makedirs(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'database'))
                try:
                    shutil.copy(PC.LAMPREPODIR+'/RepoLamp_'+filt+'.fits',os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'RepoLamp_'+filt+'.fits'))
                    shutil.copy(PC.LAMPREPODIR+'/database/idRepoLamp_'+filt,os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'database','idRepoLamp_'+filt))
                except IOError as e:
                    print(e)
                    print("ERROR: Cannot find already identified lines of this filter %s from Repository.."%(filt))
                    print("Before you proceed to next step, do copy the identified line spectra of this filter.")
                    print(" Or remove this image from "+os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'Spectrastoextract_Lamp_BandFilter.txt'))
                    print('-'*10)

            else:
                print("More than 26 images (or no images!!) in single filter? Sorry, you should consider combining them somehow. Or split into two nights directory")
        #End of all filters of this night.
        outlog.close()
    print('All nights over...')                    
        

def NoOfDaophotApertures(aperture_str):
    """ Returns the number of appertures daophot generates from daotphot's multiple apperture syntax """
    NoOfAper = 0
    for aper in aperture_str.split(','):
        if ':' in aper:
            # Generate list of apertures from daophot's closed interval range 
            # syntax ap1:apN:apstep
            ap1, apN, apstep = (float(i) for i in aper.split(':'))
            NoOfAper += int((apN - ap1)/apstep) +1
        else:
            NoOfAper += 1

    return(NoOfAper)

def Photometry(PC):
    """ Does the photometry of images in PC.MOTHERDIR/PC.OUTDIR/Images4Photo.in """
    iraf.noao(_doprint=0)     #Loading packages noao digiohot apphot daophot
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.daophot(_doprint=0)

    iraf.images(_doprint=0) 
    iraf.immatch(_doprint=0) #Loading for xyxymatch, geomap, geotran of coords
    iraf.imfilter(_doprint=0)  #Loading packages for convolution of Gauss

    iraf.apphot.unlearn()
    iraf.daophot.unlearn()
    iraf.phot.unlearn()   #Setting everything to default
    iraf.psf.unlearn()
    iraf.allstar.unlearn()
    iraf.daopars.unlearn()
    iraf.datapar.unlearn()
    iraf.fitskypar.unlearn()
    iraf.photpar.unlearn()
    iraf.findpars.unlearn()

    try:
        imgfile = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'Images4Photo.in'),'r')
    except IOError as e:
        print('Cannot open Images4Photo.in file. Run Task #6 ')
        print(e)
        print('-'*60)
        traceback.print_exc(file=sys.stdout)
        print('-'*60)
        sys.exit(1)

    # Setting flag by checking whether the size of qphotinput.txt is Zero or not.
    if os.stat(os.path.join(PC.MOTHERDIR,PC.OUTDIR,"qphotinput.txt"))[6]!=0 : QPHOT_todo='Y'
    else : QPHOT_todo='N'


    imgNo = 0
    FirstImageName = None
    for imgline in imgfile :
        imgline = imgline.rstrip()
        imgline = shlex.split(imgline)
        if FirstImageName is None : FirstImageName = os.path.abspath(imgline[0])  #Saving the first image name of this photometry run
        wdir = imgline[0].split('/')[-2]
        if wdir != os.getcwd().split('/')[-1] : #If we are not already in img dir
            iraf.cd(os.path.join(PC.MOTHERDIR,PC.OUTDIR))  #Going back to output directory
            DIRtogo = "/".join(imgline[0].split('/')[:-1]) #Now going to dir of img
            iraf.cd(DIRtogo)
            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,PC.OUTPUTFILE),'a') as foo:    #Appending into tbe log file The beginning of a new directory
                foo.write(wdir+' ---------------------------------------- \n')  # '-'*40  To mark beginning of a Directory

        
        img = imgline[0].split('/')[-1]
        filterr = imgline[1]
        intime = float(imgline[2])
        threshold = float(imgline[3])
        try :
            DateStr = imgline[4]
            Year,month,day = DateStr.split('-')
        except KeyError:    #Old data didn't have this header
            if '-' in wdir:
                Year,month,day = wdir.split('-')
            else:
                Year,month,day=wdir[:4],wdir[4:6],wdir[6:8]  # Will work only if the directory was named YYYYMMDD

        StartUTStr = imgline[5]
        h,m,s = StartUTStr.split(':')
        StartUT = int(h)*60*60+int(m)*60+int(float(s))   #Converting to seconds

        obstime = Time('-'.join([Year,month,day])+' '+StartUTStr, format='iso', scale='utc')
        #Calculate the effective epadu gain for daophot's photometry error estimation
        try :
            ImagesAveraged = int(fits.convenience.getval(img,'NCOMBINE'))  #Reading from Headers
        except KeyError:  #This image is not combined by anything.
            ImagesAveraged = 1
        ImageScaleFactor= 1  # If the image was scaled give the factor here.
        EffectiveGain=PC.EPADU*ImagesAveraged*ImageScaleFactor
        print(wdir, img)
        TrueSigma = 0.01      #Will be updated later for every img. Here it should be very small quantity
        yxdim = fits.getdata(img).shape  #Storing the (Ymax,Xmax) of image
        leftover = glob.glob(img+'?*')
        if len(leftover) > 0 :
            if not os.path.isdir('Leftover'): os.makedirs('Leftover') 
            for lft in leftover :
                shutil.move(lft,'Leftover/')

        #Calling Sextracter And find coordinates
        if not os.path.isfile(os.path.join(PC.MOTHERDIR,PC.OUTDIR,"sextractor.sex")) :  #Incase the parameter file is not already created
            Sextractor_subrout(PC,img=imgline[0])


        # Calling sextractor for this new image
        N = 30  #Number of bright stars to take
        Sextractor_subrout(PC,img=img,N=N,OutFilePrefix='Bright',OutDir='.')

        #Runing xyxymatch and geo map and geotran to create new SourceT.coo , GoodStarsT.coo, BlankSky.coo
        Nmatch = 32
        num_lines = 0
        while num_lines < PC.XYMATCHMIN :       # the number of stars it should atlest mach is set in .conf file XYMATCHMIN= 6....
            if os.path.isfile(img+"xymatch.out") :os.remove(img+"xymatch.out")
            Nmatch = Nmatch-2
            iraf.xyxymatch.unlearn()
            iraf.xyxymatch(input='Bright{0}.coo'.format(N),reference=os.path.join(PC.MOTHERDIR,PC.OUTDIR,"FirstImageTop{0}.coo".format(N)),output=img+"xymatch.out", toler=3, matching="triangles",nmatch=Nmatch)
            # Making sure atleast a few stars were matched. otherwise geoxytran will exit with error.
            os.system("awk '{if ($1 >0){print $0}}' "+img+"xymatch.out > matchedstars.txt") #Removing headers
            num_lines = sum(1 for line in open("matchedstars.txt"))  #Counting number of lines in the xymatch output. it should be more than 15 (size of header)
            print("Number of stars Matched= "+str(num_lines))
            if Nmatch < PC.XYMATCHMIN : 
                print("Failed to find the coordinates for "+img)
                print("We need to find the transformation interactively")
                Inpfirstimg = raw_input("Enter the full path to first image using which coords was generated with sextractor (default: {0}) : ".format(FirstImageName)).strip(' ')
                if Inpfirstimg : FirstImageName=Inpfirstimg
                if os.path.isfile(FirstImageName): 
                    print("Running the xyxy match interactively... Select three same points from both images by pressing a")
                    print("IMPORTANT: First select 3 stars in Frame 1 of ds9. Second image to select is in Frame 2")
                    iraf.display(img,2)
                    iraf.display(FirstImageName,1)
                    if os.path.isfile(img+"xymatch.out") :os.remove(img+"xymatch.out")
                    iraf.xyxymatch.unlearn()
                    iraf.xyxymatch(input="Bright30.coo",reference=os.path.join(PC.MOTHERDIR,PC.OUTDIR,"FirstImageTop30.coo"),output=img+"xymatch.out", toler=5, matching="tolerance",nmatch=3*PC.XYMATCHMIN,interactive="yes")
                    os.system("awk '{if ($1 >0){print $0}}' "+img+"xymatch.out > matchedstars.txt") #Removing headers
                    num_lines = sum(1 for line in open("matchedstars.txt"))  #Counting number of lines in the xymatch output. it should be more than 15 (size of header)
                    print("Number of stars Matched= "+str(num_lines))
                    break
                else:
                    print("ERROR: Cannot find the file :"+FirstImageName)
                    print("Enter the correct full path to the file again after this attempt.")
        
        GeoMapfitgeometry = "general"
        if num_lines < 6 : GeoMapfitgeometry="rotate"  # Just XY shift and rotate
        iraf.geomap(input=img+"xymatch.out", database=img+"rtran.db", xmin=1, xmax=yxdim[1], ymin=1, ymax=yxdim[0], interactive=0,fitgeometry=GeoMapfitgeometry)
        iraf.geoxytran(input=os.path.join(PC.MOTHERDIR,PC.OUTDIR,"GoodStars.coo"), output=img+"GoodStarsT.coo",database=img+"rtran.db",transforms=img+"xymatch.out")
        iraf.geoxytran(input=os.path.join(PC.MOTHERDIR,PC.OUTDIR,"Source.coo"), output=img+"SourceT.coo",database=img+"rtran.db",transforms=img+"xymatch.out")
        iraf.geoxytran(input=os.path.join(PC.MOTHERDIR,PC.OUTDIR,"BlankSky.coo"), output=img+"BlankSky.coo",database=img+"rtran.db",transforms=img+"xymatch.out")
        if QPHOT_todo == 'Y' :
            iraf.geoxytran(input=os.path.join(PC.MOTHERDIR,PC.OUTDIR,"qphotinput.txt"), output=img+"qphotinput.txt",database=img+"rtran.db",transforms=img+"xymatch.out")


        # Sanity check: To remove any new coordinates calculated lying outside image in *.coo 
        coofileLIST = ['GoodStarsT.coo','SourceT.coo','BlankSky.coo']
        if QPHOT_todo == 'Y' : coofileLIST.append('qphotinput.txt') 
        for coofile in coofileLIST :
            fooIN = open(img+coofile,'r')
            fooOUT = open(img+coofile+'TEMP','w')
            for star in fooIN:
                if float(star.split()[0]) > 1 and float(star.split()[0]) < yxdim[1] and float(star.split()[1]) > 1 and float(star.split()[1]) < yxdim[0] : fooOUT.write(star)
                else: print(star +": Outside the image field \n")
            fooIN.close()
            fooOUT.close()
            os.rename(img+coofile+'TEMP',img+coofile)

        #---------------------------------
        # Due to small error in calculation of star position, we need to create a more accurate GoodStars.coo and Source.coo
        # Plus we have to remove any saturated stars also.
        imx = iraf.imexam(input=img,frame=1,use_display=0,defkey='a',imagecur=img+'GoodStarsT.coo',Stdout=1)
        foo = open(img+'GoodStars.coo','w')    #Creating good stars coords files
        starlist = " "
        i = 2
        while i < len(imx) :
            if (imx[i+1].split()[4] != 'INDEF') and (min(float(imx[i+1].split()[10]),float(imx[i+1].split()[9])) > np.abs(float(imx[i+1].split()[10])-float(imx[i+1].split()[9]))) and (float(imx[i+1].split()[4])+float(imx[i+1].split()[3])) < float(PC.DATAMAX) : # (Peak-sky) is not INDEF and min(MoffetFWHM,directFWHM) > abs(directFWHM-MoffetFWHM) and  Peak is not saturated
                foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                starlist = starlist+str(i/2)+"_"   #Saving the string of good stars survived.
            else : print('Discarded: '+str(i/2)+' th number star not good of '+DIRtogo+' '+img)
            i = i+2
        foo.close()

        try :
            imx = iraf.imexam(input=img,frame=1,use_display=0,defkey='a',imagecur=img+'SourceT.coo',Stdout=1)
            Xprim = eval(imx[2].split()[0])  
            Yprim = eval(imx[2].split()[1])
            with open(img+'Source.coo','w') as foo :    #Creating text file containing coords of primary interest
                i = 2
                while i < len(imx) :
                    foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                    i = i+2

        except iraf.IrafError as e :
            print('Iraf Error: going forward with first estimated Source.coo')
            shutil.copy(img+'SourceT.coo',img+'Source.coo')
    
        #---------------------------------END of recalculation of coordinates------


        #Creating all the gauss convolved images if the CONVOLVEIMG variable is _not_ set to NO
        OriginalIMG = img
        convIMGS = [img]
        if PC.CONVOLVEIMG != 'NO' :  # If the PC.CONVOLVEIMG variable is not set to NO
            for si in eval(PC.CONVOLVEIMG) :
                iraf.gauss(input=img,output=img+'_'+str(si)+'.fits',sigma=si)
                convIMGS.append(img+'_'+str(si)+'.fits')       #List of all convolved images on which we have to do photometry
        # Now the loop of  doing the photometry for all the convolved) images
        for img in convIMGS :
            imx = iraf.imexam(input=img,frame=1,use_display=0,defkey='a',imagecur=OriginalIMG+'GoodStars.coo',Stdout=1)
            print('\n'.join(imx)+'\n')           #DEBUGGING---------------------------------------**
            #Calculating median FWHM, ellipticity and position angle
            fwhmlist = []
            ellipticitylist = []
            positionanglelist = []
            i = 3
            while i < len(imx) :               
                fwhmlist.append(eval(imx[i].split()[10]))
                if imx[i].split()[5] != 'INDEF' : ellipticitylist.append(eval(imx[i].split()[5]))
                if imx[i].split()[6] != 'INDEF' : positionanglelist.append(eval(imx[i].split()[6]))
                i = i+2
            #Median FWHM is
            fwhm = np.median(fwhmlist)  
            print('Setting value of FWHM =' + str(fwhm))
            ellipticity = np.median(ellipticitylist)
            print('Setting value of ellipticity =' + str(ellipticity))
            positionangle = np.median(positionanglelist)%180   #addition modulo 180 to keep values in 0 to 180 range
            print('Setting value of positionangle %180 =' + str(positionangle))
            #Calculating sky mean and stdev
            imx = iraf.imexam(input=img,frame=1,use_display=0,defkey='m',imagecur=OriginalIMG+'BlankSky.coo',Stdout=1)
            print('\n'.join(imx)+'\n')            #DEBUGGING--------------------------------------**
            #Calculating average Sigma and mean of sky
            skylist = []
            sigmalist = []
            i = 1
            while i < len(imx) :               
                skylist.append(eval(imx[i].split()[2]))
                sigmalist.append(eval(imx[i].split()[4]))
                i = i+1
            #Average mean,sigma and datamin are
            mean = np.median(skylist)
            sigma = np.median(sigmalist)
            datamin = mean - 10*max(sigma,TrueSigma)  #Setting lower limit to sigma from going less than TrueSigma of original img
            print('Mean sky = '+str(mean))
            print('Mean sigma = '+str(sigma))
            print('Datamin = '+str(datamin))

            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,PC.OUTPUTFILE),'a') as foo:    #Appending into the log file to write output of photometry
                foo.write(img +'  '+str(round(fwhm,2)) + '  "'+filterr+'"  '+str(intime) +'  '+str(StartUT)+'  '+ str(round(mean,3)) +'  ' + str(round(sigma,3)) +'  '+str(round(datamin,3)) + ' | ')

            #Creating a full table with same information
            FullOutputTable = table.Table([[obstime.iso.split()[0]], [obstime.iso.split()[1]], [obstime.jd]], names=('Date', 'UT', 'JD'),masked=True)
            FullOutputTable['Filter'] = filterr
            FullOutputTable['ExpTime'] = intime
            FullOutputTable['FWHM'] = round(fwhm,2)
            FullOutputTable['Ellipticity'] = round(ellipticity,3)
            FullOutputTable['PositionAngle'] = round(positionangle,1)
            FullOutputTable['SkyMean'] = round(mean,3)
            FullOutputTable['SkySigma'] = round(sigma,3)
            FullOutputTable['Datamin'] = round(datamin,3)
            FullOutputTable['ReadNoise'] = PC.READNOISE
            FullOutputTable['GainEffective'] = EffectiveGain
            FullOutputTable['SkyAnnulus'] = eval(PC.ANNULUS)
            FullOutputTable['SkyDannulus'] = eval(PC.DANNULUS)
            FullOutputTable['Threshold'] = threshold

            #Now starts photometry processes....
            # First is daopar, Press :w then :q to continue if everything is fine
            psfradi = 4*fwhm +1

#            aperture = eval(PC.APERTURE)  # 4*fwhm
            aperture = ','.join([':'.join([str(eval(x)) for x in y.split(':')]) for y in PC.APERTURE.split(',')])  #For multiple apperture syntax

            iraf.daopars.setParam('matchrad',fwhm)
            iraf.daopars.setParam('psfrad',psfradi)
            iraf.daopars.setParam('fitrad',fwhm)

            iraf.datapar.setParam('fwhmpsf',fwhm)
            iraf.datapar.setParam('sigma',sigma)
            iraf.datapar.setParam('datamin',datamin)
            iraf.datapar.setParam('datamax',PC.DATAMAX)
            iraf.datapar.setParam('readnoi',PC.READNOISE)
            iraf.datapar.setParam('epadu',EffectiveGain)
            iraf.datapar.setParam('itime',intime)
#            iraf.datapar.setParam('ifilter',filterr)

            iraf.fitskypar.setParam('annulus',eval(PC.ANNULUS))
            iraf.fitskypar.setParam('dannulu',eval(PC.DANNULUS))
            
            iraf.photpar.setParam('apertur',aperture)

            iraf.findpars.setParam('threshold',threshold)
            iraf.findpars.setParam('ratio',1-ellipticity)  # ellipticity = 1 - b/a
            iraf.findpars.setParam('theta',positionangle)

            if OriginalIMG == img :
                TrueSigma = sigma   #Setting the correct sigma of sky
                iraf.daofind(image=img,output="default",verify=PC.VER)
                # Check GoodStars and Source are identified in daofind. else add forcefully
                cootable = ascii.read(img+'.coo.1',format='daophot')
                with open(OriginalIMG+'GoodStars.coo','r') as goodstarsFILE :
                    for SiD,goodstarXY in zip(starlist.strip().strip('_').split('_'), goodstarsFILE):
                        gsX,gsY = goodstarXY.rstrip().split()
                        gsX = eval(gsX)
                        gsY = eval(gsY)
                        NoOfStarcoo = len(cootable[((cootable['XCENTER']-gsX)**2<4) & ((cootable['YCENTER']-gsY)**2<4)])
                        if NoOfStarcoo == 0: # Star missing in output coo file
                            print('**ERROR** : Something wrong with Good Star No:{0}'.format(SiD))
                            print('It is missing from the daofind output .coo.1')
                            print('Adding it forcefully into the .coo.1 file')
                            with open(img+'.coo.1','a') as CooOfstarsFILE :
                                CooOfstarsFILE.write(' {0} {1} \n'.format(gsX,gsY))
                        elif NoOfStarcoo > 1 :
                            print('**Warning** : {0} Multiple detection of Good Star No:{1}'.format(NoOfStarcoo,SiD))
                # Simillarly also add any missing Primary Source
                with open(OriginalIMG+'Source.coo','r') as SourcestarsFILE :
                    for SiD,SourcestarXY in enumerate(SourcestarsFILE):
                        psX,psY = SourcestarXY.rstrip().split()
                        psX = eval(psX)
                        psY = eval(psY)
                        NoOfStarcoo = len(cootable[((cootable['XCENTER']-psX)**2<4) & ((cootable['YCENTER']-psY)**2<4)])
                        if NoOfStarcoo == 0: # Star missing in output coo file
                            print('**ERROR** : Something wrong with Primary Source No:{0}'.format(SiD+1))
                            print('It is missing from the daofind output .coo.1')
                            print('Adding it forcefully into the .coo.1 file')
                            with open(img+'.coo.1','a') as CooOfstarsFILE :
                                CooOfstarsFILE.write(' {0} {1} \n'.format(psX,psY))
                        elif NoOfStarcoo > 1 :
                            print('**Warning** : {0} Multiple detection of Primary Source No:{1}'.format(NoOfStarcoo,SiD+1))
                            
                
            else :
                shutil.copy(OriginalIMG+'.coo.1',img+'.coo.1')
            #Going forward to do phot
            iraf.phot(image=img,coords="default",output="default",interactive='no',verify=PC.VER)

            magtable = ascii.read(img+'.mag.1',format='daophot')
            aperheader = magtable.meta['keywords']['APERTURES']['value']
            NoOfAper = NoOfDaophotApertures(aperheader)  #Number of daophot appertures
            FullOutputTable['NoOfAper'] = NoOfAper
            with open(OriginalIMG+'GoodStars.coo','r') as goodstarsFILE :
                tablelist = []
                for SiD,goodstarXY in zip(starlist.strip().strip('_').split('_'), goodstarsFILE):
                    gsX,gsY = goodstarXY.rstrip().split()
                    gsX = eval(gsX)
                    gsY = eval(gsY)
                    Starsmagtable = magtable[((magtable['XCENTER']-gsX)**2<4) & ((magtable['YCENTER']-gsY)**2<4)]
                    if len(Starsmagtable) > 1 :
                        print('**WARNING** : More than two stars detected for Good Star No:{0} in magfile'.format(SiD))
                        print('Taking only the first detection to Output table')
                        Starsmagtable = Starsmagtable[0:1]

                    tablelist.append(Starsmagtable['XCENTER','YCENTER','MAG','ID'])
                    #Adding to Full output table
                    FullOutputTable['GoodStarS'+SiD+'_XCENTER']= Starsmagtable['XCENTER']
                    FullOutputTable['GoodStarS'+SiD+'_YCENTER']= Starsmagtable['YCENTER']
                    FullOutputTable['GoodStarS'+SiD+'_RAPERT1']= Starsmagtable['RAPERT']
                    FullOutputTable['GoodStarS'+SiD+'_MAG1']= Starsmagtable['MAG']
                    FullOutputTable['GoodStarS'+SiD+'_MERR1']= Starsmagtable['MERR']
                    for i in range(2,NoOfAper+1): #For appertures more than 1
                        FullOutputTable['GoodStarS'+SiD+'_RAPERT'+str(i)]= Starsmagtable['RAPERT'+str(i)]
                        FullOutputTable['GoodStarS'+SiD+'_MAG'+str(i)]= Starsmagtable['MAG'+str(i)]
                        FullOutputTable['GoodStarS'+SiD+'_MERR'+str(i)]= Starsmagtable['MERR'+str(i)]

            if len(tablelist) > 1 : 
                goodstarsTable = table.vstack(tablelist)
            else : 
                goodstarsTable = tablelist[0]

            if PC.DOPSF == 'YES' :  # IF PSF photometry has to be done...
                #Creating the imcommands file by finding Star IDs
#                os.system(PC.MOTHERDIR+'/Finding_StarID_Curser_File.sh ' + img +' '+OriginalIMG+'GoodStars.coo' )  #OLD Way...
                with open('icommands.in','w') as icomFILE :
                    icomFILE.write(':a '+'\n:a '.join([str(sid) for sid in goodstarsTable['ID']])+'\n')
                    icomFILE.write('f \nw \nq \n') # Adding lines f w q at the end.

                print("Doing psf, Non-interactively.. Using Coords of good star")
                iraf.psf(image=img, pstfile="", photfile="default", psfimage="default", opstfile="default", groupfil="default", icommands='icommands.in', verify=PC.VER)
            #    print ("Doing psf, Interactively.. Use the a a ... f w q  sequence..")
            #    iraf.psf(image=img, pstfile="", photfile="default", psfimage="default", opstfile="default", groupfil="default")

                iraf.allstar(image=img, photfile="default", psfimage="default", allstarf="default", rejfile="default", subimage="default" ,verify=PC.VER )
                print("Psf photometry over")
                print("--------------------------------------")

            #Doing the phot again on Source, Just in case Daofind didn't detect it and Good stars...
            iraf.datapar.setParam('datamin',mean-5*max(sigma,TrueSigma))
            iraf.phot(image=img,coords=OriginalIMG+'Source.coo',output="default",interactive='no',verify=PC.VER)
            Sourcemagtable = ascii.read(img+'.mag.2',format='daophot')
            aperheader = Sourcemagtable.meta['keywords']['APERTURES']['value']
            NoOfAper = NoOfDaophotApertures(aperheader)  #Number of daophot appertures
            #Adding to Full output table
            for SiD, rows in enumerate(Sourcemagtable):  #Source Mags
                SiD = str(SiD +1)
                FullOutputTable['SourceS'+SiD+'_XCENTER']= rows['XCENTER']
                FullOutputTable['SourceS'+SiD+'_YCENTER']= rows['YCENTER']
                FullOutputTable['SourceS'+SiD+'_RAPERT1']= rows['RAPERT']
                FullOutputTable['SourceS'+SiD+'_MAG1']= rows['MAG']
                FullOutputTable['SourceS'+SiD+'_MERR1']= rows['MERR']
                for i in range(2,NoOfAper+1): #For appertures more than 1
                    FullOutputTable['SourceS'+SiD+'_RAPERT'+str(i)]= rows['RAPERT'+str(i)]
                    FullOutputTable['SourceS'+SiD+'_MAG'+str(i)]= rows['MAG'+str(i)]
                    FullOutputTable['SourceS'+SiD+'_MERR'+str(i)]= rows['MERR'+str(i)]

#            iraf.phot(image=img,coords=OriginalIMG+'GoodStars.coo',output="default",verify=PC.VER)
#            SecondPhotresults=iraf.txdump(textfiles=img+".mag.2",fields="XCENTER,YCENTER,MAG",expr="yes",Stdout=1)
#            SecondPhotresults.extend(iraf.txdump(textfiles=img+".mag.3",fields="XCENTER,YCENTER,MAG",expr="yes",Stdout=1))

            iraf.hedit(img, "intime", intime, add=1, ver=0)
            #Doing qphot at all the points in qphotinput.txt with the corresponding parameters.
            NoOfQphots = 0
            if QPHOT_todo == 'Y' :  #If there exist some qphot sources
                foo = open(OriginalIMG+"qphotinput.txt",'r')
                for qphotobj in foo:
                    qphotobj = qphotobj.rstrip()
                    obj = qphotobj.split()
                    with open('qphotSource.Tcoo','w') as foo2 :
                        foo2.write(obj[0]+'  '+obj[1])

                    iraf.qphot(image=img , coords='qphotSource.Tcoo', cbox=5, annulus=obj[3], dannulus=obj[4], aperture=obj[2], exposur="intime", epadu=EffectiveGain ,interactive=0 )
                    NoOfQphots += 1
                
                foo.close()
            #Now, Writing the Mag to output file
            foo = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,PC.OUTPUTFILE),'a')
#            os.system(PC.MOTHERDIR+'/Creating_Log_File.sh '+img+' '+OriginalIMG+'GoodStars.coo'+' '+OriginalIMG+'Source.coo'+' '+PC.MOTHERDIR+'/'+PC.OUTPUTFILE ) 

            #First append the qphot magnitudes to the Photometry output file
            for i in range(NoOfQphots):
                filemag = img+'.mag.'+str(i+3)
                qphottable = ascii.read(filemag,format='daophot')
                foo.write(' {0}'.format(qphottable['MAG'][-1])) #Writing the last Mag in file

                aperheaderQphot = qphottable.meta['keywords']['APERTURES']['value']
                NoOfQphotAper = NoOfDaophotApertures(aperheaderQphot)  #Number of qphot appertures
                #Adding to Full output table
                FullOutputTable['QphotS'+str(i+1)+'_SKYANNUL'] = float(qphottable.meta['keywords']['ANNULUS']['value'])
                FullOutputTable['QphotS'+str(i+1)+'_SKYDANNUL'] = float(qphottable.meta['keywords']['DANNULUS']['value'])
                FullOutputTable['QphotS'+str(i+1)+'_XCENTER'] = qphottable['XCENTER'][-1]
                FullOutputTable['QphotS'+str(i+1)+'_YCENTER'] = qphottable['YCENTER'][-1]
                FullOutputTable['QphotS'+str(i+1)+'_RAPERT1'] = qphottable['RAPERT'][-1]
                FullOutputTable['QphotS'+str(i+1)+'_MAG1'] = qphottable['MAG'][-1]
                FullOutputTable['QphotS'+str(i+1)+'_MERR1'] = qphottable['MERR'][-1]
                for j in range(2,NoOfQphotAper+1): #For appertures more than 1
                    FullOutputTable['QphotS'+str(i+1)+'_RAPERT'+str(j)] = qphottable['RAPERT'+str(j)][-1]
                    FullOutputTable['QphotS'+str(i+1)+'_MAG'+str(j)] = qphottable['MAG'+str(j)][-1]
                    FullOutputTable['QphotS'+str(i+1)+'_MERR'+str(j)] = qphottable['MERR'+str(j)][-1]
                

#            magfiles=glob.glob(img+'.mag.*')
#            magfiles.sort()
#            for filemag in magfiles:
#                if eval(filemag.split('.')[-1]) > 2 : # The qphot output files
#                    qphottable=ascii.read(filemag)
#                    foo.write(' {0}'.format(qphottable['MAG'][-1])) #Writing the last Mag in file

            foo.write(' | ')  #Adding a separator after qphot mags
            # If PSF photometry as done, adding those mags to the file.
            if PC.DOPSF == 'YES' :  # IF PSF photometry was done...
                #First the mags of Source stars
                alstable = ascii.read(img+'.als.1',format='daophot')
                FullOutputTable['PSFRAD'] = float(alstable.meta['keywords']['PSFRAD']['value'])
                FullOutputTable['PSFFITRAD'] = float(alstable.meta['keywords']['FITRAD']['value'])
                with open(OriginalIMG+'Source.coo','r') as SourcestarsFILE :
                    tablelist = []
                    for SiD, sourstarXY in enumerate(SourcestarsFILE):
                        SiD = str(SiD +1)
                        ssX,ssY = sourstarXY.rstrip().split()
                        ssX = eval(ssX)
                        ssY = eval(ssY)
                        Starsalstable = alstable[((alstable['XCENTER']-ssX)**2<4) & ((alstable['YCENTER']-ssY)**2<4)]
                        tablelist.append(Starsalstable['XCENTER','YCENTER','MAG'])
                        #Adding to Full output table
                        FullOutputTable['SourceS'+SiD+'_PSFXCENTER'] = Starsalstable['XCENTER']
                        FullOutputTable['SourceS'+SiD+'_PSFYCENTER'] = Starsalstable['YCENTER']
                        FullOutputTable['SourceS'+SiD+'_PSFMAG'] = Starsalstable['MAG']
                        FullOutputTable['SourceS'+SiD+'_PSFMERR'] = Starsalstable['MERR']
                        FullOutputTable['SourceS'+SiD+'_PSFCHI'] = Starsalstable['CHI']
                        

                if len(tablelist) > 1 : 
                    SourcestarsALSTable = table.vstack(tablelist)
                else:
                    SourcestarsALSTable = tablelist[0]

                for rows in SourcestarsALSTable: 
                    foo.write(' %f %f %f'%(rows['XCENTER'],rows['YCENTER'],rows['MAG']))
                foo.write(' | ')  #Adding a separator after Source Mags
                #Now the psf magnitudes of the good stars
                with open(OriginalIMG+'GoodStars.coo','r') as goodstarsFILE :
                    tablelist = []
                    for SiD,goodstarXY in zip(starlist.strip().strip('_').split('_'), goodstarsFILE):
                        gsX,gsY = goodstarXY.rstrip().split()
                        gsX = eval(gsX)
                        gsY = eval(gsY)
                        Starsalstable = alstable[((alstable['XCENTER']-gsX)**2<4) & ((alstable['YCENTER']-gsY)**2<4)]
                        tablelist.append(Starsalstable['XCENTER','YCENTER','MAG'])
                        #Adding to Full output table
                        FullOutputTable['GoodStarS'+SiD+'_PSFXCENTER'] = Starsalstable['XCENTER']
                        FullOutputTable['GoodStarS'+SiD+'_PSFYCENTER'] = Starsalstable['YCENTER']
                        FullOutputTable['GoodStarS'+SiD+'_PSFMAG'] = Starsalstable['MAG']
                        FullOutputTable['GoodStarS'+SiD+'_PSFMERR'] = Starsalstable['MERR']
                        FullOutputTable['GoodStarS'+SiD+'_PSFCHI'] = Starsalstable['CHI']

                if len(tablelist) > 1 : 
                    goodstarsALSTable = table.vstack(tablelist)
                else :
                    goodstarsALSTable = tablelist[0]

                for rows in goodstarsALSTable:
                    foo.write(' %f %f %f'%(rows['XCENTER'],rows['YCENTER'],rows['MAG']))
                foo.write(' | ')  #Adding a separator after Good stars X Y Mags
                
            else:
                foo.write(' | | ')
            
            # Writing the pure phot results we calculated into the list before closing the line
            for rows in Sourcemagtable:  #Source Mags
                foo.write(' %f %f %f'%(rows['XCENTER'],rows['YCENTER'],rows['MAG']))
            foo.write(' | ')  #Adding a separator after Source Mags
            for rows in goodstarsTable:  #Good Stars Mags
                foo.write(' %f %f %f'%(rows['XCENTER'],rows['YCENTER'],rows['MAG']))
            foo.write(' | ')  #Adding a separator after Good Star Mags

            foo.write(' '+starlist+' \n') # Ending this image line with Good star's ID.
            foo.close()

            # Now also append the Full output table also to an ascii file.
            FullOutputTable['Image'] = img
            try :
                PreviousFullTable = ascii.read(os.path.join(PC.MOTHERDIR,PC.OUTDIR,PC.OUTPUTFILE+'_FullOutput.txt'),delimiter='|',format='commented_header',fill_values=('--','0'))
            except IOError :
                print('No previous photometry output found, hence we will be creating a new file.')
                OutputTableToWrite = FullOutputTable
            else :
                OutputTableToWrite = table.vstack([PreviousFullTable,FullOutputTable], join_type='outer')
            # Writing the final appended table
            ascii.write(OutputTableToWrite, os.path.join(PC.MOTHERDIR,PC.OUTDIR,PC.OUTPUTFILE+'_FullOutput.txt'),delimiter='|',format='commented_header')
           
            print ("Photometry of "+img+" over. \n Now proceeding to next image")
            #END of the photometry of convolved images set..
        imgNo = imgNo+1
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,PC.OUTPUTFILE),'a') as foo :    #Appending into the log file to write output of photometry
            foo.write('-------------------------------------------- \n')  # '-'*44  To mark end of an image


    #All photometry over
    imgfile.close()
    print("Great...Photometry of all "+str(imgNo)+" images are over...")
    print("Enjoy!!! ---------------------indiajoe@gmail.com ")

def is_number(s):   # A function to check whether string s is a number or not.
    try:
        float(s)
        return True
    except ValueError:
        return False

def Sextractor_subrout(PC,img=None,N=30,OutFilePrefix='FirstImageTop',OutDir=None):
    """ Calls the Sextractor and create the sextractor parameter files if it doesn't already exists. And also create coord file of the brightest N=30 number of stars."""
    try :
        subprocess.call(['sex','--version'])
        SEXTRACTOR = 'sex'
    except OSError:
        try :
            subprocess.call(['sextractor','--version'])
            SEXTRACTOR = 'sextractor'
        except OSError:            
            print('ERROR: Cannot find the command: sex')
            print('SExtractor needs to be installed before running this task')
            sys.exit(1)

    if OutDir is None:
        OutDir = os.path.join(PC.MOTHERDIR,PC.OUTDIR)


    if not os.path.isfile(os.path.join(PC.MOTHERDIR,PC.OUTDIR,"sextractor.sex")) : #If a config file doesn't exist already
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'sextractor.sex'),'w') as sexConfigFile:
            subprocess.call([SEXTRACTOR,'-d'],stdout=sexConfigFile)
        #Change PIXEL_SCALE to 0.3
        subprocess.call(['sed','-i',r's/^\(PIXEL_SCALE\s*\)\([0-9]*\.[0-9]*\)/\10.3/',os.path.join(PC.MOTHERDIR,PC.OUTDIR,'sextractor.sex')])
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'default.conv'),'w') as convolutionFile:
            convolutionFile.write("""CONV NORM\n# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.\n1 2 1\n2 4 2\n1 2 1\n""")
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'default.param'),'w') as sexCatParamFile:
            sexCatParamFile.write('\n'.join(['NUMBER','FLUXERR_ISO','FLUX_AUTO','FLUXERR_AUTO','X_IMAGE','Y_IMAGE','FLAGS'])+'\n')
        
        print("Sextractor Config file sextractor.sex, default.parm and default.conv created. \n If required u can edit it before calling Photometry")

    if img is None : # If No img is given, then using the first image in Images4Photo.in file
        try:
            imgfile = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'Images4Photo.in'),'r')
        except IOError,e:
            print('Cannot open Images4Photo.in file. Run Task #6 ')
            print(e)
            print('-'*60)
            traceback.print_exc(file=sys.stdout)
            print('-'*60)
            sys.exit(1)

        imgline = imgfile.readline()  #First line only
        imgline = imgline.rstrip()
        img = imgline.split()[0]
        imgfile.close()

    subprocess.call([SEXTRACTOR,img,"-c",os.path.join(PC.MOTHERDIR,PC.OUTDIR,"sextractor.sex"),"-PARAMETERS_NAME",os.path.join(PC.MOTHERDIR,PC.OUTDIR,"default.param"),"-FILTER_NAME",os.path.join(PC.MOTHERDIR,PC.OUTDIR,"default.conv"),'-CATALOG_NAME',os.path.join(OutDir,'SextractorOutput.cat')])

    SExtractCat = ascii.read(os.path.join(OutDir,'SextractorOutput.cat'),format='sextractor')
    # Selecting only good stars without any major problems
    GoodStarCat = SExtractCat[SExtractCat['FLAGS']<2]  # Flag 0 is good and 1 is contaminated by less than 10% in flux by neighbor
    #Sort in descending order of Flux
    GoodStarCat.sort('FLUX_AUTO')
    GoodStarCat.reverse()
    #Write X and Y coordinates of First N number of brightest stars in text file
    GoodStarCat['X_IMAGE','Y_IMAGE'][:N].write(os.path.join(OutDir,OutFilePrefix+'{0}.coo'.format(N)),format='ascii.no_header')

    print("Brightest {0} stars coordinates of first image created in {1}{0}.coo".format(N,os.path.join(OutDir,OutFilePrefix)))



def Star_sky_subrout(PC,img=None) :
    """ Opens the image and create Source.coo, GoodStars.coo, BlankSky.coo, Polygon.coo files"""
    backupPWD = os.getcwd()
    iraf.cd(os.path.join(PC.MOTHERDIR,PC.OUTDIR))  #Going to output directory of this run.

    if img is None : # If No img is given, then using the first image in Images4Photo.in file
        try:
            imgfile = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'Images4Photo.in'),'r')
        except IOError as e:
            print('Cannot open Images4Photo.in file. Run Task #6 ')
            print(e)
            print('-'*60)
            traceback.print_exc(file=sys.stdout)
            print('-'*60)
            sys.exit(1)

        imgline = imgfile.readline()  #First line only
        imgline = imgline.rstrip()
        img = imgline.split()[0]
        imgfile.close()

    if not ( os.path.isfile("Source.coo") and os.path.isfile("GoodStars.coo") and os.path.isfile("BlankSky.coo") )  : #If the .coo files doesn't exist already
        iraf.display(img,1)
        print ('\n For taking coordinates of Source. Press _a_ over Primary Sources (atleast one).')
        imx = iraf.imexam(Stdout=1)
        with open('Source.coo','w') as foo :    #Creating text file containing coords of science sources of primary interest
            i = 2
            while i < len(imx) :               
                foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                i = i+2

        print ('\n For taking coordinates of good stars. Press _a_ over a few good stars. \n Non saturated among them will be used for psf fitting.')
        print ('IMP: Press coordinate of Stars in standard required order')
        imx = iraf.imexam(Stdout=1)
        with open('GoodStars.coo','w') as foo :    #Creating good stars coords files
            i = 2
            while i < len(imx) :               
                foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                i = i+2

        shutil.copy('GoodStars.coo','GoodStars.cooORIGINAL')   #Keeping BACKUP....
        shutil.copy('Source.coo','Source.cooORIGINAL')   #Keeping BACKUP....
        print ('\n For taking coordinates of good sky. Press _x_ over blank sky areas.')
        imx = iraf.imexam(Stdout=1)
        with open('BlankSky.coo','w') as foo :    #Creating blank sky coords files
            i = 0
            while i < len(imx) :               
                foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                i = i+1

        print('\n Use the ds9 to Mark circles centered at locations to do qphot.')
        print('Enter the center X Y and radius of aperture for qphot annulus and dannulus for sky')
        print('Enter values space separated in the format below. Enter "q" to exit')
        print('X   Y   Aperture  Annulus  Dannulus ')
        with open('qphotinput.txt','w') as foo:    #Creating the qphot parameter file
            qphot_inp="junk"
            while (qphot_inp != "q") :
                qphot_inp=raw_input("|> ")
                boolvar=True
                for i in qphot_inp.split() : boolvar = boolvar and is_number(i)
                if boolvar and (len(qphot_inp.split()) == 5) : foo.write(qphot_inp+' \n')
                elif (qphot_inp != "q") : print("Wrong Entry. Please enter properly the 5 values. q is to stop.")
    else :
        print("Source.coo, GoodStars.coo, BlankSky.coo, qphotinput.txt already exists in "+os.path.join(PC.MOTHERDIR,PC.OUTDIR)+". If you need to recreate, rename/remove them before calling this step.")
    #Finished all first images data collection. Now going forward..

    print("\n All required human input of coordinates taken..")

    iraf.cd(backupPWD) # Going back to the directory from were we entered this function.


def Createlist_subrout(PC):
    """ Creates the Images4Photo.in containing the image name , filter, exposure time, threshold """
    fooOUT = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'Images4Photo.in'),'w')
    directories = LoadDirectories(PC,CONF=False)
    for night in directories:
        PC.currentnight = night # Upgating the night directory for using PC.GetFullPath()
        print('Working on night: '+night)
        try:
            #First we load a dictionary of raw images to their filters
            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects.List'),'r') as FiltrFILE:
                FiltrFILElist = list(FiltrFILE)
            Filtrfiledic = dict([(filtset.split()[0],shlex.split(filtset.rstrip())[1]) for filtset in FiltrFILElist])  #Dictionary of filterset for each image.
            Exptimefiledic = dict([(filtset.split()[0],shlex.split(filtset.rstrip())[2]) for filtset in FiltrFILElist])  #Dictionary of filterset for each image.
            Datefiledic = dict([(filtset.split()[0],shlex.split(filtset.rstrip())[3]) for filtset in FiltrFILElist])  #Dictionary of filterset for each image.
            Timefiledic = dict([(filtset.split()[0],shlex.split(filtset.rstrip())[4]) for filtset in FiltrFILElist])  #Dictionary of filterset for each image.
            #Now Read and write the images to do photometry one by one.
            if PC.IMGCOMBINE == 'Y':
                ImgsFILE = open(PC.GetFullPath('AllObjects-ProcessedCombinedImg.List'),'r')
            else:
                ImgsFILE = open(PC.GetFullPath('AllObjects-ProcessedImg.List'),'r')
        except IOError as e:
            print('Cannot open the image file list.')
            print(e)
            print('So skipping this directory.')
            print('-'*60)
            continue

        for imgline in ImgsFILE:
            img = imgline.rstrip().split()[1]
            imgKey = imgline.rstrip().split()[0]
            imgfilter = Filtrfiledic[imgKey]
            imgexptime = Exptimefiledic[imgKey]
            date = Datefiledic[imgKey]
            time = Timefiledic[imgKey]
            fooOUT.write('{0}  "{1}"  {2}  {3} {4} {5}\n'.format(PC.GetFullPath(img),imgfilter,imgexptime,PC.THRESHOLD,date,time))
        ImgsFILE.close()
    fooOUT.close()
    print('All nights over...')

def GroupifiedList(inplist):
    """ A generator, which returns a sublist of input string list each call. 
    Each subgroup is identified by a blank entry in list."""
    sublist = []
    for entry in inplist:
        if len(entry.rstrip().split()) == 0:
            if len(sublist) > 0:
                yield sublist
                sublist = []  # Reset the sublist
        else:
            sublist.append(entry) # Keep appending..
    #yield final sublist if any remains
    if len(sublist) > 0:
        yield sublist
        
        

def AlignNcombine_subrout(PC,method="average"):
    """ This will align and combine the images in each paragraph in /FirstoneANDcombinedImages.List file. For photometry images it is interactive and will ask too select few images in first frame of each filter set."""
    
    directories = LoadDirectories(PC,CONF=False)
    for night in directories:
        PC.currentnight = night # Upgating the night directory for using PC.GetFullPath()
        print('Working on night: '+night)
        if PC.TODO == 'S' :
            print("Going forward to do blind combine fo each spectroscopy set noninteractively...")
            outObjectFinalFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-ProcessedCombinedImg.List'),'w')
            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-ProcessedImg.List'),'r') as Obj2CombFILE :
                for imglist in GroupifiedList(Obj2CombFILE):
                    objimgKey = imglist[0].split()[0]
                    if len(imglist) == 1 : #Single image. nothing to align and combine
                        OutCombimg = imglist[0].split()[1]
                    elif len(imglist) > 1 :
                        OutCombimg = os.path.splitext(imglist[0].split()[1])+'_align'+method+'_'+os.path.splitext(imglist[-1].split()[1])+'.fits'
                        imcombineInputfname = PC.GetFullPath(os.path.splitext(OutCombimg)[0]+'.imcomblist')
                        with open(imcombineInputfname,'w') as imcombineFILE:
                            imcombineFILE.write('\n'.join([PC.GetFullPath(imgline.split()[1]) for imgline in imglist])+'\n')
                        ImgCombineWithZeroFloating(imcombineInputfname,PC.GetFullPath(OutCombimg),cmethod=method,czero="median",creject="sigclip",cstatsection=PC.FLATSTATSECTION)
                    outObjectFinalFILE.write('{0}  {1}\n'.format(objimgKey,OutCombimg))
            #Done! Close the file, and continue to next directory.....
            outObjectFinalFILE.close()
            continue

        #Remaining is done only for Photometry......####
        iraf.imalign.unlearn()

        #Load all the X,Y coords of star indexed for every file already
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects2Combine.List'),'r') as XYFILE :
            XYfiledic = dict([(XYset.split()[0],XYset.rstrip().split()[1:]) for XYset in XYFILE if len(XYset.split()) == 3 ])  #Dictionary of XY coords for each image.

        if len(XYfiledic) == 0 : #No images this night..
            print('No images to work on this night. skipping...')
            continue

        print("Going forward to do align and combine for each image set Interactively...")
        outObjectFinalFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-ProcessedCombinedImg.List'),'w')
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-ProcessedImg.List'),'r') as Obj2CombFILE :
            for imglist in GroupifiedList(Obj2CombFILE):
                objimgKey = imglist[0].split()[0]
                if len(imglist) == 1 : #Single image. nothing to align and combine
                    OutCombimg = imglist[0].split()[1]
                elif len(imglist) > 1 :
                    OutCombimg = os.path.splitext(imglist[0].split()[1])[0]+'_align'+method+'_'+os.path.splitext(imglist[-1].split()[1])[0]+'.fits'
                    

                    OutCoofile = os.path.splitext(OutCombimg)[0]+'.GScoo'
                    RefimageKey = imglist[0].split()[0]
                    Refimage = imglist[0].split()[1]
                    Xref = float(XYfiledic[RefimageKey][0])
                    Yref = float(XYfiledic[RefimageKey][1])
                    iraf.display(PC.GetFullPath(Refimage),1) 
                    print ('Press _a_ over a few (~ 4 or 5) good stars to align, u can press s also, but DONT press r \n')
                    imx = iraf.imexam(Stdout=1)
                    with open(PC.GetFullPath(OutCoofile),'w') as foo :
                        i = 2
                        while i < len(imx) :               
                            foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                            i = i+2

                    #Now enter the crude shifts for other images from our dic in the file. And also create text files containing images to align and aligned output
                    with open(PC.GetFullPath('shifts.in'),'w') as foo2 :
                        alignInpfname = os.path.splitext(OutCombimg)[0]+'.ditherList'
                        alignOutfname = os.path.splitext(OutCombimg)[0]+'.AlignedditherList'
                        imcombineInputfname = os.path.splitext(OutCombimg)[0]+'.imcombineInputList'
                        imgs2align = open(PC.GetFullPath(alignInpfname),'w')    #Once the world migrates to Python 2.7+, these files also should be opened in the same with command above...
                        imgs2alignOUT = open(PC.GetFullPath(alignOutfname),'w')
                        imcombineInputFILE = open(PC.GetFullPath(imcombineInputfname),'w')
                        imcombineInputFILE.write(PC.GetFullPath(Refimage)+'\n') #First reference image
                        for imgline in imglist[1:]:
                            imgKey = imgline.split()[0]
                            img = imgline.split()[1]
                            Xin = float(XYfiledic[imgKey][0])  
                            Yin = float(XYfiledic[imgKey][1])
                            foo2.write(str(Xref-Xin)+'   '+str(Yref-Yin)+'\n')
                            imgs2align.write(PC.GetFullPath(img)+'\n')
                            imgs2alignOUT.write(PC.GetFullPath('s'+img)+'\n')
                            imcombineInputFILE.write(PC.GetFullPath('s'+img)+'\n')

                    imgs2align.close()
                    imgs2alignOUT.close()
                    imcombineInputFILE.close()
                    try :  #Now align and if succeeded combine those images....
                        iraf.imalign(input='@'+PC.GetFullPath(alignInpfname), reference=PC.GetFullPath(Refimage), coords=PC.GetFullPath(OutCoofile), output='@'+PC.GetFullPath(alignOutfname), shifts=PC.GetFullPath('shifts.in'), interp_type="nearest",boundary_type="constant",trimimages="no")
                        ImgCombineWithZeroFloating(PC.GetFullPath(imcombineInputfname),PC.GetFullPath(OutCombimg),cmethod=method,czero="median",creject="sigclip",cstatsection=PC.FLATSTATSECTION)
                    except iraf.IrafError as e :
                        print ('IRAF ERROR : Some image might be having problem. Remove it and try later')
                        print (e)
                        print('-'*60)
                        traceback.print_exc(file=sys.stdout)
                        print('-'*60)
            
                outObjectFinalFILE.write('{0}  {1}\n'.format(objimgKey,OutCombimg))
        outObjectFinalFILE.close()
    print('All nights over...')             


def DivideSmoothGradient(PC,inputimg,outputimg):
    """ This will divide the smooth gradient in image based on median filter smooting parameters in PC.
        In the end it will return the output filename."""
    hdulist = fits.open(inputimg)
    inputimgdata = hdulist[0].data
    print('Calculating median filtered gradient background using size : {0}'.format(PC.DVDMEDSMOOTHSIZE))
    try:
        smoothGrad = filters.median_filter(inputimgdata,size=PC.DVDMEDSMOOTHSIZE)
    except MemoryError:
        print('*** MEMORY ERROR : Skipping median filter Division ***')
        print('Try giving a smaller smooth size for median filter insted of {0}'.format(PC.DVDMEDSMOOTHSIZE))
        outputimg = inputimg  # Returning mack the input file name since no subtraction was done.
    else:
        prihdr = hdulist[0].header
        hdulist[0].data = inputimgdata / smoothGrad
        prihdr.add_history('Divided median filter Size:{0}'.format(PC.DVDMEDSMOOTHSIZE))
        hdulist.writeto(outputimg)
    finally:
        hdulist.close()
    # Return the name of the output filename
    return outputimg
                                
def SubtractSmoothGradient(PC,inputimg,outputimg):
    """ This will subract smooth gradients in image based on median filter smooting parameters in PC.
        In the end it will return the output filename."""
    hdulist = fits.open(inputimg)
    inputimgdata = hdulist[0].data
    print('Calculating median filtered gradient background using size : {0}'.format(PC.MEDSMOOTHSIZE))
    try:
        smoothGrad = filters.median_filter(inputimgdata,size=PC.MEDSMOOTHSIZE)
    except MemoryError:
        print('*** MEMORY ERROR : Skipping median filter Subtraction ***')
        print('Try giving a smaller smooth size for median filter')
        outputimg = inputimg  # Returning mack the input file name since no subtraction was done.
    else:
        prihdr= hdulist[0].header
        hdulist[0].data = inputimgdata - smoothGrad
        prihdr.add_history('Subtracted median filter Size:{0}'.format(PC.MEDSMOOTHSIZE))
        hdulist.writeto(outputimg)
    finally:
        hdulist.close()
    # Return the name of the output filename
    return outputimg
                
def FixBadPixels(PC,images,nightdir):
    """ This will run iraf task proto.fixpix to interpolate badpixels """
    if PC.TODO=='P' : PixelMask=nightdir+'/'+PC.PhotBadPixelMaskName
    elif PC.TODO=='S' : PixelMask=nightdir+'/'+PC.SpecBadPixelMaskName
    else : 
        print("What are you doing? (S or P)")
        return
    if not os.path.isfile(PixelMask):
        print("No Bad Pixel Mask file found by the name "+ PixelMask)
        print("Hence skipping Bad pixel interpolation")
        return

    iraf.proto(_doprint=0)
    iraf.fixpix.unlearn()
    iraf.fixpix(images=images,masks=PixelMask)

def ImgCombineWithZeroFloating(imglistfname,outputfile,cmethod="median",czero="median",creject="avgsigclip",cstatsection='[200:800,200:800]'):
    """ Returns the combined image with actuall average median flux, It does zero scaleing only for sigma rejection of stars. This is needed to remove faint stars in rejection algorithm when the background sky itself is varying from frame to frame. """
    iraf.imcombine.unlearn()
    Xmin=float(cstatsection[1:-1].split(',')[0].split(':')[0])  #Everything now in fits coordinates
    Xmax=float(cstatsection[1:-1].split(',')[0].split(':')[1])
    Ymin=float(cstatsection[1:-1].split(',')[1].split(':')[0])
    Ymax=float(cstatsection[1:-1].split(',')[1].split(':')[1])

    if czero == "median" : statfunction = np.median
    elif czero == "average" : statfunction = np.mean
    else : 
        print('Error: czero should be median or average. Unknown option {0}'.format(czero))
        raise

    with open(imglistfname,'r') as imgfile:
        statlist=[]
        for img in imgfile:
            img = img.rstrip()
            statlist.append(statfunction(fits.getdata(img)[Ymin-1:Ymax,Xmin-1:Xmax]))
    print('{0} of images: {1}'.format(czero,str(statlist)))
    statAvg=np.mean(statlist)
    Zeroshifts= statAvg - np.array(statlist)
    print('Zeroshifts of images: {0} :: ImgAvg ={1}'.format(str(Zeroshifts),statAvg))
    with open(outputfile+'_zeroshifts.txt','w') as zeroshiftFILE:
        for shift in Zeroshifts: 
            zeroshiftFILE.write('{0} \n'.format(shift))
    # Now call iraf imcombine with zero scaling
    iraf.imcombine(input='@'+imglistfname, output=outputfile, combine=cmethod, reject=creject, statsec=cstatsection, zero='@'+outputfile+'_zeroshifts.txt')
    
def Flat_basicCorrections_subrout(PC):
    """ This will create corresponding normalized flats and divide for flat correction. Also do mask pixel masking in end."""
    iraf.imcombine.unlearn()
    iraf.imstatistics.unlearn()
    directories = LoadDirectories(PC,CONF=False)

    for night in directories:
        PC.currentnight = night # Upgating the night directory for using GetFullPath()
        print('Working on night: '+night)
        
        #Load all the names of each object image and its Bias Subtracted version
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-BSImg.List'),'r') as ObjBSImgFILE :
            ObjBSimgdic = dict([(shlex.split(objimg)[0],shlex.split(objimg)[1]) for objimg in ObjBSImgFILE])  #Dictionary of objeckKey to BSimg

        #Load all the FilterSet indexing file data
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects.List'),'r') as FiltrFILE :
            FiltrFILElist = list(FiltrFILE)
        Filtrfiledic = dict([(shlex.split(filtset)[0],shlex.split(filtset.rstrip())[1]) for filtset in FiltrFILElist])  #Dictionary of filterset for each image.
        Exptimefiledic = dict([(shlex.split(filtset)[0],shlex.split(filtset.rstrip())[2]) for filtset in FiltrFILElist])  #Dictionary of exptime for each image.
        
        MasterFlatDic = dict()
        ObjectFlatDic = dict()
        ObjectFinalImgDic = {}  ## To store the FINAL Image corresponding to each object frame

        for objimgKey,objimg in ObjBSimgdic.items():
            OutMasterFlat = os.path.splitext(objimg)[0]+'_flat.fits'
            InpFlatFileList = os.path.splitext(objimgKey)[0]+'_Flat.list'
            #Load the flatnames to create a unique flatkey for each combination of flats
            with open(PC.GetFullPath(InpFlatFileList),'r') as flatlistFILE:
                listofflats = [flatimg.rstrip() for flatimg in flatlistFILE]
            if len(listofflats) != 0 : # Atlest one flat exist for this object
                FlatKey = ''.join(sorted(listofflats))  # a Key to uniquely identify the flat combination.
                try:  # symlink if the flat is already created
                    os.symlink(PC.GetFullPath(MasterFlatDic[FlatKey]),PC.GetFullPath(OutMasterFlat))
                except KeyError:
                    # Create the Master Flat and add to the dictionary
                    # First combine these Bias subtracted flats to  unnormalised flat
                    outflatname = os.path.splitext(OutMasterFlat)[0]+'_unNorm.fits'
                    if len(listofflats) == 1:
                        #Nothing to combine, simply creat a simbilic link to it.
                        os.symlink(listofflats[0],PC.GetFullPath(outflatname))
                    else:
                        print('Image section used for statistics of Flat is '+PC.FLATSTATSECTION)
                        iraf.imcombine(input='@'+PC.GetFullPath(InpFlatFileList), output=PC.GetFullPath(outflatname), combine="median", scale="median",reject="sigclip", statsec=PC.FLATSTATSECTION)

                    if (PC.TODO == 'S') and (PC.CONTINUUMGRADREMOVE == 'Y'):
                        # We will normalise this continuum flat using its median smoothed version
                        OutMasterFlat = DivideSmoothGradient(PC,PC.GetFullPath(outflatname),PC.GetFullPath(OutMasterFlat))
                    else:
                        #We will normalise this flat with the mode of pixels in FlatStatSection
                        statout = iraf.imstatistics(PC.GetFullPath(outflatname)+PC.FLATSTATSECTION,fields='mode',Stdout=1)
                        mode = float(statout[1])
                        iraf.imarith(operand1=PC.GetFullPath(outflatname),op="/",operand2=mode,result=PC.GetFullPath(OutMasterFlat))
                    # Add to the dictionary
                    MasterFlatDic[FlatKey] = OutMasterFlat
                finally:
                    ObjectFlatDic[objimgKey] = OutMasterFlat
            else:
                ObjectFlatDic[objimgKey] = None


            ###### Apply Flat correction to object frames
            if ObjectFlatDic[objimgKey] is not None:  # Flat exists
                print('Flat correcting '+objimg)
                OutFCobjimg = os.path.splitext(objimg)[0]+'_FC.fits'
                iraf.imarith(operand1=PC.GetFullPath(objimg),op="/",operand2=PC.GetFullPath(ObjectFlatDic[objimgKey]),result=PC.GetFullPath(OutFCobjimg))
                # Update the name of output file.
                OutFinalimg = OutFCobjimg 
            else:
                OutFinalimg = objimg
                
            ###### Now interpolate the bad pixels in the final image.
            FixBadPixels(PC,PC.GetFullPath(OutFinalimg),night)
            
            #If asked to do a smooth gradient removal in image, do it after everything is over now..
            if PC.GRADREMOVE == 'Y':
                print('Removing smooth Gradient from '+objimg)
                OutGSobjimg = os.path.splitext(OutFinalimg)[0]+'_GS.fits'
                #If Filter subtraciton was sucessfull it will return the output filename, else inputname
                OutGSobjimg = SubtractSmoothGradient(PC,PC.GetFullPath(OutFCobjimg),PC.GetFullPath(OutGSobjimg))
                #Updating the final image name with new _GS appended 
                OutFinalimg = os.path.basename(OutGSobjimg)
            
            ObjectFinalImgDic[objimgKey] = OutFinalimg



        # Now write all the output files required for the next step.
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects2Combine.List'),'r') as Obj2CombFILE :
            NewFilter = 'Blah'
            NewExptime = '-999'
            OutObjectFinalfilename = 'AllObjects-ProcessedImg.List'
            outObjectFinalFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,OutObjectFinalfilename),'w')
            for imgline in Obj2CombFILE:
                if len(imgline.split()) == 0: #A blank line 
                    if PC.TODO == 'S':  # If it is spectroscopy, this gap means different dither/change/exptime in grism
                        outObjectFinalFILE.write('\n')
                    elif PC.TODO == 'P':  # Skip only if there is a change in filter or exptime in next step
                        pass # We will deal filter change in next step.
                    continue
                #Write Objectframe and the Final processed object frame
                objimgKey = imgline.rstrip().split()[0]
                #Make sure to enter balnk line if there is a change in filter or change in exptime
                OldFilter = NewFilter
                NewFilter = Filtrfiledic[objimgKey]
                OldExptime = NewExptime
                NewExptime = Exptimefiledic[objimgKey]
                if (PC.TODO == 'P') and ((OldFilter != NewFilter) or (float(OldExptime) != float(NewExptime))) : 
                    outObjectFinalFILE.write('\n')
                # Write the output file
                outObjectFinalFILE.write('{0}  {1}\n'.format(objimgKey,ObjectFinalImgDic[objimgKey]))
                
            outObjectFinalFILE.close()

        if PC.IMGCOMBINE == 'Y' : print('Edit the spaces (if required) between image sets in file '+os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,OutObjectFinalfilename)+' to align and combine them in next step.')
    print('All nights over...')             



def Bias_Subtraction_subrout(PC,method="median"):
    """ This will combine (default=median) with avsigclip the bias for each image and also create corresponding file list for flats of each image. """
    iraf.imcombine.unlearn()
    iraf.imred(_doprint=0)
    iraf.bias(_doprint=0)
    iraf.colbias.unlearn()  

    iraf.ccdred(_doprint=0)
    iraf.zerocombine.unlearn()  

    directories = LoadDirectories(PC,CONF=False)

    if PC.USEALLFLATS == 'Y':
        SuperMasterFilterFlatdic = dict() # Dictionary to store all the final flats for using together form all night
        SuperMasterFlatListsFileNamesnFiltr = dict() # Dictionary to store all the final flat list to combine and its filter

    for night in directories:
        PC.currentnight = night # Upgating the night directory for using GetFullPath()
        print('Working on night: '+night)
        #Load all the Bias indexing file data
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-FinalBias.List'),'r') as BiasFILE :
            Biasfiledic = dict([(shlex.split(biasset)[0],biasset.rstrip().split()[1:]) for biasset in BiasFILE])  #Dictionary of bias list for each image.
        #Load all the Flat indexing file data
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-FinalFlat.List'),'r') as FlatFILE :
            Flatfiledic = dict([(shlex.split(flatset)[0],flatset.rstrip().split()[1:]) for flatset in FlatFILE])  #Dictionary of flats list for each image.
        if PC.TODO == 'S':  #Load all the Lamp indexing file data
            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-FinalLamp.List'),'r') as LampFILE :
                Lampfiledic = dict([(shlex.split(lampset)[0],lampset.rstrip().split()[1:]) for lampset in LampFILE])  #Dictionary of Lamps for each image.

        if PC.SEPARATESKY == 'Y' :
            #Load all the Sky files indexing file data
            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-FinalSky.List'),'r') as SkyFILE :
                Skyfiledic = dict([(shlex.split(skyset)[0],skyset.rstrip().split()[1:]) for skyset in SkyFILE])  #Dictionary of Sky list for each image.
        if PC.USEALLFLATS == 'Y':
            #Load all the Bias indexing file of each flat image
            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllFlat-FinalBias.List'),'r') as flatBiasFILE :
                FlatBiasfiledic = dict([(shlex.split(biasset)[0],biasset.rstrip().split()[1:]) for biasset in flatBiasFILE])  #Dictionary of bias list for each flat.
            

        #Load all the FilterSet indexing file data
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects.List'),'r') as FiltrFILE :
            FiltrFILElist = list(FiltrFILE)
        Filtrfiledic = dict([(shlex.split(filtset)[0],shlex.split(filtset.rstrip())[1]) for filtset in FiltrFILElist])  #Dictionary of filterset for each image.
        Exptimefiledic = dict([(shlex.split(filtset)[0],shlex.split(filtset.rstrip())[2]) for filtset in FiltrFILElist])  #Dictionary of exptime for each image.

        # Also the filter of all flats if we have a seperate master flat plan
        if PC.USEALLFLATS == 'Y':
            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllFilter-FinalFlat.List'),'r') as FiltrFlatFILE :
                for filtrline in FiltrFlatFILE:
                    filtrline = filtrline.rstrip()
                    filtr = shlex.split(filtrline)[0]
                    for flatimg in shlex.split(filtrline)[1:]:
                        Filtrfiledic[flatimg] = filtr
            
            
        if PC.OVERSCAN != 'N': #Overscan Subtracting to be done
            print('Doing Overscan subtraction from: {0}'.format(PC.OVERSCAN))
            OSsuffix = '_OS.fits'
            BiasfiledicOS = {}
            FlatfiledicOS = {}
            LampfiledicOS = {}
            SkyfiledicOS = {}
            FlatBiasfiledicOS = {} # Will be used only if PC.USEALLFLATS == 'Y'
            for objimg in Biasfiledic:
                print('DEBUG: '+objimg)
                OSobjimg = PC.GetFullPath(os.path.splitext(objimg)[0]+OSsuffix)
                iraf.colbias(input=PC.GetFullPath(objimg),output=OSobjimg, bias=PC.OVERSCAN, interactive='no') # Main object image
                BiasfiledicOS[objimg] = []
                for biasimg in Biasfiledic[objimg]:
                    outputimg = PC.GetFullPath(os.path.splitext(biasimg)[0]+OSsuffix)
                    BiasfiledicOS[objimg].append(os.path.basename(outputimg))
                    if not os.path.isfile(outputimg): # If the overscan subtracted image already doesn't exist
                        iraf.colbias(input=PC.GetFullPath(biasimg),output=outputimg, bias=PC.OVERSCAN, interactive='no') # Bias  image
                FlatfiledicOS[objimg] = []
                for flatimg in Flatfiledic[objimg]:
                    # no need to specifically ask it to take only basename for outputimg 
                    outputimg = PC.GetFullPath(os.path.splitext(os.path.basename(flatimg))[0]+OSsuffix) 
                    FlatfiledicOS[objimg].append(os.path.basename(outputimg))
                    if not os.path.isfile(outputimg): # If the overscan subtracted image already doesn't exist
                        iraf.colbias(input=PC.GetFullPath(flatimg),output=outputimg, bias=PC.OVERSCAN, interactive='no') # Flat image
                        
                if PC.TODO == 'S' :
                    LampfiledicOS[objimg] = []
                    for lampimg in Lampfiledic[objimg]:
                        outputimg = PC.GetFullPath(os.path.splitext(lampimg)[0]+OSsuffix)
                        LampfiledicOS[objimg].append(os.path.basename(outputimg))
                        if not os.path.isfile(outputimg): # If the overscan subtracted image already doesn't exist
                            iraf.colbias(input=PC.GetFullPath(lampimg),output=outputimg, bias=PC.OVERSCAN, interactive='no') # Lamp image

                if PC.SEPARATESKY == 'Y' :
                    SkyfiledicOS[objimg] = []
                    for skyimg in Skyfiledic[objimg]:
                        outputimg = PC.GetFullPath(os.path.splitext(skyimg)[0]+OSsuffix)
                        SkyfiledicOS[objimg].append(os.path.basename(outputimg))
                        if not os.path.isfile(outputimg): # If the overscan subtracted image already doesn't exist
                            iraf.colbias(input=PC.GetFullPath(skyimg),output=outputimg, bias=PC.OVERSCAN, interactive='no') # sky image

            if PC.USEALLFLATS == 'Y':
                for flatimg in FlatBiasfiledic:
                    OSflatimg = PC.GetFullPath(os.path.splitext(flatimg)[0]+OSsuffix)
                    if not os.path.isfile(OSflatimg): # If the overscan subtracted image already doesn't exist
                        iraf.colbias(input=PC.GetFullPath(flatimg),output=OSflatimg, bias=PC.OVERSCAN, interactive='no') # Main object image
                    FlatBiasfiledicOS[flatimg] = []
                    for biasimg in FlatBiasfiledic[flatimg]:
                        outputimg = PC.GetFullPath(os.path.splitext(biasimg)[0]+OSsuffix)
                        FlatBiasfiledicOS[flatimg].append(os.path.basename(outputimg))
                        if not os.path.isfile(outputimg): # If the overscan subtracted bias already doesn't exist
                            iraf.colbias(input=PC.GetFullPath(biasimg),output=outputimg, bias=PC.OVERSCAN, interactive='no') # Bias  image
                    

            # Finnaly update all the filename dictionary lists
            Biasfiledic = BiasfiledicOS 
            Flatfiledic = FlatfiledicOS 
            Lampfiledic = LampfiledicOS 
            Skyfiledic =  SkyfiledicOS 
            FlatBiasfiledic = FlatBiasfiledicOS # Not empty Only if PC.USEALLFLATS == 'Y'

        MasterBiasDic = {}  # To temperorly store created master frames
        ObjectBiasDic = {}  # To store the BIAS corresponding to each object frame
        ObjectLampDic = {}  # To store the LAMP corresponding to each object frame
        MasterSkyDic = {}   # To temperorly store created sky frames             
        ObjectSkyDic = {}   # To store the SKY corresponding to each object frame
        ObjectFinalImgDic = {}  ## To store the FINAL Image corresponding to each object frame

        if PC.USEALLFLATS == 'Y': # Subtract bias from all flats first
            for flatimgKey in FlatBiasfiledic:
                flatimg = flatimgKey if (PC.OVERSCAN == 'N') else os.path.splitext(flatimgKey)[0]+OSsuffix
                ####### Combine Bias frames.
                print('Combining Bias frames for '+flatimg)
                OutMasterBias = os.path.splitext(flatimg)[0]+'_Bias.fits'
                BiasKey = ''.join(sorted(FlatBiasfiledic[flatimgKey]))  # a Key to uniquely identify the bias combination.
                try:  # symlink if the bias is already created
                    os.symlink(PC.GetFullPath(MasterBiasDic[BiasKey]),PC.GetFullPath(OutMasterBias))
                except KeyError:
                    # Create the Master bias and add to the dictionary
                    OutMasterBiasList = os.path.splitext(flatimg)[0]+'_Bias.list'
                    with open(PC.GetFullPath(OutMasterBiasList),'w') as biaslistFILE:
                        biaslistFILE.write('\n'.join([PC.GetFullPath(biasimg) for biasimg in FlatBiasfiledic[flatimgKey]])+'\n')
                    iraf.zerocombine(input= "@"+PC.GetFullPath(OutMasterBiasList), output=PC.GetFullPath(OutMasterBias), combine="median", ccdtype="")
                    # Add to the dictionary
                    MasterBiasDic[BiasKey] = OutMasterBias
                finally:
                    # Subtract bias form this flat image
                    BSflatimg = os.path.splitext(os.path.basename(flatimg))[0]+'_BS.fits'
                    if not os.path.isfile(PC.GetFullPath(BSflatimg)): #If the bias subtracted flat already doesn't exist.
                        iraf.imarith(operand1=PC.GetFullPath(flatimg),op="-",operand2=PC.GetFullPath(OutMasterBias),result=PC.GetFullPath(BSflatimg))
                    
                    SuperMasterFilterFlatdic.setdefault(Filtrfiledic[flatimgKey], []).append(PC.GetFullPath(BSflatimg))


        for objimgKey in Biasfiledic:
            objimg = objimgKey if (PC.OVERSCAN == 'N') else os.path.splitext(objimgKey)[0]+OSsuffix
            ####### Combine Bias frames.
            print('Combining Bias frames for '+objimg)
            OutMasterBias = os.path.splitext(objimg)[0]+'_Bias.fits'
            BiasKey = ''.join(sorted(Biasfiledic[objimgKey]))  # a Key to uniquely identify the bias combination.
            try:  # symlink if the bias is already created
                os.symlink(PC.GetFullPath(MasterBiasDic[BiasKey]),PC.GetFullPath(OutMasterBias))
            except KeyError:
                # Create the Master bias and add to the dictionary
                OutMasterBiasList = os.path.splitext(objimg)[0]+'_Bias.list'
                with open(PC.GetFullPath(OutMasterBiasList),'w') as biaslistFILE:
                    biaslistFILE.write('\n'.join([PC.GetFullPath(biasimg) for biasimg in Biasfiledic[objimgKey]])+'\n')
                iraf.zerocombine(input= "@"+PC.GetFullPath(OutMasterBiasList), output=PC.GetFullPath(OutMasterBias), combine="median", ccdtype="")
                # Add to the dictionary
                MasterBiasDic[BiasKey] = OutMasterBias
            finally:
                ObjectBiasDic[objimgKey] = OutMasterBias

            ####### Subtract Bias frames from Lamp frames
            if PC.TODO == 'S':
                BSlampimgs = [os.path.splitext(os.path.basename(lampimg))[0]+'_BS.fits'  for lampimg in Lampfiledic[objimgKey]]
                for lampimg,bslampimg in zip(Lampfiledic[objimgKey],BSlampimgs):
                    if not os.path.isfile(PC.GetFullPath(bslampimg)): #If the bias subtracted lamp already doesn't exist.
                        iraf.imarith(operand1=PC.GetFullPath(lampimg),op="-",operand2=PC.GetFullPath(ObjectBiasDic[objimgKey]),result=PC.GetFullPath(bslampimg))
                # Keep only the first lamp in the list, since we need only one for Wavelength calibration
                ObjectLampDic[objimgKey] = BSlampimgs[0]


            ###### Subtract Bias from Flat frames for each image
            # Also Create the list of flat images to combine for each image
            OutMasterFlatList = os.path.splitext(objimgKey)[0]+'_Flat.list'

            if Flatfiledic[objimgKey]:  # If flats exist!
                print('Subtracting Bias from flats for '+objimg)
                BSflatimgs = [os.path.splitext(os.path.basename(flatimg))[0]+'_BS.fits'  for flatimg in Flatfiledic[objimgKey]]
                for flatimg,bsflatimg in zip(Flatfiledic[objimgKey],BSflatimgs):
                    if not os.path.isfile(PC.GetFullPath(bsflatimg)): #If the bias subtracted flat already doesn't exist.
                        iraf.imarith(operand1=PC.GetFullPath(flatimg),op="-",operand2=PC.GetFullPath(ObjectBiasDic[objimgKey]),result=PC.GetFullPath(bsflatimg))

                if PC.USEALLFLATS != 'Y': # Create the list of flats to combine form this directory itself.
                    with open(PC.GetFullPath(OutMasterFlatList),'w') as flatlistFILE:
                        flatlistFILE.write('\n'.join([PC.GetFullPath(bsflatimgs) for bsflatimgs in BSflatimgs])+'\n')
            else:
                print('ALERT: No Flats for doing flat correction of {0} from same night'.format(objimg))
            if PC.USEALLFLATS == 'Y': # Save the name of the list to create flat from all other directories in the end of this function
                SuperMasterFlatListsFileNamesnFiltr[PC.GetFullPath(OutMasterFlatList)] = Filtrfiledic[objimgKey]

            ####### Bias/Sky subtract object frame
            BiasRemovedObjimg = os.path.splitext(objimg)[0]+'_BS.fits'
            if PC.SEPARATESKY=='Y': 
                print('Subtracting sky from '+objimg)
                ####### Combine Sky frames and subtract it instead of bias
                OutMasterSky = os.path.splitext(objimg)[0]+'_Sky.fits'
                SkyKey = ''.join(sorted(Skyfiledic[objimgKey]))  # a Key to uniquely identify the bias combination.
                try:  # symlink if the bias is already created
                    os.symlink(PC.GetFullPath(MasterSkyDic[SkyKey]),PC.GetFullPath(OutMasterSky))
                except KeyError:
                    # Create the Master sky and add to the dictionary
                    OutMasterSkyList = os.path.splitext(objimg)[0]+'_Sky.list'
                    with open(PC.GetFullPath(OutMasterSkyList),'w') as skylistFILE:
                        skylistFILE.write('\n'.join([PC.GetFullPath(skyimg) for skyimg in Skyfiledic[objimgKey]])+'\n')
                    ImgCombineWithZeroFloating(PC.GetFullPath(OutMasterSkyList),PC.GetFullPath(OutMasterSky),cmethod="median",czero="median",creject="pclip",cstatsection=PC.FLATSTATSECTION)
                    # Add to the dictionary
                    MasterSkyDic[SkyKey] = OutMasterSky
                finally:
                    ObjectSkyDic[objimgKey] = OutMasterSky

                #Now subtract the sky form the science frame
                iraf.imarith(operand1=PC.GetFullPath(objimg),op="-",operand2=PC.GetFullPath(ObjectSkyDic[objimgKey]),result=PC.GetFullPath(BiasRemovedObjimg))
                                             
            else:  # Subtract the Bias directly
                print('Subtracting Bias from '+objimg)
                iraf.imarith(operand1=PC.GetFullPath(objimg),op="-",operand2=PC.GetFullPath(ObjectBiasDic[objimgKey]),result=PC.GetFullPath(BiasRemovedObjimg))

            ObjectFinalImgDic[objimgKey] = BiasRemovedObjimg

        # Now write all the output files required for the next step.
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects2Combine.List'),'r') as Obj2CombFILE :
            OutObjectFinalfilename = 'AllObjects-BSImg.List'
            outObjectFinalFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,OutObjectFinalfilename),'w')
            if PC.TODO == 'S':
                outObjectLampFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-BSFinalLamp.List'),'w')
            for imgline in Obj2CombFILE:
                if len(imgline.split()) == 0: #A blank line 
                    continue
                #Write Objectframe and the Final processed object frame
                objimgKey = imgline.rstrip().split()[0]
                outObjectFinalFILE.write('{0}  {1}\n'.format(objimgKey,ObjectFinalImgDic[objimgKey]))
                if PC.TODO == 'S': # Also create the log of image and its new Lamp
                    outObjectLampFILE.write('{0}  {1}\n'.format(objimgKey,ObjectLampDic[objimgKey]))
                
            outObjectFinalFILE.close()
            if PC.TODO == 'S': outObjectLampFILE.close()

    ### Now write the Flats to cobine list if we are using SuperMasterFlat
    if PC.USEALLFLATS == 'Y':
        for flatslistfilename,filtr in SuperMasterFlatListsFileNamesnFiltr.items():
            with open(flatslistfilename,'w') as flatlistFILE:
                flatlistFILE.write('\n'.join(SuperMasterFilterFlatdic[filtr])+'\n')

    print('All nights over...')             
                

def Manual_InspectCalframes_subrout(PC):
    """ This will display Flats, Sky and Lamps one after other, and based on user input select/reject """
    directories = LoadDirectories(PC,CONF=True)
    filelist = ['AllObjects-Bias.List','AllObjects-Flat.List']
    outfilelist = ['AllObjects-FinalBias.List','AllObjects-FinalFlat.List']
    if PC.USEALLFLATS == 'Y':
        filelist += ['AllFilter-Flat.List', 'AllFlat-Bias.List']
        outfilelist +=['AllFilter-FinalFlat.List', 'AllFlat-FinalBias.List']
    if PC.TODO == 'S' :
        filelist.append('AllObjects-Lamp.List')
        outfilelist.append('AllObjects-FinalLamp.List')
    if PC.SEPARATESKY == 'Y' :
        filelist.append('AllObjects-Sky.List')
        outfilelist.append('AllObjects-FinalSky.List')

    AcceptAllThisNight = False  #Flags to skip this step
    AcceptAllEveryNight = False
    for night in directories:
        print("Working on night :\033[91m {0} \033[0m ".format(night))
        AlwaysRemove = []
        AlwaysAccept = []
        for inpfile,outfile in zip(filelist,outfilelist):
            print('-*-'*8)
            print('Files in: \033[91m {0} \033[0m'.format(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,inpfile)))
            print('-*-'*8)
            inFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,inpfile),'r')
            ouFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,outfile),'w')
            print('-'*5)
            print('To skip all the remaining verifications, you can enter following two options')
            print('acceptall       # Accept all remaining images in this night')
            print('acceptallnights # Accept all images in all remaining nights')
            print('Use the above two options only when you are sure all the images are good. Do not use them if you are not sure.')
            print('-'*5)
            for inpline in inFILE:
                inplinelist = shlex.split(inpline.rstrip())
                if len(inplinelist) > 0 : ScienceImg = inplinelist[0]
                else : continue   #Skipping to next line
                CalImgs = [imgs for imgs in inplinelist[1:] if imgs not in AlwaysRemove]
                FinalCalImgs = CalImgs[:]
                print('For the science image: '+ScienceImg)
                if not AcceptAllThisNight :
                    for img in CalImgs:
                        if img not in AlwaysAccept:
                            iraf.display(night+'/'+img,1)
                            print(night+'/'+img)
                            verdict=''
                            verdict=raw_input('Enter "r" to reject, "ra" to reject always in future, "aa" to always accept in future:')
                            if verdict == 'r' :
                                FinalCalImgs.remove(img)
                                print("Removing this image : "+img)
                            elif verdict == 'ra' :
                                FinalCalImgs.remove(img)
                                AlwaysRemove.append(img)
                                print("Removing this image forever: "+img)
                            elif verdict == 'aa' :
                                AlwaysAccept.append(img)
                                print("Always accept this image forever this night : "+img)
                            elif verdict == 'acceptall' :
                                AcceptAllThisNight = True
                                print("Accepting every single remainging images of this night (Dangerous). ")
                                break
                            elif verdict == 'acceptallnights' :
                                AcceptAllEveryNight = True
                                AcceptAllThisNight = True
                                print("Accepting all remainging images from all remaining nights (Super Dangerous). ")
                                break
                if not FinalCalImgs : print('\033[91m ALERT: \033[0m No Calibration Images for {0} {1}'.format(night,ScienceImg))
                #Writing the final surviving calibration files to output file
                ouFILE.write('"{0}" {1}\n'.format(ScienceImg,' '.join(FinalCalImgs)))
            ouFILE.close()
            inFILE.close()
        if not AcceptAllEveryNight : AcceptAllThisNight = False
    print('All nights over...') 
               

def Manual_InspectObj_subrout(PC):
    """ This will display one image after other, and based on user input classify images of each dither position """
    directories=LoadDirectories(PC,CONF=True)
    if PC.TODO == 'P': print("Press _a_ and then _q_ over one good central star for selecting image")
    if PC.TODO == 'S': print("Press _j_ and then _q_ over one good position on dispersed spectrum for selecting image \n IMP: Press j on some good part of star spectrum, not on the sky region around.")
    for night in directories:
        print("Working on night : "+night)
        ObjFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects.List'),'r')
        Obj2CombFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects2Combine.List'),'w')
        newX = 0
        newY = 0
        newsX = 0
        FWHM = 4.0  #Not important what number you give here...
        newfilter = 'Blah'
        newexptime = '-999'
        for objline in ObjFILE:
            try:
                img = objline.split()[0]
            except IndexError:
                print('Blank line in '+os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects.List'))
                continue
            try :
                iraf.display(night+'/'+img,1)
            except iraf.IrafError as e:
                # ds9 might not be open, hence open it and try again
                print('No ds9 seems to be open, I am opening a ds9..')
                subprocess.Popen('ds9')
                time.sleep(4)  # Give 4 seconds for ds9 to startup
                iraf.display(night+'/'+img,1)
                
            print(objline)
            if PC.TODO == 'P':print('\n To discard this image press _q_ without pressing _a_')
            if PC.TODO == 'S':print('\n To discard this image press _q_ without pressing _j_')
            try:
                imx = iraf.imexam(Stdout=1)
            except iraf.IrafError as e :
                print('IRAF ERROR : This image %s might be having problem, still choosing it'%(night+'/'+img))
                print(e)
                if len(imx) < 1 :imx=['center= %f  peak fwhm= %f bkg'%(newsX,FWHM)]  #A dummy entry..
                
            if len(imx) < 1 : #Then discard this image
                print('Discarding image :'+night+'/'+img)
                continue
            #Else, continue below
            if PC.TODO == 'P': #If we were doing photometry
                oldX = newX
                oldY = newY
                oldfilter = newfilter
                newfilter = shlex.split(objline)[1]
                oldexptime = newexptime
                newexptime = shlex.split(objline)[2]
                FWHM = float(imx[3].split()[-1])
                newX = float(imx[2].split()[0])
                newY = float(imx[2].split()[1])
                #Print blank enter in the output file if the star has shifted
                StarShifted = np.sqrt((newX-oldX)**2 +(newY-oldY)**2) > 1*FWHM
            elif PC.TODO == 'S' : #If doing spectroscopy
                oldsX = newsX
                oldfilter = newfilter
                newfilter = shlex.split(objline)[1]
                oldexptime = newexptime
                newexptime = shlex.split(objline)[2]
                s = imx[-1]
                FWHM = float(s[s.rfind('fwhm=')+5:s.rfind('bkg')])
                newsX = float(s[s.rfind('center=')+7:s.rfind('peak')])
                #Print blank enter in the output file if the star has shifted
                StarShifted = np.abs(newsX-oldsX) > 1*FWHM
                
            # or filter wheels or exptime have changed.
            FiltersChanged = newfilter != oldfilter 
            ExptimeChanged = float(newexptime) != float(oldexptime)
            if StarShifted or FiltersChanged or ExptimeChanged : Obj2CombFILE.write('\n')

            #Now, add this img name to dither image list
            Obj2CombFILE.write(img+' '+str(newX)+' '+str(newY)+'\n')
        Obj2CombFILE.close()
        ObjFILE.close()
        print('We have made the selected list of images in '+os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects2Combine.List')+' \n Add blank lines between file names to prevent them from median combining. \n Remove the blank line between file names, which you want to combine.')
        raw_input("Press Enter to continue...")
        subprocess.call(PC.TEXTEDITOR.split()+[os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects2Combine.List')])

    print('All nights over...')
        
             
def SelectionofFrames_subrout(PC):
    """ Selects the images to reduce and create tables of corresponding Flat, Sky and Lamps """
    directories = LoadDirectories(PC,CONF=True)
    LogFilename = PC.NIGHTLOGFILE
    FiltREdic = dict()
    LampREdic = dict()
    SkyREdic = dict()
    # Load Instrument
    Instrument = InstrumentObject(PC)
    
    # Filter/Grism Column index in log file
    if PC.TODO == 'P' :
        FiltColumn = 3
    elif PC.TODO == 'S' :
        FiltColumn = 4


    print('-'*10) 
    print('For Regular Expression rules See: http://docs.python.org/2/howto/regex.html#regex-howto')
    print('Some examples of typical input are shown below')
    print(' .*M31.*   is the regular expression to select all the objects lines which has "M31" in it.')
    print(' .*M31.*sky.*   is to select all the object lines which has both "M31" and then "sky" in it.')
    print('While selecting Lamp, Flat etc, you can give a range of filenumbers to uniquely choose specific files')
    print(' .*Flat.* 7 18   is to select all filenames of _same_ filter which has "Continu" in it and has a filenumber in the range of 7 to 18')
    print('-'*10)
    ObjRE = ' '
    #Generating list of objects frames
    for night in directories:
        print("Working on night : "+night)
        print("Obs log file: file://{0}".format(os.path.join(PC.MOTHERDIR,night,LogFilename)))
        InpObjRE = raw_input("Enter Regular Expression to select Science object frames (default: {0}): ".format(ObjRE)).strip(' ')
        if InpObjRE:
            ObjRE = InpObjRE
        regexpObj = re.compile(r''+ObjRE)

        with open(os.path.join(night,LogFilename),'r') as imglogFILE :
            # Skip blank lines and Commented out lines with #
            imglogFILElines = [imageLINE.rstrip() for imageLINE in imglogFILE if ((imageLINE.strip() is not '') and (imageLINE[0] !='#'))]

        ObjListT = [imgline for imgline in imglogFILElines if regexpObj.search(' '.join(shlex.split(imgline)[0:2])) is not None ] # Search both Filename and OBJECT Name
        ObjList = []
        # Remove any obvious non Object files which creeps in through regexp
        for Objline in ObjListT:
            ImgFrame = Instrument.IdentifyFrame(Objline)
            if ImgFrame not in ['BIAS','LAMP_SPEC','FLAT_SPEC']:  # Just in case any of these creaped in ObJlist
                if not (PC.TODO=='P' and ('SPEC' in ImgFrame)) and not (PC.TODO=='S' and ('SPEC' not in ImgFrame)) :
                    ObjList.append(Objline)

        FiltList = set()  #Set to store the list of filters needed to find flat/Lamp for
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects.List'),'w') as ObjFILE :
            for Objline in ObjList:
                Name = shlex.split(Objline)[0]
                FiltOrGrism = shlex.split(Objline)[FiltColumn]
                Exptime = shlex.split(Objline)[2]
                Date = shlex.split(Objline)[7]
                Time = shlex.split(Objline)[8]
                FiltList.add(FiltOrGrism)
                ObjFILE.write('{0}   "{1}"  {2}  {3} {4}\n'.format(Name,FiltOrGrism,Exptime,Date,Time))

        # Create a complete list of filters is we need to use flats from all nights
        CompleteFiltList = set()
        if PC.USEALLFLATS == 'Y':
            for imgline in imglogFILElines:
                ImgFrame = Instrument.IdentifyFrame(imgline)
                if ImgFrame in ['BIAS','LAMP_SPEC']:  # Skip Bobvious Bias, Wavelength Calibration Lamps etc
                    continue    #Skip these and go to the next object.
                CompleteFiltList.add(shlex.split(imgline)[FiltColumn])
            

        if (not FiltList) and (not CompleteFiltList)  : #No files in this directory
            print('\033[91m ALERT: \033[0m No Images to reduce found in directory : {0}'.format(night))
            print('Please remove {0} from directory list next time.'.format(night))
            continue

        #Create dictionary of image sizes and available bias frames.
        BiasSizeDic = {}
        for imgline in imglogFILElines:
            if Instrument.IdentifyFrame(imgline) == 'BIAS':
                # Append to the list of biases stored in dictinary with the key (X,Y)
                BiasSizeDic.setdefault((imgline.split()[-3],imgline.split()[-2]), []).append(imgline.split()[0])

        #Now ask for flats in each filters
        Flatlistdic = dict()
        Flatsizedic = dict()
        print("Below in addition to regexp, if needed you can enter the starting and ending filenumbers separated by space also.")
        for filt in FiltList|CompleteFiltList:
            filenumbregexp=re.compile(r'.*')
            if filt not in FiltREdic.keys() : 
                if PC.TODO=='P': FiltREdic[filt] = '.*[Ff]lat.*'  #Setting default to *[Ff]lat*
                elif PC.TODO=='S': FiltREdic[filt] = '.*Halogen.*'
            #Ask user again to confirm or change if he/she needs to
            InpfiltRE = raw_input("Enter Regular Expression for the flat of filters %s (default: %s) : "%(str(filt),FiltREdic[filt])).strip(' ')
            if InpfiltRE :
                FiltREdic[filt] = InpfiltRE
                if len(InpfiltRE.split())==3 and is_number(InpfiltRE.split()[1]) and is_number(InpfiltRE.split()[2]):
                    filenumbregexp = re.compile('|'.join([str(i) for i in range(int(InpfiltRE.split()[1]), int(InpfiltRE.split()[2])+1)]))
                    FiltREdic[filt] = InpfiltRE.split()[0]
            regexpFilt = re.compile(r''+FiltREdic[filt])

            FlatList = []
            for imgline in imglogFILElines:
                ImgFrame = Instrument.IdentifyFrame(imgline)
                if ImgFrame in ['BIAS','LAMP_SPEC']:  # Skip Bobvious Bias, Wavelength Calibration Lamps etc
                    continue    #Skip these and go to the next object.
                if regexpFilt.search(' '.join(shlex.split(imgline)[0:2])) is not None :
                    if filt == shlex.split(imgline)[FiltColumn]: # Matching filter
                        if filenumbregexp.search(imgline.split()[-1]): # Last column is filenumber
                            FlatList.append(imgline.split()[0])
                            Flatsizedic[imgline.split()[0]] = tuple(shlex.split(imgline)[-3:-1])  # X, Y size of each flat
            Flatlistdic[filt] = FlatList  #Saving flat list for this filter set

        #Now if Separate sky is being used to subtract, ask for each filter
        if PC.SEPARATESKY == 'Y':
            Skylistdic = dict()
            print("Below in addition to regexp, if needed you can enter the starting and ending filenumbers separated by space also.")
            for filt in FiltList:
                filenumbregexp=re.compile(r'.*')
                if filt not in SkyREdic.keys() : 
                    SkyREdic[filt]=ObjRE+'_[Ss]ky.*'  #Setting default to object_sky
                #Ask user again to confirm or change if he/she needs to
                InpfiltRE=raw_input("Enter Regular Expression for the Sky of filters %s (default: %s) : "%(str(filt),SkyREdic[filt])).strip(' ')
                if InpfiltRE :
                    SkyREdic[filt]=InpfiltRE
                    if len(InpfiltRE.split())==3 and is_number(InpfiltRE.split()[1]) and is_number(InpfiltRE.split()[2]):
                        filenumbregexp = re.compile('|'.join([str(i) for i in range(int(InpfiltRE.split()[1]), int(InpfiltRE.split()[2])+1)]))
                        SkyREdic[filt] = InpfiltRE.split()[0]
                        
                regexpSky = re.compile(r''+SkyREdic[filt])
                SkyList = []
                for imgline in imglogFILElines:
                    ImgFrame = Instrument.IdentifyFrame(imgline)
                    if ImgFrame in ['BIAS','LAMP_SPEC']:  # Skip obvious Bias, Wavelength Calibration Lamps etc
                        continue    #Skip these and go to the next object.
                    if regexpSky.search(' '.join(shlex.split(imgline)[0:2])) is not None :
                        if filt == shlex.split(imgline)[FiltColumn]: # Matching filter
                            if filenumbregexp.search(imgline.split()[-1]): # Last column is filenumber
                                SkyList.append(imgline.split()[0])
                Skylistdic[filt] = SkyList  #Saving Sky list for this filter set

        
        #Now if We are doing Spectroscopy, Find the corresponding Lamp lamps also
        if PC.TODO == 'S':
            Lamplistdic = dict()
            print("Below in addition, if needed you can enter the starting and ending filenumbers separated by space.")
            for filt in FiltList:
                filenumbregexp = re.compile(r'.*')
                if filt not in LampREdic.keys() : 
                    LampREdic[filt] = '.*Fe-.*'  #Setting default to *Lamp*
                #Ask user again to confirm or change if he/she needs to
                InpfiltRE = raw_input("Enter Regular Expression for the Lamp of filters %s (default: %s) : "%(str(filt),LampREdic[filt])).strip(' ')
                if InpfiltRE :
                    LampREdic[filt] = InpfiltRE
                    if len(InpfiltRE.split())==3 and is_number(InpfiltRE.split()[1]) and is_number(InpfiltRE.split()[2]):
                        filenumbregexp = re.compile('|'.join([str(i) for i in range(int(InpfiltRE.split()[1]), int(InpfiltRE.split()[2])+1)]))
                        LampREdic[filt] = InpfiltRE.split()[0]

                regexpLamp = re.compile(r''+LampREdic[filt])
                LampList = []
                for imgline in imglogFILElines:
                    ImgFrame = Instrument.IdentifyFrame(imgline)
                    if ImgFrame != 'LAMP_SPEC':  # Skip if it is not a Wavelength Calibration Lamp.
                        continue    #Skip these and go to the next object.
                    if regexpLamp.search(' '.join(shlex.split(imgline)[0:2])) is not None :
                        if filt == shlex.split(imgline)[FiltColumn]: # Matching filter
                            if filenumbregexp.search(imgline.split()[-1]): # Last column is filenumber
                                LampList.append(imgline.split()[0])
                Lamplistdic[filt]=LampList  #Saving Lamp list for this filter set
            
        #Now, load the Object list and write to a file the Obj and corresponding flats/Lamps
        ObjFlatFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-Flat.List'),'w')
        ObjBiasFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-Bias.List'),'w')
        if PC.TODO=='S': ObjLampFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-Lamp.List'),'w')
        if PC.SEPARATESKY=='Y': ObjSkyFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-Sky.List'),'w')
        for Objline in ObjList:
            Name=shlex.split(Objline)[0]
            FiltOrGrism = shlex.split(Objline)[FiltColumn]
            ObjFlatFILE.write(Name+'  '+' '.join(Flatlistdic[FiltOrGrism])+'\n')
            ImgSizeX, ImgSizeY = tuple(shlex.split(Objline)[-3:-1])
            ObjBiasFILE.write(Name+'  '+' '.join(BiasSizeDic[(ImgSizeX, ImgSizeY)])+'\n')
            if PC.TODO=='S' :ObjLampFILE.write(Name+'  '+' '.join(Lamplistdic[FiltOrGrism])+'\n')
            if PC.SEPARATESKY=='Y': ObjSkyFILE.write(Name+'  '+' '.join(Skylistdic[FiltOrGrism])+'\n')
        ObjFlatFILE.close()
        ObjBiasFILE.close()
        # Create the files for all flat files
        if PC.USEALLFLATS == 'Y':
            FiltrFlatFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllFilter-Flat.List'),'w')
            FlatBiasFILE = open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllFlat-Bias.List'),'w')
            for filt in FiltList|CompleteFiltList:
                FiltrFlatFILE.write('"{0}" {1}\n'.format(filt,' '.join(Flatlistdic[filt])))
                for flatimg in Flatlistdic[filt]:
                    ImgSizeX, ImgSizeY = Flatsizedic[flatimg]
                    FlatBiasFILE.write('{0} {1}\n'.format(flatimg,' '.join(BiasSizeDic[(ImgSizeX, ImgSizeY)])))
            FiltrFlatFILE.close()
            FlatBiasFILE.close()
            
        print('Edit, save (if required) and close the Flat,Bias/Lamp/Sky list associations for this night :'+night)
        raw_input("Press Enter to continue...")
        subprocess.call(PC.TEXTEDITOR.split()+[os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-Flat.List')])
        subprocess.call(PC.TEXTEDITOR.split()+[os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-Bias.List')])
        if PC.TODO == 'S': 
            ObjLampFILE.close()
            subprocess.call(PC.TEXTEDITOR.split()+[os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-Lamp.List')])
        if PC.SEPARATESKY == 'Y': 
            ObjSkyFILE.close()
            subprocess.call(PC.TEXTEDITOR.split()+[os.path.join(PC.MOTHERDIR,PC.OUTDIR,night,'AllObjects-Sky.List')])

    print('All nights over...')    
    
def CreateLogFilesFromFits_subrout(PC,hdu=0):
    """ Creates a log file of all the fits file in each directory 
        hdu = 0 is the index of hdu in hdulist to read in fits file"""

    directories = LoadDirectories(PC,CONF=True)
    LogFilename = PC.NIGHTLOGFILE
    # List the set of coulmn entry after the filename in the log. Last three extra columns will be Xsize Ysize and filenumber.
    LogColumns = [PC.OBJECTHDR, PC.EXPTIMEHDR, PC.FILTERHDR, PC.GRISMHDR, PC.LAMPHDR, PC.SLITHDR, PC.DATEHDR, PC.UTHDR, PC.RAHDR, PC.DECHDR, PC.COMMENTHDR] 
    #Final column heads will be: Filename LogColumns Xsize Ysize FileNumber

    RowString = ' "{'+'}" "{'.join(LogColumns)+'}"'  # String with  header keywords in curly brackets " "
    # Load Instrument
    Instrument = InstrumentObject(PC)

    for night in directories:
        print("Working on night : "+night)
        if os.path.isfile(os.path.join(PC.MOTHERDIR,night,LogFilename)):
            print('WARNING: Log file {0} already exist.'.format(os.path.join(PC.MOTHERDIR,night,LogFilename)))
            print('Skipping this directory now. If you want to create new log file, delete existing one or change filename NIGHTLOGFILE= in config file.')
            continue
        os.chdir(os.path.join(PC.MOTHERDIR,night))
        listOFimgs = [f for f in os.listdir('.') if re.match(r'^.*\.fits$', f)]
        listOFimgs = sorted(listOFimgs)#, key=lambda k: int(k[:-5].split('-')[-1]))
        with open(LogFilename,'w') as outfile:
            for i,img in enumerate(listOFimgs):
                hdulist = fits.open(img)
                ImgSizeX, ImgSizeY = hdulist[hdu].data.shape
                prihdr = hdulist[hdu].header
                prihdr = Instrument.StandardiseHeader(prihdr)
                for hkeys in LogColumns :   #Capture if these keywords are missing in header, and replace with -NA-
                    if hkeys not in prihdr : prihdr[hkeys]='-NA-'
                
                outfile.write(img+' '+RowString.format(**prihdr)+' {0} {1} {2}\n'.format(ImgSizeX,ImgSizeY,i))
                hdulist.close()
    print("{0} saved in each night's directory. Edit in manually for errors like ACTIVE filter.".format(LogFilename))

def LoadDirectories(PC,CONF=False):
    """ Loads the directories and return the list of directories to do analysis """
    try :
        directoriesF=open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'directories'),'r')
    except IOError :
        #Creating a text file containing the directories which has fits files in it
        directories = [dirs for dirs in os.walk(PC.MOTHERDIR).next()[1] if glob.glob(os.path.join(PC.MOTHERDIR,dirs,'*.fits'))]
        directories.sort()
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'directories'),'w') as directoriesF : #Creating directories file
            directoriesF.write('\n'.join(directories)+'\n')
        #Now reopening the file to read and proceed
        directoriesF=open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'directories'),'r')

    #Load directories list from the file
    directories = [dirs.rstrip().strip(' ').rstrip('/') for dirs in directoriesF if dirs.strip()] #Removing spaces or trailing / and ignore Blank empty lines
    directoriesF.close()

    if CONF == True :
        #Ask user again to confirm or change if he/she needs to
        InpList=raw_input('Enter the directories to analyse (default: %s) :'%','.join(directories)).strip(' ')
        if InpList : 
            directories=[dirs.rstrip().strip(' ').rstrip('/') for dirs in InpList.split(',')] #Removing spaces or trailing /
            for dirs in list(directories): # iterating over a copy of the list
                if not os.path.isdir(os.path.join(PC.MOTHERDIR,dirs)):
                    print('Cannot find the data directory: {0} in the current directory {1}'.format(dirs,PC.MOTHERDIR))
                    print('WARNING: Removing the non-existing directory : {0} from the list'.format(dirs))
                    directories.remove(dirs)
                          
            with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'directories'),'w') as directoriesF : #Updating directories file
                directoriesF.write('\n'.join(directories)+'\n')


    for dirs in directories:
        #Create a corresponding night directory in OUTPUT directory also if not already present.
        try:
            os.makedirs(os.path.join(PC.MOTHERDIR,PC.OUTDIR,dirs))
        except OSError:
            if os.path.isdir(os.path.join(PC.MOTHERDIR,PC.OUTDIR,dirs)) :
                print("Output directory "+os.path.join(PC.MOTHERDIR,PC.OUTDIR,dirs)+" already exists.\n Everything inside it will be overwritten.")
            else:
                raise
        else:
            print('Created directory :'+os.path.join(PC.MOTHERDIR,PC.OUTDIR,dirs))
        
    if len(directories) == 0 : 
        print('ERROR: No valid directories to reduce data found.')
        print('Atleast one directory containing data should be given as input.')
        sys.exit(1)
    else :
        return directories
    
# To be removed. Nobody uses this to backup!!
def Backup_subrout(PC):
    """ Copies all the files in present directory to the ../PC.BACKUPDIR """
    os.makedirs('../'+PC.BACKUPDIR)
    print("Copying files to ../"+PC.BACKUPDIR)
    os.system('cp -r * ../'+PC.BACKUPDIR)


def KeyboardInterrupt_handler():
    print('\nYou pressed Ctrl+C!')
    print('Stoping the reduction abruptly...')
    sys.exit(2)

def InitialTest(PC):
    """ Few initial tests to see all settings are correct . Input PC is the PipelineConfiguration object"""
    #Check login.cl file exist in ~/iraf/ directory.
    if not os.path.isfile(os.path.expanduser('~/iraf/login.cl')):
        print('login.cl file not found in {0} directory. Copy your login.cl file into {0}'.format(os.path.expanduser('~/iraf/')))
        raise IOError(os.path.expanduser('~/iraf/login.cl'))

    #Check the "set     stdimage        = imt1024" is present in login.cl file.
    with open(os.path.expanduser('~/iraf/login.cl')) as clfile:
        clstdimage = None
        for line in clfile:
            line=line.rstrip().split()
            if len(line) >= 3:
                if line[0]=='set' and line[1]=='stdimage' and line[2] =='=' :
                    clstdimage=line[3]
        if clstdimage != 'imt2048':
            print('WARNING: Image size 2048 not set in login.cl file')
            print('Please set the following line in {0} file.'.format(os.path.expanduser('~/iraf/login.cl')))
            print('set     stdimage        = imt2048')
            print('This is to avoid a certain bug in new version of IRAF and ds9')
            print('This will be an issue in ds9 >7.1 while you open images from pyraf')
            
    #Check for SExtractor installation
    if PC.TODO == 'P':
        with open('/dev/null','w') as devnull:
            try :
                subprocess.call(['sex','--version'],stdout=devnull)
            except OSError:
                try :
                    subprocess.call(['sextractor','--version'],stdout=devnull)
                except OSError:            
                    print('Cannot find the command: sex OR sextractor')
                    print('SExtractor needs to be installed from http://www.astromatic.net/software/sextractor')
                    raise

    #Check LampRepo directory exists.
    if PC.TODO == 'S':
        if not os.path.exists(PC.LAMPREPODIR):
            print('Lamp Repository directory not found: {0}'.format(PC.LAMPREPODIR))
            print('You can obtain {0} LampRepo directory by extracting data.tar.gz'.format(PC.INSTRUMENT))
            print('Please add the correct path to Lamp Repository before proceeding.')
            raise IOError(PC.LAMPREPODIR)



class PipelineConfig(object):
    """ This class is just a collection of variables required to run the pipeline """
    def __init__(self,ConfigFilename = None):
        if ConfigFilename is not None:
            self.LoadFromFile(ConfigFilename)

        self.MOTHERDIR = os.getcwd()
        # variables for helper functions
        self._FrameDatabaseInitialised = False
        self.currentnight = None

    def LoadFromFile(self,ConfigFilename):
        """ Loads the configuration from the input file """
        with open(ConfigFilename,'r') as configfile:
            for con in configfile:
                con = con.rstrip()
                if len(con.split()) >= 2 :
                    
                    if con.split()[0] == "INSTRUMENT=" :
                        self.INSTRUMENT = con.split()[1]
                    elif con.split()[0] == "TODO=" :
                        self.TODO = con.split()[1]
                    elif con.split()[0] == "OUTPUTDIR=" :
                        self.OUTDIR = con.split()[1]
                    elif con.split()[0] == "VERBOSE=" :
                        self.VER = con.split()[1]
                    elif con.split()[0] == "TEXTEDITOR=" :
                        self.TEXTEDITOR = shlex.split(con)[1]
                    elif con.split()[0] == "OVERSCAN=" :
                        self.OVERSCAN = con.split()[1]
                    elif con.split()[0] == "COMBINEIMGS=" :
                        self.IMGCOMBINE = con.split()[1][0].upper()
                    elif con.split()[0] == "IMGCOMBMETHOD=" :
                        self.IMGCOMBMETHOD = con.split()[1]
                    elif con.split()[0] == "FLATSTATSECTION=" :
                        self.FLATSTATSECTION = con.split()[1]
                    elif con.split()[0] == "USEALLFLATS=" :
                        self.USEALLFLATS = con.split()[1][0].upper()

                    elif con.split()[0] == "SEPARATE_SKY=" :
                        self.SEPARATESKY = con.split()[1][0].upper()
                    elif con.split()[0] == "GRADIENT_REMOVE=" :
                        self.GRADREMOVE = con.split()[1][0].upper()
                    elif con.split()[0] == "GRAD_FILT_SIZE=" :
                        self.MEDSMOOTHSIZE = (int(con.split()[1]), int(con.split()[2]))

                    elif con.split()[0] == "BPMPHOTO=" :
                        self.PhotBadPixelMaskName = con.split()[1]
                    elif con.split()[0] == "BPMSPEC=" :
                        self.SpecBadPixelMaskName = con.split()[1]
        
                    elif con.split()[0] == "THRESHOLD=" :
                        self.THRESHOLD = float(con.split()[1])
                    elif con.split()[0] == "EPADU=" :
                        self.EPADU = float(con.split()[1])
                    elif con.split()[0] == "READNOISE=" :
                        self.READNOISE = float(con.split()[1])
                    elif con.split()[0] == "DATAMAX=" :
                        self.DATAMAX = con.split()[1]
                    elif con.split()[0] == "XYMATCHMIN=" :
                        self.XYMATCHMIN = int(con.split()[1])
            
                    elif con.split()[0] == "APERTURE=" :
                        self.APERTURE = con.split()[1]
                    elif con.split()[0] == "ANNULUS=" :
                        self.ANNULUS = con.split()[1]
                    elif con.split()[0] == "DANNULUS=" :
                        self.DANNULUS = con.split()[1]

                    elif con.split()[0] == "EXPTIME=" :
                        self.EXPTIMEHDR = con.split()[1]
                    elif con.split()[0] == "FILTER=" :
                        self.FILTERHDR = con.split()[1]
                    elif con.split()[0] == "GRISM=" :
                        self.GRISMHDR = con.split()[1]
                    elif con.split()[0] == "SLIT=" : 
                        self.SLITHDR = con.split()[1]
                    elif con.split()[0] == "LAMP=" : 
                        self.LAMPHDR = con.split()[1]

                    elif con.split()[0] == "UT=" :
                        self.UTHDR = con.split()[1]
                    elif con.split()[0] == "ODATE=" :
                        self.DATEHDR = con.split()[1]

                    elif con.split()[0] == "OBJECT=" : 
                        self.OBJECTHDR = con.split()[1]
                    elif con.split()[0] == "COMMENT=" :
                        self.COMMENTHDR = con.split()[1]
                    elif con.split()[0] == "RA_HDR=" : 
                        self.RAHDR = con.split()[1]
                    elif con.split()[0] == "DEC_HDR=" :
                        self.DECHDR = con.split()[1]


                    elif con.split()[0] == "OUTPUT=" :
                        self.OUTPUTFILE = con.split()[1]
                    elif con.split()[0] == "BACKUP=" :
                        self.BACKUPDIR = con.split()[1]

                    elif con.split()[0] == "CONVOLVEIMG=" :
                        self.CONVOLVEIMG = con.split()[1]
                    elif con.split()[0] == "DOPSF=" :
                        self.DOPSF = con.split()[1]

                    elif con.split()[0] == "LAMPDIRECTORY=" :
                        self.LAMPREPODIR = con.split()[1]
                    elif con.split()[0] == "SPECAPERTURE=" :
                        self.SPECAPERTURE = con.split()[1]
                    elif con.split()[0] == "SPECAPERTURE_LLIMIT=" :
                        self.SPECAPERTURE_LLIMIT = con.split()[1]
                    elif con.split()[0] == "SPECAPERTURE_ULIMIT=" :
                        self.SPECAPERTURE_ULIMIT = con.split()[1]


                    elif con.split()[0] == "BACKGROUND=" :
                        self.BACKGROUND = con.split()[1]
                    elif con.split()[0] == "TRACEFUNC=" :
                        self.TRACEFUNC = con.split()[1]
                    elif con.split()[0] == "TRACEORDER=" :
                        self.TRACEORDER = con.split()[1]
                    elif con.split()[0] == "NORMFUNC=" :  # Not used yet
                        self.NORMFUNC = con.split()[1]
                    elif con.split()[0] == "NORMORDER=" :  # Not used yet
                        self.NORMORDER = con.split()[1]
                    elif con.split()[0] == "SCOMBINE=" :
                        self.SCOMBINE = con.split()[1]
                    elif con.split()[0] == "DISPAXIS=" :
                        self.DISPAXIS = con.split()[1]
                    elif con.split()[0] == "REMOVE_CONTINUUM_GRAD=" :
                        self.CONTINUUMGRADREMOVE = con.split()[1][0].upper()
                    elif con.split()[0] == "CONT_GRAD_FILT_SIZE=" :
                        self.DVDMEDSMOOTHSIZE = (int(con.split()[1]), int(con.split()[2]))

                    elif con.split()[0] == "NIGHTLOGFILE=" :
                        self.NIGHTLOGFILE = con.split()[1]


    # Some helper functions for obtaining the full path to files.
    # They are not necessory for many steps.
    def GetFullPath(self,filename):
        """ Returns full path to filename.
        If the filename is a raw image, then it will give adress to the raw data directory.
        If the filename is new, then it will give adress to the output directory of PC.currentnight""" 
        if not self._FrameDatabaseInitialised:
            #Initialise the database. # Don't call this if the Logfiles are not already created (in 0th step)
            self.FrameDatabase = {}
            for night in LoadDirectories(self,CONF=False):
                with open(os.path.join(self.MOTHERDIR,night,self.NIGHTLOGFILE),'r') as imglogFILE:
                    # Skip blank lines and Commented out lines with #
                    images = [imageLINE.split()[0] for imageLINE in imglogFILE if \
                              ((imageLINE.strip() is not '') and (imageLINE[0] !='#'))]
                self.FrameDatabase[night] = set(images)
            self._FrameDatabaseInitialised = True

        if self.currentnight is None:
            # Raise error
            print('ERROR: Current directory not defined PC.currentnight')
            raise ValueError('ERROR: PC.currentnight not defined')

        # All fine now,  give full path intelegently.
        if filename in self.FrameDatabase[self.currentnight]:
            return os.path.join(self.MOTHERDIR,self.currentnight,filename)
        else:
            return os.path.join(self.MOTHERDIR,self.OUTDIR,self.currentnight,filename)


################################# Instrument Definitions  ###########
# Details as well as the name of each function calls of each instrument
# each step in pipeline is a function call with one input argument PC

# Ascii art courtesy : http://patorjk.com/software/taag/

class InstrumentObject(object):
    """ All the definitions and procedures of each instrument is defined inside this class."""
    def __init__(self,PC):
        self.PC = PC
        self.Name = self.PC.INSTRUMENT
        self.Banner = InstrumentDictionaries[self.Name]['Banner']
        self.About = InstrumentDictionaries[self.Name]['About']
        self.Steps = InstrumentDictionaries[self.Name]['Steps']
        self.StepsToRun = self.GetListofStepsToRun()

    def GetListofStepsToRun(self):
        """ Returns the list of steps to run based onPC, the PipelineConfiguration object """
        StepsToRun = []
        if self.Name in ['HFOSC','IFOSC']:
            StepsToRun += [0, 1, 2, 3, 4, 5]
            if self.PC.IMGCOMBINE == 'Y':
                StepsToRun += [ 6 ]
            if self.PC.TODO == 'P':
                StepsToRun += [7, 8, 9, 10]
            elif self.PC.TODO == 'S':
                StepsToRun += [11, 12]
        else:
            print('Unknown Instrument')

        return StepsToRun

    def StandardiseHeader(self,prihdr):
        """ Return the heared object after standardising the values in it, specific to instrument."""
        if self.Name == 'HFOSC':
            ut_sec = float(prihdr[self.PC.UTHDR])
            prihdr[self.PC.UTHDR] = str(datetime.timedelta(seconds=ut_sec))  # Converting to HH:MM:SS
        elif self.Name == 'IFOSC':
            prihdr[self.PC.EXPTIMEHDR] = float(prihdr[self.PC.EXPTIMEHDR])/1000.0  # Conver ms to sec for IFOSC
        return prihdr

    def IdentifyFrame(self,ObjectLine):
        """ Returns what object the input frame is based on the ObjectLine 
        'Filename, PC.OBJECTHDR, PC.EXPTIMEHDR, PC.FILTERHDR, PC.GRISMHDR, PC.LAMPHDR, PC.SLITHDR, PC.DATEHDR, PC.UTHDR, PC.RAHDR, PC.DECHDR, PC.COMMENTHDR, Xsize, Ysize, FileNumber' """
        Frame = 'UNKNOWN'
        if self.Name == 'HFOSC':
            if float(shlex.split(ObjectLine)[2]) == 0:
                Frame = 'BIAS'
            elif 'GRISM' in shlex.split(ObjectLine)[4].upper():
                if 'FE-' in shlex.split(ObjectLine)[5].upper():
                    Frame = 'LAMP_SPEC'
                elif 'HALOGEN' in shlex.split(ObjectLine)[5].upper():
                    Frame = 'FLAT_SPEC'
                else:
                    Frame = 'OBJECT_SPEC'
            else:
                Frame = 'OBJECT_IMG'
            
        if self.Name == 'IFOSC':
            if float(shlex.split(ObjectLine)[2]) == 0:
                Frame = 'BIAS'
            elif 'GRISM' in shlex.split(ObjectLine)[4].upper():
                if 'FE-' in shlex.split(ObjectLine)[5].upper():
                    Frame = 'LAMP_SPEC'
                elif 'HALOGEN' in shlex.split(ObjectLine)[5].upper():
                    Frame = 'FLAT_SPEC'
                else:
                    Frame = 'OBJECT_SPEC'
            else:
                Frame = 'OBJECT_IMG'

        return Frame
        
        
                             

InstrumentDictionaries = {'HFOSC':{
    'Banner' :"""
.---.  .---.  ________     ,-----.       .-'''-.     _______    
|   |  |_ _| |        |  .'  .-,  '.    / _     \   /   __  \   
|   |  ( ' ) |   .----' / ,-.|  \ _ \  (`' )/`--'  | ,_/  \__)  
|   '-(_{;}_)|  _|____ ;  \  '_ /  | :(_ o _).   ,-./  )        
|      (_,_) |_( )_   ||  _`,/ \ _/  | (_,_). '. \  '_ '`)      
| _ _--.   | (_ o._)__|: (  '\_/ \   ;.---.  \  : > (_)  )  __  
|( ' ) |   | |(_,_)     \ `"/  \  ) / \    `-'  |(  .  .-'_/  ) 
(_{;}_)|   | |   |       '. \_/``".'   \       /  `-'`-'     /  
'(_,_) '---' '---'         '-----'      `-...-'     `._____.'   
                                                                
                                          Data Reduction Pipeline...
""",
    'About' : "HFOSC Optical Imaging and Spectrometer on 2-m HCT, IAO, IIA",
    'Steps' : {
        0 : {'Menu': 'Generate Log files of fits files in each directory.',
             'RunMessage': "RUNNING TASK:0  Generating log files of fits files in each directory..",
             'function': CreateLogFilesFromFits_subrout },
        1 : {'Menu': 'Selection of object frames, Flats/Sky/Lamps etc. to reduce',
             'RunMessage': "RUNNING TASK:1 Selecting object frames, Flats/Sky/Lamps etc..",
             'function': SelectionofFrames_subrout },
        2 : {'Menu': 'Visually inspect and/or reject object images one by one.',
             'RunMessage': "RUNNING TASK:2  Visual inspection and/or rejection of object frames..",
             'function': Manual_InspectObj_subrout},
        3 : {'Menu': 'Visually inspect and/or reject Bias/Flats/Sky/Lamps one by one.',
             'RunMessage': "RUNNING TASK:3  Visual inspection and/or rejection of Bias/Flats/Sky/Lamps frames..",
             'function': Manual_InspectCalframes_subrout},
        4 : {'Menu': 'Subtract Overscan,Bias/sky',
             'RunMessage': "RUNNING TASK:4  Subtracting biases or sky..",
             'function': Bias_Subtraction_subrout},
        5 : {'Menu': 'Apply Flat Correction and/or Bad pixel interpolation',
             'RunMessage': "RUNNING TASK:5  Flat correction etc..",
             'function': Flat_basicCorrections_subrout},
        6 : {'Menu': 'Align and combine images.',
             'RunMessage': "RUNNING TASK:6  Aligning and combining images..",
             'function': AlignNcombine_subrout },
        7 : {'Menu': 'Make the list of images, Images4Photo.in to do Photometry.',
             'RunMessage': "RUNNING TASK:7  Makeing list of images, Images4Photo.in to do Photometry..",
             'function': Createlist_subrout },
        8 : {'Menu': 'Select Stars and Sky region of the field on first image',
             'RunMessage': "RUNNING TASK:8  Selecting Stars and Sky region of the field from first image..",
             'function': Star_sky_subrout },
        9 : {'Menu': 'Create Sextracter config file & coordinate output of first image.',
             'RunMessage': "RUNNING TASK:9  Create Sextracter config file & coordinate output of first image..",
             'function': Sextractor_subrout },
        10 : {'Menu': 'Do Photometry',
             'RunMessage': "RUNNING TASK:10 Doing Photometry..",
             'function': Photometry },
        11: {'Menu': 'Input Spectrum pair subtraction and filenames.',
             'RunMessage': "RUNNING TASK:11  Inputing Spectrum pair subtraction and filenames...",
             'function': SpectralPairSubtraction_subrout },
        12: {'Menu': 'Extract wavelength calibrated 1D spectra from image.',
             'RunMessage': "RUNNING TASK:12  Extracting wavelength calibrated 1D spectra..",
             'function': SpectralExtraction_subrout }
    }
},
'IFOSC':{
    'Banner' :"""
.-./`)  ________     ,-----.       .-'''-.     _______    
\ .-.')|        |  .'  .-,  '.    / _     \   /   __  \   
/ `-' \|   .----' / ,-.|  \ _ \  (`' )/`--'  | ,_/  \__)  
 `-'`"`|  _|____ ;  \  '_ /  | :(_ o _).   ,-./  )        
 .---. |_( )_   ||  _`,/ \ _/  | (_,_). '. \  '_ '`)      
 |   | (_ o._)__|: (  '\_/ \   ;.---.  \  : > (_)  )  __  
 |   | |(_,_)     \ `"/  \  ) / \    `-'  |(  .  .-'_/  ) 
 |   | |   |       '. \_/``".'   \       /  `-'`-'     /  
 '---' '---'         '-----'      `-...-'     `._____.'   
                                                                
                                          Data Reduction Pipeline...
""",
    'About' : "IFOSC Optical Imaging and Spectrometer on 2-m IGO, IUCAA",
    'Steps' : {
        0 : {'Menu': 'Generate Log files of fits files in each directory.',
             'RunMessage': "RUNNING TASK:0  Generating log files of fits files in each directory..",
             'function': CreateLogFilesFromFits_subrout },
        1 : {'Menu': 'Selection of object frames, Flats/Sky/Lamps etc. to reduce',
             'RunMessage': "RUNNING TASK:1 Selecting object frames, Flats/Sky/Lamps etc..",
             'function': SelectionofFrames_subrout },
        2 : {'Menu': 'Visually inspect and/or reject object images one by one.',
             'RunMessage': "RUNNING TASK:2  Visual inspection and/or rejection of object frames..",
             'function': Manual_InspectObj_subrout},
        3 : {'Menu': 'Visually inspect and/or reject Bias/Flats/Sky/Lamps one by one.',
             'RunMessage': "RUNNING TASK:3  Visual inspection and/or rejection of Bias/Flats/Sky/Lamps frames..",
             'function': Manual_InspectCalframes_subrout},
        4 : {'Menu': 'Subtract Overscan,Bias/sky',
             'RunMessage': "RUNNING TASK:4  Subtracting biases and/or I fringe subtraction etc..",
             'function': Bias_Subtraction_subrout},
        5 : {'Menu': 'Apply Flat Correction and/or Bad pixel interpolation',
             'RunMessage': "RUNNING TASK:5  Flat correction etc..",
             'function': Flat_basicCorrections_subrout},
        6 : {'Menu': 'Align and combine images.',
             'RunMessage': "RUNNING TASK:6  Aligning and combining images..",
             'function': AlignNcombine_subrout },
        7 : {'Menu': 'Make the list of images, Images4Photo.in to do Photometry.',
             'RunMessage': "RUNNING TASK:7  Makeing list of images, Images4Photo.in to do Photometry..",
             'function': Createlist_subrout },
        8 : {'Menu': 'Select Stars and Sky region of the field on first image',
             'RunMessage': "RUNNING TASK:8  Selecting Stars and Sky region of the field from first image..",
             'function': Star_sky_subrout },
        9 : {'Menu': 'Create Sextracter config file & coordinate output of first image.',
             'RunMessage': "RUNNING TASK:9  Create Sextracter config file & coordinate output of first image..",
             'function': Sextractor_subrout },
        10 : {'Menu': 'Do Photometry',
             'RunMessage': "RUNNING TASK:10 Doing Photometry..",
             'function': Photometry },
        11: {'Menu': 'Input Spectrum pair subtraction and filenames.',
             'RunMessage': "RUNNING TASK:11  Inputing Spectrum pair subtraction and filenames...",
             'function': SpectralPairSubtraction_subrout },
        12: {'Menu': 'Extract wavelength calibrated 1D spectra from image.',
             'RunMessage': "RUNNING TASK:12  Extracting wavelength calibrated 1D spectra..",
             'function': SpectralExtraction_subrout }
    }
}
} # END of Instrument Definitions..#######

                        
        
#-----Main Program Calls Begins here........................
# If you are using the functions seperately from module, see the Class PipelineConfig 
# to create the PipelineConfiguration instance which has to be passed on to many function calls.

def main():
    """ The Pipeline main funtion to run while run from terminal as a standalone pipeline """

    if len(sys.argv) < 2 :
        print('-'*10)
        print('Usage : {0} ScriptSettings.conf'.format(sys.argv[0]))
        print('where,')
        print('     ScriptSettings.conf is the configuration file for this run of reduction pipeline')
        print(' ')
        print("Note: This script should be run from the directory containing all the night's data directories.")
        print('-'*10)
        sys.exit(1)

    try : 
        PC = PipelineConfig(ConfigFilename = sys.argv[1])
    except IOError :
        print ("Error: Cannot open the file "+sys.argv[1]+". Setup the config file based on ScriptSettings.conf file correctly, before running the script.")
        sys.exit(1)

    # Load Instrument
    Instrument = InstrumentObject(PC)
    # Print Welcome Banner.. 
    print("Time: {0}".format(time.strftime("%c")))
    print(Instrument.Banner)
    print(Instrument.About)

    print('Running Initial Tests...')
    try :
        InitialTest(PC)  # Some important Initial tests for Sanity before running the code.
    except Exception as e :
        print(e)
        sys.exit(1)


    print('-'*10)
    print('IMP: All outputs of this run will be written to the directory '+os.path.join(PC.MOTHERDIR,PC.OUTDIR)+'\n')
    try:
        os.makedirs(os.path.join(PC.MOTHERDIR,PC.OUTDIR))
    except OSError:
        if os.path.isdir(os.path.join(PC.MOTHERDIR,PC.OUTDIR)) :
            print("WARNING : Output directory "+os.path.join(PC.MOTHERDIR,PC.OUTDIR)+" already exists.\n Everything inside it will be overwritten. Be warned...")
        else:
            raise
    else:
        print('Created directory :'+os.path.join(PC.MOTHERDIR,PC.OUTDIR))
        
        
    if PC.TODO == 'P' : todoinwords = 'Photometry'
    elif PC.TODO == 'S' : todoinwords = 'Spectroscopy'


    print("\n Very Very Important: Backup your files first. Don't proceed without backup.\n")
    print(" ---------------- Welcome to \033[91m {0} {1} \033[0m Pipeline --------------- \n".format(PC.INSTRUMENT,todoinwords))
    print("Enter the Serial numbers (space separated if more than one task in succession) \n")
    for i in Instrument.StepsToRun:
        print("{0}  {1} \n".format(i,Instrument.Steps[i]['Menu']))
    print("--------------------------------------------------------------- \n")

    try:
        with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'StepsFinished'),'r') as stepsoverFLS:
            StepsOver = stepsoverFLS.read()
    except IOError:
        StepsOver = 'Nothing...'
    print('Steps you have already finished : '+StepsOver)    

    try:
        todo = raw_input('Enter the list : ')
        todo = todo.split()

        for task in todo :
            print("Time now: {0}".format(time.strftime("%c")))
            CalledTheTask=True
            try:
                task = int(task)
                print(Instrument.Steps[task]['RunMessage'])
            except ValueError, IndexError :
                print('Cannot understand the input task: {0}'.format(task))
                print('Skipping task '+task)
                CalledTheTask=False
            else: # Run the task..
                Instrument.Steps[task]['function'](PC)

            if CalledTheTask :
                with open(os.path.join(PC.MOTHERDIR,PC.OUTDIR,'StepsFinished'),'a') as stepsoverFLS:
                    stepsoverFLS.write(str(task)+' ')
    except KeyboardInterrupt :
        KeyboardInterrupt_handler()

    print("All tasks over....Enjoy!!!")
            
if __name__ == "__main__":
    main()
