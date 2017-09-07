import os, re, time, datetime, sys, math, fnmatch
from os.path import join, getsize
from datetime import date
from collections import deque
import Utilities
#from multiprocessing import Process	# No longer needed now SERPent is parallel
#from multiprocessing import Pool
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import math, time, datetime
from numpy import *
import itertools
from time import gmtime, strftime, localtime
import pyfits
import numpy as np
### inputs ###
AIPS.userno = 1002
i=1
auto_rms = True
rms = 4.73189515179e-05
edge = 10
rms_box=250
postfix = 'JVLA'
useSAD = False ## If True SAD will be used otherwise blobcat is used (C. Hales+12)
ds9 = True ## Writes out ds9 region file
write_blobs = True ## Writes new blob images

if write_blobs == True:
    write_blobs = '--write'
else:
    write_blobs = ''
if ds9 == True:
    ds9 = '--ds9'
else:
    ds9 = ''

def SAD_fit_remove(files,postfix):
    if os.path.isfile('catalogue_%s.csv' % postfix) == False:
        s = 'Catalog_name rms_{0} BMAJ_{0} BMIN_{0} BPA_{0} #_{0}      Peak_{0}    Dpeak_{0}     Flux_{0}    Dflux_{0}    RA---SIN_{0}   DEC--SIN_{0}  Dx_{0}      Dy_{0}       Maj_{0}     Min_{0}      PA_{0}    Dmaj_{0}    Dmin_{0}    Dpa_{0} #_{0}  MAJ-fit_{0} MIN-fit_{0} PA-fit_{0}    MAJ-dec_{0} MIN-dec_{0}  PA-dec_{0}  R_{0} MAJ-low_{0} MIN-low_{0}  PA-low_{0}    MAJ-hi_{0}  MIN-hi_{0}  PA-hi_{0}    Xpix_{0}   Ypix_{0}   MAXresid_{0}\n'.format(postfix)
        s = ' '.join(s.split())+'\n'
        s = s.replace(' ',',')
        os.system('touch catalogue_%s.csv' % postfix)
        text_file = open('catalogue_%s.csv' % postfix,'a')
        text_file.write(s)
    for j in files:
        with open(j) as f:
            x = f.read().splitlines()
        remove = ['#','Component']
        for i in remove:
            x = [y for y in x if not i in y]
        x = [y for y in x if y != ' ']
        x = [y.lstrip() for y in x]
        x = ' '.join(x).replace('(',' ').replace(')',' ')
        x = ' '.join(x.split())
        text_file.write(x.replace(' ',',')+'\n')

def blobcat_fit_remove(files,postfix):
    if os.path.isfile('catalogue_%s.csv' % postfix) == False:
        s = 'Catalog_name rms_{0} BMAJ_{0} BMIN_{0} BPA_{0} ID_{0}  npix_{0} x_p_{0} y_p_{0}           RA_p_{0}  Dec_p_{0}     RA_err_{0}    Dec_err_{0}   x_c_{0}   y_c_{0} RA_c_{0}           Dec_c_{0}      cFlag_{0}   x_wc_{0}   y_wc_{0}   RA_wc_{0} Dec_wc_{0}   wcFlag_{0}       xmin_{0}       xmax_{0}       ymin_{0}       ymax_{0}        rms_{0}    BWScorr_{0} M_{0}    SNR_OBS_{0}     SNR_FIT_{0}         SNR_{0}        S_p_OBS_{0}        S_p_FIT_{0}            S_p_{0}         S_p_CB_{0}      S_p_CBBWS_{0}  S_p_CBBWS_err_{0}      S_int_OBS_{0}    S_int_OBSCB_{0}          S_int_{0}       S_int_CB_{0}   S_int_CB_err_{0}    R_EST_{0}     VisArea_{0}\n'.format(postfix)
        s = ' '.join(s.split())+'\n'
        s = s.replace(' ',',')
        os.system('touch catalogue_%s.csv' % postfix)
        text_file = open('catalogue_%s.csv' % postfix,'a')
        text_file.write(s)
    for j in files:
        with open(j) as f:
            i = 0
            prefix = []
            for x in f:
                i = i +1
                if i>6:
                    x = [',,,,']+[x]
                    x = ' '.join(x).split()
                    x = ' '.join(x)
                    text_file.write(x.replace(' ',',')+'\n')
                elif i == 6:
                    x = prefix +[x]
                    x = ' '.join(x).split()
                    x = ' '.join(x)
                    text_file.write(x.replace(' ',',')+'\n')
                else:
                    prefix = prefix + [x]

os.system('rm catalogue_%s.csv detections.txt' % postfix)

detections = []

if useSAD == True:
    for file in os.listdir('./'):
        if file.endswith('_casa.fits'):
            fitld = AIPSTask('FITLD')
            hduheader = pyfits.open(file)[0].header
            print file
            try:
                data = np.array(pyfits.open(file)[0].data[0,0,edge:edge+rms_box,edge:edge+rms_box])
            except IndexError:
                data = np.array(pyfits.open(file)[0].data[edge:edge+rms_box,edge:edge+rms_box])
            if auto_rms == True:
                rms = float(np.sqrt(np.mean(data**2)))
                print rms
            fitld.datain = 'PWD:%s' % file
            fitld.outname = str(i)
            fitld.outclass = 'IM'
            fitld.go()
            image = AIPSImage(str(i),'IM',1,1)
            sad = AIPSTask('SAD')
            sad.cparm[1:] = 8*rms, 6*rms, 5*rms
            sad.indata = image
            sad.in2data = image
            sad.blc[1:] = edge,edge
            sad.trc[1:] = int(hduheader['NAXIS1'])-edge,int(hduheader['NAXIS2'])-edge
            sad.dparm[1] = 6*rms
            sad.fitout = 'PWD:%s.fitout' % file
            sad.go()
            image.zap()
            lines = open('%s.fitout' % file).readlines()
            try:
                BMAJ = hduheader['BMAJ']/hduheader['CDELT2'] ## assuming cell is same size on both axes
                BMIN = hduheader['BMIN']/hduheader['CDELT2']
                BPA = hduheader['BPA']
            except KeyError:
                print 'Run casa_convert.py first to get beam parameters into header'
                sys.exit()
            if len(lines) > 24:
                detections = detections + [file]
                open('%s_r.fitout' % file, 'w').writelines(file[:8]+'\n')
                open('%s_r.fitout' % file, 'a').writelines(str(rms)+'\n')
                open('%s_r.fitout' % file, 'a').writelines(str(BMAJ)+'\n')
                open('%s_r.fitout' % file, 'a').writelines(str(BMIN)+'\n')
                open('%s_r.fitout' % file, 'a').writelines(str(BPA)+'\n')
                open('%s_r.fitout' % file, 'a').writelines(lines[18:])
            os.system('rm %s.fitout' % file)
    catalog_list = []
    for file in os.listdir('./'):
        if file.endswith('fitout'):
            catalog_list = catalog_list + [file]
    SAD_fit_remove(catalog_list,postfix)
    os.system('rm *fitout')
else:
    for file in os.listdir('./'):
        if file.endswith('_casa.fits'):
            fitld = AIPSTask('FITLD')
            hduheader = pyfits.open(file)[0].header
            print file
            try:
                data = np.array(pyfits.open(file)[0].data[0,0,edge:edge+rms_box,edge:edge+rms_box])
            except IndexError:
                data = np.array(pyfits.open(file)[0].data[edge:edge+rms_box,edge:edge+rms_box])
            if auto_rms == True:
                rms = float(np.sqrt(np.mean(data**2)))
                print rms
            os.system('python blobcat.py --dSNR=6 --rmsval=%f --edgemin=%d %s %s %s' % (rms,int(edge),ds9,write_blobs,file))
            print 'python blobcat.py --dSNR=6 --rmsval=%f --edgemin=%d %s %s %s' % (rms,int(edge),ds9,write_blobs,file)
            lines = open('%s_blobs.txt' % file[:-5]).readlines()
            try:
                BMAJ = hduheader['BMAJ']/hduheader['CDELT2'] ## assuming cell is same size on both axes
                BMIN = hduheader['BMIN']/hduheader['CDELT2']
                BPA = hduheader['BPA']
            except KeyError:
                print 'Run casa_convert.py first to get beam parameters into header'
                sys.exit()
            if len(lines) > 30:
                detections = detections + [file]
                open('%s_r.blobs' % file, 'w').writelines(file[:8]+'\n')
                open('%s_r.blobs' % file, 'a').writelines(str(rms)+'\n')
                open('%s_r.blobs' % file, 'a').writelines(str(BMAJ)+'\n')
                open('%s_r.blobs' % file, 'a').writelines(str(BMIN)+'\n')
                open('%s_r.blobs' % file, 'a').writelines(str(BPA)+'\n')
                open('%s_r.blobs' % file, 'a').writelines(lines[30:])
            os.system('rm %s_blobs.txt' % file[:-5])
    catalog_list = []
    for file in os.listdir('./'):
        if file.endswith('.blobs'):
            catalog_list = catalog_list + [file]
    blobcat_fit_remove(catalog_list,postfix)
    os.system('rm *.blobs')

os.system('touch detections.txt')
thefile = open('detections.txt', 'w')
for item in detections:
  thefile.write("%s\n" % item)
