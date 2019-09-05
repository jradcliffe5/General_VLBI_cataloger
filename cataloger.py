import os, re, time, datetime, sys, math, fnmatch
from os.path import join, getsize
from datetime import date
from collections import deque
import math, time, datetime
from numpy import *
import itertools
from time import gmtime, strftime, localtime
try:
    import astropy.io.fits as pyfits
except ImportError:
    try:
        print('Using pyfits')
        import pyfits
    except ImportError:
        print('No pyfits or astropy installed... exiting!')
        exit()
import numpy as np
from inspect import getsourcefile
from os.path import abspath

file_path = abspath(getsourcefile(lambda:0)).split('cataloger.py')[0]

def headless(inputfile):
	''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
	INPUTFILE = open(inputfile, "r")
	control = {}
	# a few useful regular expressions
	newline = re.compile(r'\n')
	space = re.compile(r'\s')
	char = re.compile(r'\w')
	comment = re.compile(r'#.*')
	# parse the input file assuming '=' is used to separate names from values
	for line in INPUTFILE:
		if char.match(line):
			line = comment.sub(r'', line)
			line = line.replace("'", '')
			(param, value) = line.split('=')
			param = newline.sub(r'', param)
			param = param.strip()
			param = space.sub(r'', param)
			value = newline.sub(r'', value)
			value = value.replace(' ','').strip()
			valuelist = value.split(',')
			if len(valuelist) == 1:
				if valuelist[0] == '0' or valuelist[0]=='1' or valuelist[0]=='2':
					control[param] = int(valuelist[0])
				else:
					control[param] = str(valuelist[0])
			else:
				control[param] = ','.join(valuelist)
	return control

inputs = headless('catalog_inputs.txt')

auto_rms = inputs['auto_rms']
rms = float(inputs['rms'])
edge = int(inputs['edge'])
rms_box = int(inputs['rms_box'])
S_N_ratio = float(inputs['S_N_ratio'])
postfix = str(inputs['postfix'])
shorthand = inputs['shorthand']
useSAD = inputs['useSAD']
ds9 = inputs['ds9']
write_blobs = inputs['write_blobs']
run_BANE = inputs['run_BANE']
use_BANE_rms = inputs['use_BANE_rms']
SNR_flood = float(inputs['SNR_flood'])
pmep = float(inputs['pmep'])                  ## Max estimate pixellation error
ppe = float(inputs['ppe'])
cpeRA = float(inputs['cpeRA'])                ## Error in phase cal RA (arcsec)
cpeDec = float(inputs['cpeDec'])             ## Error in phase cal Dec (arcsec)
pasbe = float(inputs['pasbe'])                ## Surface Brightness error from calibration in per cent
split_catalogues = str(inputs['split_catalogues'])

### Test for running with Parseltongue
try:
	from AIPS import AIPS, AIPSDisk
	from AIPSTask import AIPSTask, AIPSList
	from AIPSData import AIPSUVData, AIPSImage, AIPSCat
	from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
	AIPS.userno = int(inputs['AIPS_user'])
except ImportError:
	print('No AIPS/Parseltongue available, using BLOBCAT instead')
	useSAD = 'False'

'''
# Inputs
## General Inputs
auto_rms = True            ## Determines the rms automatically
rms_box= 250               ## Size of rms box if auto_rms is True (from TLC currently)
rms = 4.4e-05              ## Constant rms value if auto_rms is False
S_N_ratio = 6.             ## Peak S/N ratio
edge = 10                  ## Edge number of pixels not to consider for cataloging
shorthand = False          ## If true, catalog names will be appended to the first 8 characters
useSAD = False             ## If True SAD will be used, otherwise blobcat is used (C. Hales+12)
postfix = 'natural_weight' ## Append this to the catalog name and column rows

## Inputs for SAD
AIPS_user = 1002

## Inputs for blobcat
### Cataloguing options
SNR_flood = 3.              ## SNR ratio to flood down to 2.5-3 x rms is usually ok
pmep = 1.                  ## Max estimate pixellation error
cpeRA = 0.005              ## Error in phase cal RA (arcsec)
cpeDec = 0.005             ## Error in phase cal Dec (arcsec)
pasbe = 0.2                ## Surface Brightness error from calibration in per cent
run_BANE = True            ## Runs BANE (aegean) to create a rms map before cataloging
use_BANE_rms = True        ## Takes rms of each file (needs to be appended with _rms.fits)

### Output options
ds9 = True                 ## Writes out ds9 region file
write_blobs = True         ## Writes new blob images

'''
i=1

if write_blobs == 'True':
	write_blobs = '--write'
else:
	write_blobs = ''
if ds9 == 'True':
	ds9 = '--ds9'
else:
	ds9 = ''

def SAD_fit_remove(files,postfix):
	if os.path.isfile('catalogue_SAD_%s.csv' % postfix) == False:
		s = 'Catalog_name_{0} rms_{0} BMAJ_{0} BMIN_{0} BPA_{0} #_{0}      Peak_{0}    Dpeak_{0}     Flux_{0}    Dflux_{0}    RA---SIN_{0}   DEC--SIN_{0}  Dx_{0}      Dy_{0}       Maj_{0}     Min_{0}      PA_{0}    Dmaj_{0}    Dmin_{0}    Dpa_{0} #_{0}  MAJ-fit_{0} MIN-fit_{0} PA-fit_{0}    MAJ-dec_{0} MIN-dec_{0}  PA-dec_{0}  R_{0} MAJ-low_{0} MIN-low_{0}  PA-low_{0}    MAJ-hi_{0}  MIN-hi_{0}  PA-hi_{0}    Xpix_{0}   Ypix_{0}   MAXresid_{0}\n'.format(postfix)
		s = ' '.join(s.split())+'\n'
		s = s.replace(' ',',')
		os.system('touch catalogue_SAD_%s.csv' % postfix)
		text_file = open('catalogue_SAD_%s.csv' % postfix,'a')
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
	if os.path.isfile('catalogue_BLOBCAT_%s.csv' % postfix) == False:
		s = 'Catalog_name_{0} rms_{0} BMAJ_{0} BMIN_{0} BPA_{0} ID_{0}  npix_{0} x_p_{0} y_p_{0}           RA_p_{0}  Dec_p_{0}     RA_err_{0}    Dec_err_{0}   x_c_{0}   y_c_{0} RA_c_{0}           Dec_c_{0}      cFlag_{0}   x_wc_{0}   y_wc_{0}   RA_wc_{0} Dec_wc_{0}   wcFlag_{0}       xmin_{0}       xmax_{0}       ymin_{0}       ymax_{0}        rms_{0}    BWScorr_{0} M_{0}    SNR_OBS_{0}     SNR_FIT_{0}         SNR_{0}        S_p_OBS_{0}        S_p_FIT_{0}            S_p_{0}         S_p_CB_{0}      S_p_CBBWS_{0}  S_p_CBBWS_err_{0}      S_int_OBS_{0}    S_int_OBSCB_{0}          S_int_{0}       S_int_CB_{0}   S_int_CB_err_{0}    R_EST_{0}     VisArea_{0}\n'.format(postfix)
		s = ' '.join(s.split())+'\n'
		s = s.replace(' ',',')
		os.system('touch catalogue_BLOBCAT_%s.csv' % postfix)
		text_file = open('catalogue_BLOBCAT_%s.csv' % postfix,'a')
		text_file.write(s)
	for j in files:
		with open(j) as f:
			i = 0
			text_file = open('catalogue_BLOBCAT_%s.csv' % postfix,'a')
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

if useSAD == 'True':
	os.system('rm catalogue_SAD_%s.csv' % postfix)
	print('RUNNING AIPS TASK SAD')
	for file in os.listdir('./'):
		if file.endswith('_casa.fits'):
			fitld = AIPSTask('FITLD')
			hduheader = pyfits.open(file)[0].header
			print(file)
			try:
				data = np.array(pyfits.open(file)[0].data[0,0,edge:edge+rms_box,edge:edge+rms_box])
			except IndexError:
				data = np.array(pyfits.open(file)[0].data[edge:edge+rms_box,edge:edge+rms_box])
			if auto_rms == 'True':
				rms = float(np.sqrt(np.mean(data**2)))
				print(rms)
			fitld.datain = 'PWD:%s' % file
			fitld.outname = str(i)
			fitld.outclass = 'IM'
			fitld.go()
			image = AIPSImage(str(i),'IM',1,1)
			sad = AIPSTask('SAD')
			sad.cparm[1:] = (S_N_ratio+2)*rms, S_N_ratio*rms, (S_N_ratio-1.)*rms
			sad.indata = image
			sad.in2data = image
			sad.blc[1:] = edge,edge
			sad.trc[1:] = int(hduheader['NAXIS1'])-edge,int(hduheader['NAXIS2'])-edge
			sad.dparm[1] = S_N_ratio*rms
			sad.fitout = 'PWD:%s.fitout' % file
			sad.go()
			image.zap()
			lines = open('%s.fitout' % file).readlines()
			if shorthand == 'True':
				names = file[:8]
			else:
				names = file
			try:
				BMAJ = hduheader['BMAJ']
				BMIN = hduheader['BMIN']
				BPA = hduheader['BPA']
			except KeyError:
				print('Run casa_convert.py first to get beam parameters into header')
				sys.exit()
			if len(lines) > 24:
				detections = detections + [file]
				open('%s_r.fitout' % file, 'w').writelines(names+'\n')
				open('%s_r.fitout' % file, 'a').writelines(str(rms)+'\n')
				open('%s_r.fitout' % file, 'a').writelines(str(BMAJ)+'\n')
				open('%s_r.fitout' % file, 'a').writelines(str(BMIN)+'\n')
				open('%s_r.fitout' % file, 'a').writelines(str(BPA)+'\n')
				open('%s_r.fitout' % file, 'a').writelines(lines[18:])
			os.system('rm %s.fitout' % file)
	catalog_list = []
	for file in os.listdir('./'):
		if file.endswith('.fitout'):
			catalog_list = catalog_list + [file]
	if split_catalogues == 'True':
		for k in catalog_list:
			SAD_fit_remove([k],k.split('_casa.fits_r.fitout')[0])
	else:
		SAD_fit_remove(catalog_list,postfix)
	os.system('rm *.fitout')
	print('COMPLETE...')
else:
	os.system('rm catalogue_BLOBCAT_%s.csv' % postfix)
	print('RUNNING BLOBCAT with parameters')
	print('--dSNR=%.2f --fSNR=%.2f --pmep=%.4f --ppe=%.4f --pasbe=%.4f --cpeRA=%.6f --cpeDec=%.6f --edgemin=%d %s %s' % (S_N_ratio,SNR_flood,pmep,ppe,pasbe,cpeRA,cpeDec,int(edge),ds9,write_blobs))
	for file in os.listdir('./'):
		if file.endswith('_casa.fits'):
			print('Cataloguing %s' % file)
			if run_BANE == 'True':
				rms_map = file[:-5]+'_rms.fits'
				os.system('rm %s %s_bkg.fits' % (rms_map,file[:-5]))
				os.system('BANE %s' % file)
			hduheader = pyfits.open(file)[0].header
			try:
				data = np.array(pyfits.open(file)[0].data[0,0,edge:edge+rms_box,edge:edge+rms_box])
			except IndexError:
				data = np.array(pyfits.open(file)[0].data[edge:edge+rms_box,edge:edge+rms_box])
			if auto_rms == 'True':
				rms = float(np.sqrt(np.mean(data**2)))
				print(rms)
			if use_BANE_rms == 'True':
				print('Using BANE rms')
				rms_map = file[:-5]+'_rms.fits'
				os.system("python \"%sblobcat.py\" --dSNR=%.2f --fSNR=%.2f --pmep=%.4f --ppe=%.4f --pasbe=%.4f --cpeRA=%.6f --cpeDec=%.6f --rmsmap=%s --edgemin=%d %s %s %s" % (file_path,S_N_ratio,SNR_flood,pmep,ppe,pasbe,cpeRA,cpeDec,rms_map,int(edge),ds9,write_blobs,file))
			else:
				os.system('python %sblobcat.py --dSNR=%.2f --fSNR=%.2f --pmep=%.4f --ppe=%.4f --pasbe=%.4f --cpeRA=%.6f --cpeDec=%.6f --rmsval=%f --edgemin=%d %s %s %s' % (file_path,S_N_ratio,SNR_flood,pmep,ppe,pasbe,cpeRA,cpeDec,rms,int(edge),ds9,write_blobs,file))
			lines = open('%s_blobs.txt' % file[:-5]).readlines()
			try:
				BMAJ = hduheader['BMAJ'] ## assuming cell is same size on both axes
				BMIN = hduheader['BMIN']
				BPA = hduheader['BPA']
			except KeyError:
				print('Run casa_convert.py first to get beam parameters into header')
				sys.exit()
			if shorthand == 'True':
				names = file[:8]
			else:
				names = file
			if len(lines) > 30:
				detections = detections + [file]
				open('%s_r.blobs' % file, 'w').writelines(names+'\n')
				open('%s_r.blobs' % file, 'a').writelines(str(rms)+'\n')
				open('%s_r.blobs' % file, 'a').writelines(str(BMAJ)+'\n')
				open('%s_r.blobs' % file, 'a').writelines(str(BMIN)+'\n')
				open('%s_r.blobs' % file, 'a').writelines(str(BPA)+'\n')
				open('%s_r.blobs' % file, 'a').writelines(lines[30:])
			elif ds9 == '--ds9':
				os.system('rm %s_ds9.reg' % file[:-5])
			os.system('rm %s_blobs.txt' % file[:-5])
	print('COMPLETE...')
	catalog_list = []
	for file in os.listdir('./'):
		if file.endswith('.blobs'):
			catalog_list = catalog_list + [file]
	if split_catalogues == 'True':
		for k in catalog_list:
			blobcat_fit_remove([k],k.split('_casa.fits_r.blobs')[0])
	else:
		blobcat_fit_remove(catalog_list,postfix)
	os.system('rm *.blobs')

os.system('touch detections.txt')
thefile = open('detections.txt', 'w')
for item in detections:
	thefile.write("%s\n" % item)
