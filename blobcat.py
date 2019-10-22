#!/usr/bin/env python
#
#
#
# ABOUT:
#   BLOBCAT is an application to catalogue blobs (islands comprising individual
#   or multiple blended sources) in a FITS image of total intensity (Stokes I)
#   or linear polarization (L or L_RM). To download this application and some
#   test data to illustrate its use, see http://blobcat.sourceforge.net/ . For
#   reference and intended usage, see README.TXT. If this application helps
#   your work, please cite us. BLOBCAT is released under a BSD 3-Clause
#   License, as follows.
#
#
#
# Copyright (c) 2012, C. A. Hales, T. Murphy, J. R. Curran, E. Middelberg
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of THE UNIVERSITY OF SYDNEY, nor RUHR-UNIVERSITAT
#       BOCHUM, nor the NAMES OF ITS CONTRIBUTORS may be used to endorse or
#       promote products derived from this software without specific prior
#       written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL C. A. HALES, T. MURPHY, J. R. CURRAN, OR E.
# MIDDELBERG BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
#
#
# USAGE: Type the following to get a list of program arguments and an example:
#        ./blobcat.py -h
#
#
#
# ACKNOWLEDGEMENTS:
#   This software was originally developed by C. A. Hales at The University
#   of Sydney. Thanks to Tara Murphy and James Curran who made their core
#   flood fill algorithm used in http://arxiv.org/abs/0708.3092 (Murphy T.
#   et al. 2007, MNRAS, 382, 382) available for use here. Thanks to Enno
#   Middelberg for pointing out compatibility issues with AIPS and for
#   contributing the --plot function.
#
#
#
# HISTORY:
#   1.0  25May2012  Initial version.
#   1.1  11Feb2013  Fixed rad/deg bug in pix2aitwcs; thanks to John Morgan
#                   for spotting this. Corrected DS9 overlay coordinates.
#                   Modified bpaSET user check to prevent hangup at 0deg.
#   1.2  21Jul2014  Fixed PyFits 3.x incompability; thanks to Stefan Blex.
#
#
#



import os
import sys
import time
try:
    import astropy.io.fits as pyfits
except ImportError:
    try:
        print('Using pyfits')
        import pyfits
    except ImportError:
        print('No pyfits or astropy installed... exiting!')
        exit()
from optparse import OptionParser
from scipy.special import erf
import numpy as np
from math import *

# Set some bitwise logic
PEAKED  = 1   # 001
QUEUED  = 2   # 010
VISITED = 4   # 100

# Set code timer variables
t0 = t1 = 0


def starttimer():
    global t0
    t0 = time.time()


def finishtimer():
    global t1
    t1 = time.time()


def timerseconds():
    return (t1 - t0)


def parse_cmds():
    usage = "Usage: %prog surface_brightness_img [options]\n"
    usage += "   eg: %prog sb.fits --rmsmap rms.fits --bwsval 0.92 --write --kvis\n\n"
    usage += "Notes: 1) Input FITS images must have the same pixel dimensions, have their\n"
    usage += "          primary image world coordinate system expressed in equatorial\n"
    usage += "          coordinates (i.e. RA-Dec), and be in either a zenithal equal-area\n"
    usage += "          (ZEA), Hammer-Aitoff (AIT), or (for small fields of view)\n"
    usage += "          slant orthographic (SIN) or North-celestial-pole (NCP) projection.\n"
    usage += "       2) If used, the BWS image should map apparent/true SB (range 0-1)\n"
    usage += "       3) Polarization users: No corrections applied for polarization bias"
    version1 = "%prog"
    version2 = "Version 1.2 (July 2014)"
    version3 = "C. A. Hales <chales@users.sourceforge.net>"
    version4 = "http://blobcat.sourceforge.net/"
    versionFull = '%s\n%s\n%s\n%s\n' % (version1, version2, version3, version4)
    parser = OptionParser(usage=usage, version=versionFull)
    parser.add_option("--rmsval", dest="rmsval", metavar="REAL", 
                      default=False, type="float", 
                      help="use constant rms value (Jy/beam) [default: %default]")
    parser.add_option("--rmsmap", dest="rmsmap", metavar="FILE", 
                      default=False, type="string", 
                      help="read in RMS noise FITS image [default: %default]")
    parser.add_option("--bwsval", dest="bwsval", metavar="REAL", 
                      default=1.0, type="float", 
                      help="use constant bws value (Jy/beam) [default: %default]")
    parser.add_option("--bwsmap", dest="bwsmap", metavar="FILE", 
                      default=False, type="string", 
                      help="read in bandwidth smearing FITS image [default: %default]")
    parser.add_option("--bmaj", dest="bmaj", metavar="REAL", 
                      default=False, type="float", 
                      help="set beam major axis (arcsec) [default: %default]")
    parser.add_option("--bmin", dest="bmin", metavar="REAL", 
                      default=False, type="float", 
                      help="set beam minor axis (arcsec) [default: %default]")
    parser.add_option("--bpa", dest="bpa", metavar="REAL", 
                      default=False, type="float", 
                      help="set beam PA (E from N) (deg) [default: %default]")
    parser.add_option("--dSNR", dest="dSNR", metavar="REAL", 
                      default=5.0, type="float", 
                      help="peak SNR cutoff for accepting blobs [default: %default]")
    parser.add_option("--fSNR", dest="fSNR", metavar="REAL", 
                      default=2.6, type="float", 
                      help="SNR floor to flood down to [default: %default]")
    parser.add_option("--pmep", dest="pmep", metavar="REAL", 
                      default=1.0, type="float", 
                      help="% max estimated pixellation (range 0-1) [default: %default]")
    parser.add_option("--cpeRA", dest="cpeRA", metavar="REAL", 
                      default=0.01, type="float", 
                      help="phase calibrator RA pos. error (arcsec) [default: %default]")
    parser.add_option("--cpeDec", dest="cpeDEC", metavar="REAL", 
                      default=0.01, type="float", 
                      help="phase calibrator Dec pos. error (arcsec) [default: %default]")
    parser.add_option("--SEM", dest="sem", metavar="REAL", 
                      default=0.5, type="float", 
                      help="SELFCAL SEM of gain phases (deg) [default: %default]")
    parser.add_option("--pasbe", dest="pasbe", metavar="REAL", 
                      default=0.05, type="float", 
                      help="% absolute SB error (range 0-1) [default: %default]")
    parser.add_option("--ppe", dest="ppe", metavar="REAL", 
                      default=0.02, type="float", 
                      help="% peak SB pixellation error (range 0-1) [default: %default]")
    parser.add_option("--cb", dest="cb", metavar="REAL", 
                      default=0.e-6, type="float", 
                      help="average clean bias correction (+Jy/beam) [default: %default]")
    parser.add_option("--lamfac", dest="lamfac", metavar="REAL", 
                      default=3.5, type="float", 
                      help="lambda factor for bias correction [default: %default]")
    parser.add_option("--visArea", dest="visArea", action="store_true", default=False,
                      help="compute visibility areas (slows code) [default: %default]")
    parser.add_option("--minpix", dest="minpix", metavar="INT", 
                      default=5, type="int", 
                      help="minimum blob size in pixels [default: %default]")
    parser.add_option("--maxpix", dest="maxpix", metavar="INT", 
                      default=100000, type="int", 
                      help="maximum blob size in pixels [default: %default]")
    parser.add_option("--pixdim", dest="pixdim", metavar="(RA,Dec)", nargs=2, 
                      default=(3, 3), type="int", 
                      help="minimum pixels in RA, Dec dimensions [default: %default]")
    parser.add_option("--edgemin", dest="edgemin", metavar="INT", 
                      default=10, type="int", 
                      help="edge buffer in pixels for extracting blobs [default: %default]")
    parser.add_option("--write", dest="write", action="store_true", default=False,
                      help="write flooded blobs to new FITS image [default: %default]")
    parser.add_option("--hfill", dest="hfill", metavar="REAL", 
                      default=999.0, type="float", 
                      help="output blob highlight value (Jy/beam) [default: %default]")
    parser.add_option("--ds9", dest="ds9", action="store_true", default=False,
                      help="write DS9 annotation file [default: %default]")
    parser.add_option("--kvis", dest="kvis", action="store_true", default=False,
                      help="write KVIS annotation file [default: %default]")
    parser.add_option("--plot", dest="plot", action="store_true", default=False,
                      help="produce screen plot of blobs [default: %default]")
    # In case you want to save your plot as well, use this instead:
    #parser.add_option("--plot", dest="plot", metavar="FILE",
    #                  default=False, type="string",
    #                  help="output plot of blobs (FILE.png) [default: %default]")
    parser.add_option("--pltRng", dest="pltminmax", metavar='(min,max)', nargs=2,
		      default=False, type="float",
		      help="dynamic range of screen plot [default: automagic]")
    (options, args) = parser.parse_args()
    
    if len(args) != 1:
        msg = "Incorrect number of arguments\n"
        msg += "                Rerun with --help for instructions.\n"
        parser.error(msg)
    return (options, args, version2)


def read_data(filename):
    hdulist = pyfits.open(filename)
    hdulist.info()
    data = hdulist[0].data
    header = hdulist[0].header
    keys = header.keys()
    refp = {}
    for key in keys:
        refp[key.lower()] = header[key]
    return (data, header, refp)


def write_data(data, header, filename):
    if header:
        hdu = pyfits.PrimaryHDU(data, header)
    else:
        hdu = pyfits.PrimaryHDU(data)

    if os.path.exists(filename):
        os.remove(filename)
    hdu.writeto(filename)


def sign(x):
    if x >= 0:
        return 1
    else:
        return -1


def pix2coords(pixRA, pixDec, c, p):
    # assume no rotation
    llcos = 1
    llsin = 0
    #  Note that the header values for crpix1 and crpix2 start from 1,
    #+ and that the values coming into this function start from 0.
    #+ So just fool the script by adding 1 to the pixels coming in.
    pixRAm = pixRA + 1
    pixDecm = pixDec + 1
    # get pixRA and pixDec back in degrees
    if p=='ZEA':
      (ra, dec) = pix2zeawcs(pixRAm, pixDecm, c['crpix1'], c['crpix2'],
                                 c['cdelt1'], c['cdelt2'], c['crval1'],
                                 c['crval2'], llcos, llsin)
    elif p=='AIT':
      (ra, dec) = pix2aitwcs(pixRAm, pixDecm, c['crpix1'], c['crpix2'],
                                 c['cdelt1'], c['cdelt2'], c['crval1'],
                                 c['crval2'], llcos, llsin)
    elif p=='NCP':
      (ra, dec) = pix2ncpwcs(pixRAm, pixDecm, c['crpix1'], c['crpix2'],
                                 c['cdelt1'], c['cdelt2'], c['crval1'],
                                 c['crval2'], llcos, llsin)
    elif p=='SIN':
      (ra, dec) = pix2sinwcs(pixRAm, pixDecm, c['crpix1'], c['crpix2'],
                                 c['cdelt1'], c['cdelt2'], c['crval1'],
                                 c['crval2'], llcos, llsin)
    
    return (ra, dec)


def pix2zeawcs(x, y, xpix, ypix, dx, dy, xval, yval, llcos=1, llsin=0):
    # calculates ZEA projection coordinates from pixels
    # NOTE THAT THIS VERSION EXPECTS INPUTS IN DEGREES!
    # See Calabretta & Greisen 2002, A&A, 395, 1077 for details on
    #   how to obtain these equations. In particular, see sections
    #   (in order) eq14, eq15, eq70, table13-ZEA_params, paragraph
    #   under "phi_p=LONPOLEa" p1079, and equations 2.
    # Note that arg(u,v) in C&G2002 = atan(v,u) in python!
    # 
    # the defaults for llcos and llsin are assuming no rotation
    # x,y          input pixel values (in abspix)
    # xpix,ypix    pixel coordinate of reference pixel (crpix)
    # dx,dy        pixel size in degrees (cdelt)
    # xval,yval    coordinate value in degrees at the reference pixel (crval)
    # llcos,llsin  rotation values
    syval = sin(radians(yval))
    cyval = cos(radians(yval))
    
    # Determine phi_p = LONPOLEa
    if yval >= 90.:
      phip = 0.
    else:
      phip = pi
    
    # Determine how far the pixel is from the reference pixel in radians
    # (ie calculate the "intermediate world" or "projection plane" coordinates)
    M = (y - ypix)*radians(dy)
    L = (x - xpix)*radians(dx)
    M = M*llcos + L*llsin
    L = L*llcos - M*llsin
    
    # Determine native spherical coordinates in radians
    phi = atan2( L , -M )
    R = sqrt( M**2 + L**2 )
    theta = pi/2. - 2. * asin(R/2.)
    
    stheta = sin(theta)
    ctheta = cos(theta)
    sphipi = sin(phi - phip)
    cphipi = cos(phi - phip)
    
    # Determine celestial spherical coordinates in degrees
    ra  = xval + degrees(atan2(-ctheta*sphipi , stheta*cyval - ctheta*syval*cphipi))
    dec = degrees(asin((stheta*syval + ctheta*cyval*cphipi)))
    return (ra, dec)


def pix2aitwcs(x, y, xpix, ypix, dx, dy, xval, yval, llcos=1, llsin=0):
    # calculates AIT projection coordinates from pixels
    # NOTE THAT THIS VERSION EXPECTS INPUTS IN DEGREES!
    # See Calabretta & Greisen 2002, A&A, 395, 1077 for details on
    #   how to obtain these equations. In particular, see sections
    #   (in order) eq108, eq106, eq107, table13-ZEA_params, paragraph
    #   under "phi_p=LONPOLEa" p1079, section 2.4, and equations 2.
    # Note that arg(u,v) in C&G2002 = atan(v,u) in python!
    # 
    # Native (long,lat) of the fiducial point for AIT is (0,0)
    #
    # the defaults for llcos and llsin are assuming no rotation
    # x,y          input pixel values (in abspix)
    # xpix,ypix    pixel coordinate of reference pixel (crpix)
    # dx,dy        pixel size in degrees (cdelt)
    # xval,yval    coordinate value in degrees at the reference pixel (crval)
    # llcos,llsin  rotation values
    syval = sin(radians(yval))
    cyval = cos(radians(yval))
    
    # Determine phi_p = LONPOLEa
    if yval >= 0.:
      phip = 0.
    else:
      phip = pi
    
    # Go through all the stuff in section 2.4, yikes...
    deltapP = atan2(0. , cos(phip)) + acos(syval)
    deltapM = atan2(0. , cos(phip)) - acos(syval)
    if fabs(deltapP) <= pi/2. and fabs(deltapM) <= pi/2.:	# point 4
      if yval >= 0.:		# want closest to +pi/2
        if deltapP > deltapM:
	  deltap = deltapP
	else:
	  deltap = deltapM
      else:			# want closest to -pi/2
        if deltapP < deltapM:
	  deltap = deltapP
	else:
	  deltap = deltapM
    elif fabs(deltapP) <= pi/2.:	# point 5 (note also final comment in point 3)
      deltap = deltapP
    elif fabs(deltapM) <= pi/2.:	# point 5 (note also final comment in point 3)
      deltap = deltapM
    
    sdelp = sin(deltap) 
    cdelp = cos(deltap)
    
    if deltap == pi/2.:			# point 2
      alphap = xval + degrees(phip) - 180.
    elif deltap == -pi/2.:		# point 2
      alphap = xval - degrees(phip)
    else:
      alphap = xval
    
    # Determine how far the pixel is from the reference pixel in radians
    # (ie calculate the "intermediate world" or "projection plane" coordinates)
    M = (y - ypix)*radians(dy)
    L = (x - xpix)*radians(dx)
    M = M*llcos + L*llsin
    L = L*llcos - M*llsin
    
    # Determine native spherical coordinates in radians
    Z     = sqrt(1. - (L/4.)**2 - (M/2.)**2)
    phi   = 2.*atan2( Z/2.*L , 2.*Z**2-1. )
    theta = asin(M*Z)
    
    stheta = sin(theta)
    ctheta = cos(theta)
    sppp   = sin(phi - phip)
    cppp   = cos(phi - phip)
    
    # Determine celestial spherical coordinates in degrees
    ra  = alphap + degrees(atan2(-ctheta*sppp , stheta*cdelp - ctheta*sdelp*cppp))
    dec = degrees(asin(stheta*sdelp + ctheta*cdelp*cppp))
    return (ra, dec)


def pix2ncpwcs(x, y, xpix, ypix, dx, dy, xval, yval, llcos=1, llsin=0):
    # calculates NCP projection coordinates from pixels
    # NOTE THAT THIS VERSION EXPECTS INPUTS IN DEGREES!
    # see AIPS memo 27 for details on how to obtain these equations
    # 
    # the defaults for llcos and llsin are assuming no rotation
    # x,y          input pixel values (in abspix)
    # xpix,ypix    pixel coordinate of reference pixel (crpix)
    # dx,dy        pixel size in degrees (cdelt)
    # xval,yval    coordinate value in degrees at the reference pixel (crval)
    # llcos,llsin  rotation values
    syval = sin(radians(yval))
    cyval = cos(radians(yval))
    
    # Determine how far the pixel is from the reference pixel in radians
    # (ie calculate the "intermediate world" or "projection plane" coordinates)
    M = (y - ypix)*radians(dy)
    L = (x - xpix)*radians(dx)
    M = M*llcos + L*llsin
    L = L*llcos - M*llsin
    
    ra  = degrees(atan(L / (cyval - M*syval))) + xval
    dec = sign(yval) * degrees(acos((cyval - M*syval)/cos(radians(ra-xval))))
    return (ra, dec)


def pix2sinwcs(x, y, xpix, ypix, dx, dy, xval, yval, llcos=1, llsin=0):
    # calculates SIN projection coordinates from pixels
    # NOTE THAT THIS VERSION EXPECTS INPUTS IN DEGREES!
    # see AIPS memo 27 for details on how to obtain these equations
    # 
    # the defaults for llcos and llsin are assuming no rotation
    # x,y          input pixel values (in abspix)
    # xpix,ypix    pixel coordinate of reference pixel (crpix)
    # dx,dy        pixel size in degrees (cdelt)
    # xval,yval    coordinate value in degrees at the reference pixel (crval)
    # llcos,llsin  rotation values
    syval = sin(radians(yval))
    cyval = cos(radians(yval))
    
    # Determine how far the pixel is from the reference pixel in radians
    # (ie calculate the "intermediate world" or "projection plane" coordinates)
    M = (y - ypix)*radians(dy)
    L = (x - xpix)*radians(dx)
    M = M*llcos + L*llsin
    L = L*llcos - M*llsin
    
    ra  = degrees(atan(L / (cyval*sqrt(1-L**2-M**2) - M*syval))) + xval
    dec = degrees(asin(M*cyval + syval*sqrt(1-L**2-M**2)))
    return (ra, dec)


def ra2hms(decimalRA):
    r = decimalRA / 15.
    h = floor(r)
    m = floor((r-h)*60.)
    s = (r-h-m/60.)*3600.
    
    mystring = '%i:%i:%.2f' % (h,m,s)
    return (mystring)


def dec2dms(decimalDEC):
    dec = fabs(decimalDEC)
    d = floor(dec)
    m = floor((dec-d)*60.)
    s = (dec-d-m/60.)*3600.
    
    if decimalDEC >= 0:
        mystring = '%i:%i:%.2f' % (d,m,s)
    else:
        mystring = '-%i:%i:%.2f' % (d,m,s)
    return (mystring)


def blob_border(blob):
    # Tara Murphy / James Curran
    x = map(lambda a: a[1], blob)
    min_x = min(x)
    max_x = max(x)
    y = map(lambda a: a[0], blob)
    min_y = min(y)
    max_y = max(y)
    return (min_x, max_x, min_y, max_y)


def explore(data, status, queue, bounds, cutoff, pixel):
    # Tara Murphy / James Curran
    (x, y) = pixel
    if x < 0 or y < 0:
        print '\n Uh-oh. Just found a pixel at coordinate' , pixel
	print 'Something screwy going on, edge masking should have caught this.'
	print '*** Code terminating ***'
	sys.exit()

    if x > 0:
        new = (x - 1, y)
        if not status[new] & QUEUED and data[new] >= cutoff:
	    queue.append(new)
	    status[new] |= QUEUED

    if x < bounds[0]:
        new = (x + 1, y)
        if not status[new] & QUEUED and data[new] >= cutoff:
	    queue.append(new)
	    status[new] |= QUEUED

    if y > 0:
        new = (x, y - 1)
        if not status[new] & QUEUED and data[new] >= cutoff:
	    queue.append(new)
	    status[new] |= QUEUED

    if y < bounds[1]:
        new = (x, y + 1)
        if not status[new] & QUEUED and data[new] >= cutoff:
	    queue.append(new)
	    status[new] |= QUEUED


def flood(data, status, bounds, peak, cutoff):
    # Tara Murphy / James Curran
    if status[peak] & VISITED:
        return ([], 0)

    flux = 0.
    blob = []
    queue = [peak]
    status[peak] |= QUEUED

    for pixel in queue:
        if status[pixel] & VISITED:
            continue
    
        status[pixel] |= VISITED

        blob.append(pixel)
        flux += dataimg[pixel]
        explore(data, status, queue, bounds, cutoff, pixel)

    return (blob, flux)


def beamvol(c):
    # Omega_b, Equation 4
    vol = pi/(4.*log(2.)) * c['bmaj'] * c['bmin'] / fabs(c['cdelt1'] * c['cdelt2'])
    return (vol)


def blobarea(cut,p,rms,bws,n,c):
    # R_EST, Equation 34
    # cut = cutoff SNR
    # p   = peak bias corrected peak SB of source (no CB or BWS correction)
    # rms = local rms at peak
    # bws = bandwidth smearing (attenuated/true peak value)
    # n   = number of pixels in blob
    # c   = header items
    beamsize = c['bmaj'] * c['bmin']
    beamsizeLOCAL = beamsize / bws
    beams = n / (pi/4. * log(p/(cut*rms),2.) * beamsizeLOCAL / fabs(c['cdelt1']*c['cdelt2']) )
    
    return (beams)


def getintf(minSNR,maxSNR,floodf):
    # Equation 17
    totf = floodf / (erf(sqrt(-log(minSNR/maxSNR))))**2
    
    return (totf)


def genparabCM():
    # Generate 9x6 dummy coefficient matrix for parabfit
    # Do this here to prevent having to make it for each blob
    
    global pfCM
    pfCM = np.zeros((9,6))
    
    k = 0
    for i in range(1,4):
      y = i-1
      for j in range(1,4):
	x = j-1
	pfCM[k,0] = 1.
	pfCM[k,1] = x
	pfCM[k,2] = y
	pfCM[k,3] = x*y
	pfCM[k,4] = x*x
	pfCM[k,5] = y*y
        k += 1


def parabfit(f):
    # Appendix A
    # Returns peak of 2D parabolic fit from 3x3 array about observed peak
    # Use linear least squares to solve for an overdetermined
    # (m>n) system of linear equations; here, m=9 is the number
    # of equations, and n=6 is the number of unknowns.
    #
    # 2D parabola = c0 + c1*x + c2*y + c3*x*y + c4*x**2 + c5*y**2
    # 
    # f = 9 element vector of pixel values in 3x3 array about peak
    # pfCM = 9x6 coefficient matrix, see genparabCM()
    
    # Get c coefficients
    c = np.linalg.lstsq(pfCM,f)[0]
    
    # To get the position of the peak, take df/dx=df/dy=0
    # 2 equations, 2 unknowns, can solve analytically
    xp = ( c[1] - 0.5*c[2]*c[3]/c[5] ) / ( 0.5*c[3]**2/c[5] - 2.*c[4] )
    yp = ( -c[1] - 2.*c[4]*xp ) / c[3]
    
    # FYI, the fitted peak in coordinates you know and love are:
    # xpix = peak[1] + xp - 1
    # ypix = peak[0] + yp - 1
    
    # now get the value of the parabolic fit at the fitted peak
    fittedpeak = c[0] + c[1]*xp + c[2]*yp + c[3]*xp*yp + c[4]*xp**2 + c[5]*yp**2
    
    return (fittedpeak)


def poserr(sigcal,sigfr,pcsnr,w,c):
    # Positional uncertainty - Equations 18/22
    #
    # sigcal = positional error in secondary calibrator location (asec)
    # sigfr  = frame registration error (deg)
    # pcsnr  = peak-corrected SNR for blob
    # w      = 1 for projected resolution along RA axis;
    #        = 2 for projected resolution along Dec axis
    # c      = header items
    #
    # beam position angle assumed east from north
    if w == 1:
	# ie when PA=0, this returns the BMIN
	b = c['bmaj'] * c['bmin'] / sqrt( ( c['bmaj']*cos(radians(c['bpa'])) )**2 + \
	                                  ( c['bmin']*sin(radians(c['bpa'])) )**2 )
    else:
	# ie when PA=0, this returns the BMAJ
	b = c['bmaj'] * c['bmin'] / sqrt( ( c['bmaj']*sin(radians(c['bpa'])) )**2 + \
	                                  ( c['bmin']*cos(radians(c['bpa'])) )**2 )
    
    # factor from equation 37
    errfac = 1.4
    
    # return error in deg
    perr = sqrt((sigcal/3600.)**2 + (sigfr/180./sqrt(2.)*b)**2 + (b/(errfac*pcsnr))**2)
    
    return (perr)


def genlookup():
    # populates lookup table for correctpeak function below
    # Equation 14
    
    global cpM , cpA
    
    a1 = 0.89
    a2 = 0.27
    a3 = 3.75
    a4 = -3.67
    a5 = 1.61
    
    smin = 0.
    smax = 2.
    step = 0.01
    
    # populate table for M = 1 --> ~30
    k = int((smax-smin)/step) + 1
    # cpM = M, the known x-value in correctpeak()
    # cpA = E[snr], the unknown y-value in correctpeak()
    cpM = np.zeros(k)
    cpA = np.zeros(k)
    
    # Note that range(a,b) goes from a --> b-1
    for i in range(0,k):
        cpA[i] = smin + i*step
	cpM[i] = 1. + a1*cpA[i] + a2*(cpA[i])**2 + a3*(cpA[i])**3 + \
	                          a4*(cpA[i])**4 + a5*(cpA[i])**5


def correctpeak(Aobs,m):
    # Equation 15
    # Aobs = observed snr at peak
    # m    = source extent in beams at rawpeak-lambda*sig level
    #
    # if m >= 1.1, subtract expectation value of positive bias
    #+ from lookup table
    if m >= 1.1:
        A = Aobs - cpA[abs(cpM-m).argmin()]
    else:
        A = Aobs
    
    return (A)






# main program
starttimer()
print ''

# get command line parameters
(options, args, version) = parse_cmds()

FILENAME = args[0]
if not os.path.exists(FILENAME):
  print '\n Error: FITS file "%s" does not exist\n' % FILENAME
  print '*** Code terminating ***'
  sys.exit()

RMSFILENAME = options.rmsmap
rmsV        = options.rmsval
if RMSFILENAME:
  if not os.path.exists(RMSFILENAME):
    print '\n Error: RMS noise FITS file "%s" does not exist\n' % RMSFILENAME
    print '*** Code terminating ***'
    sys.exit()
else:
  if not rmsV:
    print '\n Error: You need to specify a constant rms noise value or supply a FITS image\n'
    print '*** Code terminating ***'
    sys.exit()

#  If a bws fits image is specified, values are taken from there, otherwise take
#+ the user specified bws value or, failing that, the default value of 1.0
bwsV        = options.bwsval
BWSFILENAME = options.bwsmap
if BWSFILENAME:
  if not os.path.exists(BWSFILENAME):
    print '\n Error: Bandwidth smearing FITS file "%s" does not exist' % BWSFILENAME
    print '        You need to specify a constant bandwidth smearing value or supply a FITS image\n'
    print '*** Code terminating ***'
    sys.exit()
else:
  if not 0. <= bwsV <= 1.:
    print '\n Error: Bandwidth smearing value must be in range 0-1\n'
    print '*** Code terminating ***'
    sys.exit()


OUTFILE      = os.path.splitext(FILENAME)[0]
FLOOR_CUTOFF = options.fSNR
PEAK_CUTOFF  = options.dSNR
PEAK_DETECT  = PEAK_CUTOFF * (1. - options.pmep)
if PEAK_DETECT < FLOOR_CUTOFF:
  PEAK_DETECT  = options.fSNR
  print ' Warning: You\'ve set your max estimated pixellation error'
  print '          (pmep) too high. Modifying 1st pass peak SNR'
  print '          detection threshold to be no less than the'
  print '          flooding threshold. Don\'t worry, the worst'
  print '          that can happen is that this program takes'
  print '          a little bit longer to run.'

MIN_PIX      = options.minpix
MAX_PIX      = options.maxpix
PIX_DIM      = options.pixdim
EDGE         = options.edgemin
calRAerr     = options.cpeRA
calDECerr    = options.cpeDEC
semerr       = options.sem
abserr       = options.pasbe
pixerr       = options.ppe
CBcorr       = options.cb
lamval       = options.lamfac
BLOBBED      = options.hfill
bmajSET      = options.bmaj
bminSET      = options.bmin
bpaSET       = options.bpa
pltRng       = options.pltminmax
#  to save output plot, modify --plot, uncomment save line at
#+ the very bottom of this code, and uncomment line below
#pltFile = OUTFILE + '.png'




if CBcorr < 0:
  print '\n Error: Clean bias correction must be positive\n'
  print '*** Code terminating ***'
  sys.exit()

print '\nReading FITS image data...'
(dataimg, header, refp) = read_data(FILENAME)

# get image projection (assume RA and Dec are in same projection - they'd better be!)
proj=refp['ctype1'][-3:].upper()
# check if this projection is supported here
projVals=['ZEA','AIT','NCP','SIN']
if not proj in projVals:
  print '\n Error: Image projections must be ZEA, AIT, NCP, or SIN.\n'
  print '*** Code terminating ***'
  sys.exit()
else:
  print '\nUsing %s projection.' % proj
  if proj=='NCP' or proj=='SIN':
    print 'Warning: You\'d better not be working too far from the reference point with %s projection.' % proj

if proj=='AIT' and options.kvis:
  print '\n Warning: You can\'t view AIT projection in kvis.'
  print '          Turning off your kvis annotation request.'
  options.kvis = False

# get coordinate system epoch
if 'equinox' in refp.keys():
  epoch=refp['equinox']
elif 'epoch' in refp.keys():
  epoch=refp['epoch']
else:
  epoch=-1
  print '\n Warning: Unknown equinox for coordinate system.'

# Warn user that no attempt has been made to account for weird coordinate system rotations.
# Don't try to assume what is going on - make user manually update code to do what they want.
# The routines pix2coords and pix2###wcs will need to be modified.
if 'llrot' in refp.keys():
  # pyfits should have converted this to crota2, but run check just in case
  if not refp['llrot']==0:
    print '\n Error: Coordinate rotation not supported - please modify the code (some provision provided).\n'
    print '*** Code terminating ***'
    sys.exit()

if 'crota1' in refp.keys():
  if not refp['crota1']==0:
    print '\n Error: Coordinate rotation not supported - please modify the code (some provision provided).\n'
    print '*** Code terminating ***'
    sys.exit()

if 'crota2' in refp.keys():
  if not refp['crota2']==0:
    print '\n Error: Coordinate rotation not supported - please modify the code (some provision provided).\n'
    print '*** Code terminating ***'
    sys.exit()

# check if beam keywords exist in refp
if not ('bmaj' in refp.keys() and 'bmin' in refp.keys() and 'bpa' in refp.keys()):
  #  bmaj/bmin/bpa are not standardised FITS header keywords, so
  #+ check the 'HISTORY' card for this info to satisfy AIPS users...
  print '\n Beam information in header incomplete - trying to guess it from history card'
  # loop over header items and look for 'bmaj'
  for x in header.items():
    # extract bmaj, bmin and bpa from history info
    if type(x[1])==str and 'bmaj' in x[1].lower():
      try:
	line=x[1].lower().replace('=', '').split()
	refp['bmaj']=float(line[line.index('bmaj')+1])
	refp['bmin']=float(line[line.index('bmin')+1])
	refp['bpa']=float(line[line.index('bpa')+1])
      except:
	print ' Couldn\'t extract beam info from file. Now checking if you specified bmaj/bmin/bpa explicitly...'

# overwrite with user-specified values if required
if bmajSET:
  refp['bmaj'] = bmajSET/3600.

if bminSET:
  refp['bmin'] = bminSET/3600.

if not bpaSET is False:
  refp['bpa'] = bpaSET

if 'bmaj' in refp.keys():
  print '\nUsing bmaj = %.3e arcsec' % (refp['bmaj']*3600.0)
else:
  print '\n Error: Beam major axis could not be found in FITS image header.'
  print '        You need to set --bmaj manually.\n'
  print '*** Code terminating ***'
  sys.exit()

if 'bmin' in refp.keys():
  print 'Using bmin = %.3e arcsec' % (refp['bmin']*3600.0)
else:
  print '\n Error: Beam minor axis could not be found in FITS image header.'
  print '        You need to set --bmin manually.\n'
  print '*** Code terminating ***'
  sys.exit()

if 'bpa' in refp.keys():
  print 'Using bpa  = %.3e deg' % (refp['bpa'])
else:
  print '\n Error: Beam position angle could not be found in FITS image header.'
  print '        You need to set --bpa manually.\n'
  print '*** Code terminating ***'
  sys.exit()

#  Assume that x is the second-last axis and y is the last axis, independent
#+ of how many axes there are in dataimg when it comes back from read_data.
#+ Working assumption here is that the the slower-varying axes always come
#+ first in a FITS file, so the image pixels should be the last two axes.
dataimg = dataimg.reshape(dataimg.shape[-2],dataimg.shape[-1])
status  = np.zeros(dataimg.shape, dtype=np.uint8)
bounds  = (dataimg.shape[0] - 1, dataimg.shape[1] - 1)

if RMSFILENAME:
  print ' '
  print 'Reading FITS rms data...'
  (datarms, unused1, refpRMS) = read_data(RMSFILENAME)
  datarms = datarms.reshape(datarms.shape[-2],datarms.shape[-1])
  # check RMS projection, just in case...
  projRMS=refpRMS['ctype1'][-3:].upper()
  if not proj==projRMS:
    print '\n Error: RMS image projection inconsistent with SB image.\n'
    print '*** Code terminating ***'
    sys.exit()
else:
  datarms = dataimg*0+rmsV

# already checked that bws file exists
if BWSFILENAME:
  print ' '
  print 'Reading FITS bandwidth smearing data...'
  (dataBWS, unused1, refpBWS) = read_data(BWSFILENAME)
  dataBWS = dataBWS.reshape(dataBWS.shape[-2],dataBWS.shape[-1])
  # check BWS image projection, just in case...
  projBWS=refpBWS['ctype1'][-3:].upper()
  if not proj==projBWS:
    print '\n Error: BWS image projection inconsistent with SB image.\n'
    print '*** Code terminating ***'
    sys.exit()
else:
  dataBWS = dataimg*0+bwsV

print ' '
print 'Making signal-to-noise map...'
datasnr = dataimg.__div__(datarms)

print ' '
print 'Finding pixels above peak SNR cutoff...'
peaks = np.where(datasnr >= PEAK_DETECT)
status[peaks] = PEAKED

# reform the peaks as a list of tuples
peaks = map(lambda x, y: (datasnr[x,y], x, y), peaks[0], peaks[1])
peaks.sort(reverse=True)
peaks = map(lambda x: x[1:], peaks)
# remove peaks within EDGE pixels of border
EDGE = EDGE - 1
peaks = filter(lambda x: x[0] > EDGE and x[0] < bounds[0] - EDGE and
                         x[1] > EDGE and x[1] < bounds[1] - EDGE, peaks)

# get beam volume conversion factor
beamconv = beamvol(refp)

# calculate number of non-blank pixels in image
totpix = len(np.where(datarms >= 0.)[0])

# generate lookup table for Atrue (comes out of function correctpeak)
genlookup()

# populate dummy coefficient matrix for 2D parabolic fitting
genparabCM()

print ' '
print 'Number of pixels above peak SNR cutoff and within valid region =', len(peaks)

if options.visArea:
    print ' '
    print 'Making bandwidth smearing corrected rms map for visibility area correction...'
    datarmsBS = datarms.__div__(dataBWS)
    print ' '

pos = open('%s_blobs.txt' % OUTFILE, 'w')
print >> pos, '# %s' % os.path.basename(sys.argv[0])
print >> pos, '# %s\n#' % version
print >> pos, '# SB image = %s' % FILENAME
if RMSFILENAME:
  print >> pos, '# RMS image = %s' % RMSFILENAME
else:
  print >> pos, '# RMS (constant) = %.3e' % rmsV
if BWSFILENAME:
  print >> pos, '# Bandwidth smearing image = %s' % BWSFILENAME
else:
  print >> pos, '# Bandwidth smearing (constant) = %.3e' % bwsV
print >> pos, '# Total image size (including blanks) along dimensions (RA , Dec) = (%d , %d) pixels' % \
               (datasnr.shape[1], datasnr.shape[0])
print >> pos, '# Image projection = %s' % proj
print >> pos, '# Total non-blank pixels = %d' % totpix
print >> pos, '# Pixel dimensions (RA , Dec) = (%.3e , %.3e) arcsec' % (abs(refp['cdelt1']*3600), abs(refp['cdelt2']*3600))
print >> pos, '# Resolution (bmaj , bmin) = (%.3e , %.3e) arcsec, PA = %.3e deg' % \
               (refp['bmaj']*3600, refp['bmin']*3600, refp['bpa'])
print >> pos, '# SNR blob detection threshold (T_d) = %.3e' % PEAK_CUTOFF
print >> pos, '# SNR blob flooding threshold (T_f)  = %.3e\n#' % FLOOR_CUTOFF
print >> pos, '# Phase calibrator positional error (RA , Dec) = (%.3e , %.3e) arcsec' % (calRAerr, calDECerr)
print >> pos, '# SELFCAL standard error of the mean (SEM) of gain phases = %.3e deg' % semerr
print >> pos, '# Absolute surface brightness error = %.2f %%' % (abserr*100.)
print >> pos, '# Pixellation error (affects peak errors) = %.2f %%' % (pixerr*100.)
print >> pos, '# Clean bias correction = %.3e Jy/beam' % CBcorr
print >> pos, '# Lambda factor for resolved blob peak bias correction = %.2f' % lamval
print >> pos, '# Blob size limits (min , max) = (%d , %d) pixels' % (MIN_PIX, MAX_PIX)
print >> pos, '# Blob size limits along dimensions (RA , Dec) = (%d , %d) pixels' % PIX_DIM
print >> pos, '# Edge buffer for extracting blobs = %d pixels\n#' % (EDGE+1)
print >> pos, '# Important: pixel coordinates start from 0, not 1.'
if epoch == -1:
  print >> pos, '# Equinox (in years) for equatorial coordinate system = UNKNOWN\n#'
else:
  print >> pos, '# Equinox (in years) for equatorial coordinate system = %.2f\n#' % epoch
print >> pos, '#     C1         C2         C3         C4              C5              C6         C7         C8  ',
print >> pos, '          C9            C10             C11             C12        C13  ',
print >> pos, '         C14            C15             C16             C17        C18  ',
print >> pos, '     C19        C20        C21        C22        C23        C24      C25         C26         C27     ',
print >> pos, '   C28            C29            C30            C31            C32            C33            C34     ',
print >> pos, '      C35            C36            C37            C38            C39      C40         C41'
print >> pos, '#     ID       npix        x_p        y_p            RA_p           Dec_p     RA_err    Dec_err  ',
print >> pos, '         x_c            y_c            RA_c           Dec_c      cFlag  ',
print >> pos, '        x_wc           y_wc           RA_wc          Dec_wc     wcFlag  ',
print >> pos, '    xmin       xmax       ymin       ymax        rms    BWScorr        M     SNR_OBS     SNR_FIT     ',
print >> pos, '   SNR        S_p_OBS        S_p_FIT            S_p         S_p_CB      S_p_CBBWS  S_p_CBBWS_err     ',
print >> pos, 'S_int_OBS    S_int_OBSCB          S_int       S_int_CB   S_int_CB_err    R_EST     VisArea'
print >> pos, '#             [pix]      [pix]      [pix]           [deg]           [deg]      [deg]      [deg]  ',
print >> pos, '       [pix]          [pix]           [deg]           [deg]             ',
print >> pos, '       [pix]          [pix]           [deg]           [deg]             ',
print >> pos, '   [pix]      [pix]      [pix]      [pix]    [ Jy/b]                                                 ',
print >> pos, '               [Jy/b]         [Jy/b]         [Jy/b]         [Jy/b]         [Jy/b]         [Jy/b]     ',
print >> pos, '     [Jy]           [Jy]           [Jy]           [Jy]           [Jy]'

if options.ds9:
    annD = open('%s_ds9.reg' % OUTFILE, 'w')
    print >> annD, 'global color=green font=\"helvetica 12 normal\" edit=1 move=1 delete=1 include=1'
    print >> annD, 'image'

if options.kvis:
    annK = open('%s_kvis.ann' % OUTFILE, 'w')
    print >> annK, 'COLOUR red\nCOORD W\nFONT hershey12\nPA sky'

print '--------------------------------'
print 'N peak_loc = npixels'
nblobs = 0
nfilt_size = 0
for peak in peaks:
    (blob, flux) = flood(datasnr, status, bounds, peak, FLOOR_CUTOFF)
    npixels = len(blob)
    # filter on blob size
    if MIN_PIX <= npixels <= MAX_PIX:
        border = blob_border(blob)
        #  Recall that python stores the row (Dec) in x. However, to throw
	#+ a spanner in the works, I have set blob_border to return things
	#+ in the correct order, so things coming out of border are in
	#+ order of: RAmin RAmax Decmin Decmax
        npixRA = border[1] - border[0]
	npixDec = border[3] - border[2]

        # filter on blob dimension
        if npixRA >= PIX_DIM[0] and npixDec >= PIX_DIM[1]:
            
	    #  only attempt to get fitted peak if the pixel in
	    #+ question is the maximum within a 3x3 array
	    fitarray = np.array(dataimg[peak[0]-1:peak[0]+2,peak[1]-1:peak[1]+2].reshape(9))
	    if dataimg[peak] >= fitarray.max():
	        # get fitted peak SB using 2D parabola to 3x3 array about peak
	        peakFIT = parabfit(fitarray)
	        
	        #  if peakFIT is smaller than peak (which is possible for small values
	        #+ of N; see Fig. A1 of the BLOBCAT manuscript), then for consistency, choose
	        #+ to accept the observed peak as the fitted peak
	        if dataimg[peak] > peakFIT:
	            peakFIT = dataimg[peak]
	    else:
	        peakFIT = dataimg[peak]
	    
	    # blobs satisfying PEAK_DETECT may not satisfy PEAK_CUTOFF
	    # filter on PEAK_CUTOFF using the fitted peak SB
	    if peakFIT/datarms[peak] >= PEAK_CUTOFF:
		
		# area of source at lambda*sig below observed raw peak
		pixcount = 0
		for pixel in blob:
		    if datasnr[pixel] >= datasnr[peak]-lamval:
		        pixcount += 1
		
		# number of independent beams
		Mval = pixcount / (sqrt(12)/4/log(2.) * refp['bmaj'] * refp['bmin'] / 
		                   fabs(refp['cdelt1']*refp['cdelt2']) )
		
		# get peak bias corrected SNR
		Aval = correctpeak(peakFIT/datarms[peak],Mval)
		
		# Aval may now be smaller than PEAK_CUTOFF
		# filter on PEAK_CUTOFF again...
		if Aval >= PEAK_CUTOFF:
		    
		    # print out results
		    nblobs += 1
		    print nblobs, '(', peak[1],',',peak[0], ') =', npixels, 'pixels'
		    
		    (peakra, peakdec) = pix2coords(peak[1], peak[0], refp, proj)
		    #peakraHMS  = ra2hms(peakra)
		    #peakdecDMS = dec2dms(peakdec)
		    
		    # get unweighted centroid position
		    peakC = np.average(blob, axis=0)
		    (peakraC, peakdecC) = pix2coords(peakC[1], peakC[0], refp, proj)
		    # is unweighted centroid located in a flooded pixel?
		    if (floor(peakC[0]),floor(peakC[1])) in blob:
		        ctFlag = 1
		    else:
		        ctFlag = 0
		    
		    # get weighted centroid position
		    wts = datasnr[np.array(blob)[:,0],np.array(blob)[:,1]]
		    peakWC = np.average(blob, axis=0, weights=wts)
		    (peakraWC, peakdecWC) = pix2coords(peakWC[1], peakWC[0], refp, proj)
		    # is weighted centroid located in a flooded pixel?
		    if (floor(peakWC[0]),floor(peakWC[1])) in blob:
		        wctFlag = 1
		    else:
		        wctFlag = 0
		    
		    # get peak bias corrected peak SB
		    peakCORR = Aval*datarms[peak]
		    
		    # calculate positional errors
		    RAerr  = poserr(calRAerr,semerr,Aval,1,refp)
		    DECerr = poserr(calDECerr,semerr,Aval,2,refp)
		    
		    # Bottom left going clockwise
		    (ra1, dec1) = pix2coords(border[0], border[2], refp, proj)
		    (ra2, dec2) = pix2coords(border[1], border[2], refp, proj)
		    (ra3, dec3) = pix2coords(border[1], border[3], refp, proj)
		    (ra4, dec4) = pix2coords(border[0], border[3], refp, proj)
		    
		    # clean bias correction
		    peakCB = peakCORR + CBcorr
		    fluxCB = flux + npixels * CBcorr
		    
		    # get flux density in Jy, not surface brightness in Jy/beam
		    flux   = flux / beamconv
		    fluxCB = fluxCB / beamconv
		    
		    # correct flux based on Gaussian morphology estimate of what
		    #+ has been missed below the SNR flood cutoff
		    fluxCORR   = getintf(FLOOR_CUTOFF,Aval,flux)
		    fluxCBCORR = getintf(FLOOR_CUTOFF,peakCB/datarms[peak],fluxCB)
		    
		    # apply bandwidth smearing correction to peak SB
		    BWScorr   = 1 / dataBWS[peak]
		    peakCBBWS = peakCB * BWScorr
		    
		    # estimate measurement errors
		    peakerr = sqrt((abserr*peakCBBWS)**2+(pixerr*peakCBBWS)**2+(BWScorr*datarms[peak])**2)
		    fluxerr = sqrt((abserr*fluxCBCORR)**2+(datarms[peak])**2)
		    
		    # diagnostic
		    Rest = blobarea(FLOOR_CUTOFF,peakCORR,datarms[peak],dataBWS[peak],npixels,refp)
		    
		    #  visibility area (observable source density), taking into account both
		    #+ bandwidth smearing and image sensitivity
		    if options.visArea:
		        visA = (float(len(np.where(datarmsBS <= (peakCORR*BWScorr)/PEAK_CUTOFF)[0])) /
		                float(totpix))
		    else:
		        visA = -1.
		    
		    print >> pos, '%8d' % nblobs,
		    print >> pos, '%10d' % npixels,
		    print >> pos, '%10d' % peak[1],
		    print >> pos, '%10d' % peak[0],
		    print >> pos, '%15.9f' % peakra,
		    print >> pos, '%15.9f' % peakdec,
		    #print >> pos, '%12s' % peakraHMS,	# fiddle here if you want RA in HMS
		    #print >> pos, '%14s' % peakdecDMS,	# fiddle here if you want Dec in DMS
		    print >> pos, '%10.2e' % RAerr,
		    print >> pos, '%10.2e' % DECerr,
		    print >> pos, '%14.2f' % peakC[1],
		    print >> pos, '%14.2f' % peakC[0],
		    print >> pos, '%15.9f' % peakraC,
		    print >> pos, '%15.9f' % peakdecC,
		    print >> pos, '%10d' % ctFlag,
		    print >> pos, '%14.2f' % peakWC[1],
		    print >> pos, '%14.2f' % peakWC[0],
		    print >> pos, '%15.9f' % peakraWC,
		    print >> pos, '%15.9f' % peakdecWC,
		    print >> pos, '%10d' % wctFlag,
		    print >> pos, '%10d %10d %10d %10d' % blob_border(blob),
		    print >> pos, '%10.2e' % datarms[peak],
		    print >> pos, '%10.2e' % BWScorr,
		    print >> pos, '%8.2f' % Mval,
		    print >> pos, '%11.3e' % (dataimg[peak]/datarms[peak]),
		    print >> pos, '%11.3e' % (peakFIT/datarms[peak]),
		    print >> pos, '%11.3e' % Aval,
		    print >> pos, '%14.6e' % dataimg[peak],
		    print >> pos, '%14.6e' % peakFIT,
		    print >> pos, '%14.6e' % peakCORR,
		    print >> pos, '%14.6e' % peakCB,
		    print >> pos, '%14.6e' % peakCBBWS,
		    print >> pos, '%14.6e' % peakerr,
		    print >> pos, '%14.6e' % flux,
		    print >> pos, '%14.6e' % fluxCB,
		    print >> pos, '%14.6e' % fluxCORR,
		    print >> pos, '%14.6e' % fluxCBCORR,
		    print >> pos, '%14.6e' % fluxerr,
		    print >> pos, '%8.2f' % Rest,
		    print >> pos, '%11.3e' % visA
		    
		    (ra1, dec1) = pix2coords(border[0], border[2], refp, proj)
		    (ra2, dec2) = pix2coords(border[1], border[2], refp, proj)
		    (ra3, dec3) = pix2coords(border[1], border[3], refp, proj)
		    (ra4, dec4) = pix2coords(border[0], border[3], refp, proj)
		    
		    # Ain't nothing l33t about the annotation files...just make boxes
		    if options.ds9:
		        # ds9 pixels start from 1, not 0 like kvis
		        print >> annD, "polygon %d %d %d %d %d %d %d %d %d %d # text={%d}" \
				% (border[0]+1, border[2]+1, border[1]+1, border[2]+1, \
				   border[1]+1, border[3]+1, border[0]+1, border[3]+1, \
				   border[0]+1, border[2]+1, nblobs)
		    
		    if options.kvis:
		        print >> annK, "TEXT %.9f %.9f %d" % (peakra, peakdec, nblobs)
		        print >> annK, "CLINES %f %f %f %f %f %f %f %f %f %f" \
				% (ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, ra1, dec1)
		    
		    for pixel in blob:
		        dataimg[pixel] = BLOBBED
		
		else:
		    print 'Blob has peak below SNR cutoff due to peak bias correction. Filtering.'
		    nfilt_size += 1
            
	    # No need to accumulate filtering statistic for this loop.
	
	elif datasnr[peak] > PEAK_CUTOFF:
	    #  only print those with SNR greater than PEAK_CUTOFF
	    #+ so that pixellation rejections are not counted.
	    print 'Flooded blob dimensions too small. Filtering.'
	    nfilt_size += 1
    
    elif npixels > MAX_PIX and datasnr[peak] > PEAK_CUTOFF:
	#  only print those with SNR greater than PEAK_CUTOFF
	#+ so that pixellation rejections are not counted.
        print 'Flooded blob npixels above max range. Filtering.'
        nfilt_size += 1
    
    elif 0 < npixels < MIN_PIX and datasnr[peak] > PEAK_CUTOFF:
        #  blobs with zero pixels should not be counted, nor
	#+ should those with SNR greater than PEAK_CUTOFF
	#+ so that pixellation rejections are not counted.
        print 'Flooded blob npixels below min range. Filtering.'
        nfilt_size += 1

print '--------------------------------'
print 'number of extracted blobs =', nblobs
print 'number of filtered blobs =', nfilt_size

pos.close()
if options.ds9:
    print '\nDS9 annotation file written...'
    annD.close()

if options.kvis:
    print '\nkvis annotation file written...'
    annK.close()

if options.write:
    print '\nWriting output FITS file...'
    write_data(dataimg, header, '%s_blobs.fits' % OUTFILE)

finishtimer()

print '\nCode ran in %.6f seconds.' % timerseconds()
print ''

# make screen plot if desired, first check availability of relevant packages
if options.plot:
    # Enno Middelberg
    print 'Making diagnostic plot...'
    try:
        import matplotlib.pyplot as plt
        # switch to interactive mode
        plt.ion()
        import aplpy
    except(ImportError):
        print 'Can\'t find package matplotlib and/or aplpy - please install to enable plotting.\n'
        print '*** Code terminating ***'
        sys.exit()
    #
    # open up the input fits file and plot it as greyscale image
    hdu = pyfits.open(FILENAME)[0]
    # check if pltRng has been specified on command line and act accordingly
    if pltRng:
	fig=aplpy.FITSFigure(hdu)
	fig.show_grayscale(vmin=pltRng[0], vmax=pltRng[1], interpolation='nearest')
    else:
    	# figure out the range automagically
	fig=aplpy.FITSFigure(hdu)
	fig.show_grayscale(interpolation='nearest')
    # replace data in input fits file with floodfilled data and make a contour plot
    hdu.data=dataimg
    fig.show_contour(hdu, levels=[0.999*BLOBBED], colors='red')
    fig.add_beam(major=refp['bmaj'], minor=refp['bmin'], angle=refp['bpa'])
    #  save a copy of the image? you'll need to change the way --plot works;
    #+ you'll need to modify cmd line options and uncomment pltFile
    #fig.save(pltFile)
    raw_input('Press a key to exit.')
    fig.close()
