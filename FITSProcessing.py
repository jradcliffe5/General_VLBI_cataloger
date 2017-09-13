# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 14:33:35 2015
# Usage is Python FITSProcessing.py filename imagesize
@author: Jack
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
import aplpy
import os.path
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

matplotlib.pyplot.close("all")
filename = str(sys.argv[1])
pixsiz = int(sys.argv[2])
outfilename = filename[:-5]+'Py.fits'

def convertAIPStoPythonImage(filename,outfilename):
    hdu_list = fits.open(filename)

    head = hdu_list['PRIMARY'].header
    head['CRVAL1'] = 360 - (head['CRVAL1']*-1)
    head['NAXIS'] = 2
    del head['NAXIS3']
    del head['NAXIS4'] 
    del head['CTYPE3'],  head['CRVAL3'], head['CDELT3'], head['CRPIX3'], head['CROTA3'] 
    del head['CTYPE4'], head['CRVAL4'], head['CDELT4'], head['CRPIX4'], head['CROTA4']
        #if filename.endswith('NA.fits') == False:
        #head.set('BMIN', float(hdu_list[1].data.field('BMIN')))
        #head.set('BMAJ', float(hdu_list[1].data.field('BMAJ')))
        #head.set('BPA', float(hdu_list[1].data.field('BPA')))
    image_data = hdu_list[0].data
    image_data = image_data[0,0,:,:]
    hdu = fits.PrimaryHDU(image_data, header=head)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(outfilename)
    return outfilename

if os.path.isfile(outfilename) == False: 
    hdu_list = fits.open(convertAIPStoPythonImage(filename,outfilename))
else:
    hdu_list = fits.open(outfilename)

print hdu_list.info()

image_data = hdu_list[0].data
    #if filename.endswith('NA.fits') == False:
    #bmaj = hdu_list[0].header['BMAJ']/hdu_list[0].header['CDELT2']
    #bmin = hdu_list[0].header['BMIN']/hdu_list[0].header['CDELT2']
#bpa = hdu_list[0].header['BPA']
hdu_list.close()

image_data = image_data*1000000
plt.figure(1)
plt.imshow(image_data, cmap='gray', origin='lower', interpolation='nearest')
plt.colorbar()

print('Min:', np.min(image_data))
print('Max:', np.max(image_data))
print('Mean:', np.mean(image_data))
print('Stdev:', np.std(image_data))

plt.figure(2)
NBINS = 100
histogram = plt.hist(image_data.flat, NBINS)
plt.yscale('log')
plt.show()

f = plt.figure(3)
gs = gridspec.GridSpec(2, 1,width_ratios=[1],height_ratios=[4,1])
gs.update(wspace=0.1, hspace=0.1)
hist, bin_edges = np.histogram(image_data.flat, bins=NBINS, density=False)
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
p0 = [1000., 1, 7]

coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

# Get the fitted curve
hist_fit = gauss(bin_centres, *coeff)
ax1 = f.add_subplot(gs[0])
histogram = plt.hist(image_data.flat, NBINS)
#plt.plot(bin_centres, hist, label='Test data')
list1, = plt.plot(bin_centres, hist_fit, label='$\sigma=$'+str(round(abs(coeff[2]),2))+',$\mu=$'+str(round(coeff[1],3)),color='y', linewidth=2.0)
ax1.set_title('Non-detection')

# Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
print 'Fitted mean = ', coeff[1]
print 'Fitted standard deviation = ', coeff[2]
for i in range(1,8):
    if i & 1:
        plt.axvline(i*coeff[2], color='k', linestyle='--')
        plt.text(i*abs(coeff[2]),np.amax(hist)*5,str(i)+'$\sigma$',rotation=90)
        plt.axvline(-i*coeff[2], color='k', linestyle='--')
#plt.text(-i*coeff[2],np.amax(hist)*5,str(-i)+'$\sigma$',rotation=90)
    if i == 6:
        plt.axvline(i*coeff[2], color='r', linestyle='-.')
        plt.text(i*abs(coeff[2]),np.amax(hist)*5,str(i)+'$\sigma$',rotation=90)
        plt.axvline(-i*coeff[2], color='r', linestyle='-.')
plt.ylim(0.1, np.amax(hist)*10)
plt.yscale('log')
plt.ylabel('Number of Pixels')
blue_patch = mpatches.Patch(color='blue', label='Data')
plt.legend(handles=[blue_patch,list1],loc=2,prop={'size':10})

ax2 = f.add_subplot(gs[1],sharex=ax1)
list2, = plt.plot(bin_centres,abs(hist-hist_fit)/hist,label='Residuals')
plt.xlabel('Peak Flux Density $\mu$Jy/beam')
plt.ylabel('(hist - model)/hist')
plt.legend(handles=[list2],loc=2,prop={'size':10})
plt.setp(ax1.get_xticklabels(), visible=True)
plt.setp([ax1, ax2])
gs.tight_layout(f, rect=[0, 0, 1, 0.97])
plt.show()

from astropy import wcs
    
from wcsaxes import WCSAxes
mywcs = wcs.WCS(hdu_list[0].header)
fig = plt.figure(6)
ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=mywcs)
lon = ax.coords['ra']
lat = ax.coords['dec']
lon.set_major_formatter('hh:mm:ss.ss')
lat.set_major_formatter('dd:mm:ss.ss')
lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
lat.set_axislabel('Declination (J2000)', minpad=1)
maxpixcoord = np.argwhere(image_data==np.max(image_data))
ax.set_xlim(maxpixcoord[0,1]-pixsiz,maxpixcoord[0,1]+pixsiz)
ax.set_ylim(maxpixcoord[0,0]-pixsiz,maxpixcoord[0,0]+pixsiz)
fig.add_axes(ax)
levs = [-1*np.std(image_data)]
levs = np.append(levs,np.linspace(np.std(image_data),np.max(image_data),8))
cont = ax.contour(image_data, levels=levs, cmap='gray_r', alpha=0.5)
im = ax.imshow(image_data, origin='lower',cmap="magma")
divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="5%", pad=0.00,
                          axes_class=matplotlib.axes.Axes)
cb = plt.colorbar(orientation="horizontal",mappable=im, cax=cax)
cb.ax.xaxis.set_ticks_position('top')
#cb = plt.colorbar(mappable=im, cax=cax,  ticks=np.linspace(0,1,2).astype(int))
cb.add_lines(cont)
cax.set_xlabel("Flux Density ($\mu$Jy)")
'''
oldlabels = cb.ax.get_yticklabels()
oldlabels[-1] ='900+'
newlabels = oldlabels
print(newlabels)
cb.ax.set_yticklabels(newlabels)
'''
#cb.ax.xaxis.set_label_position('bottom')

from matplotlib.patches import Ellipse
#if filename.endswith('NA.fits') == False:
#   circ=Ellipse((maxpixcoord[0,1]-(pixsiz-15),maxpixcoord[0,0]-(pixsiz-15)), width=bmin, height=bmaj, angle=bpa, fill=False, color='k', hatch='xxxxx')
#else:
circ=Ellipse((maxpixcoord[0,1]-(pixsiz-18),maxpixcoord[0,0]-(pixsiz-18)), width=14.2, height=15.6, angle=0, fill=False, color='w', hatch='xxx')
ax.add_patch(circ)
ax.set_zorder(20)
plt.show()
#plt.savefig(filename[:-4] + 'pdf')



'''
F = aplpy.FITSFigure(hdu_list[0], figure=plt.figure(7))
F.show_grayscale()
F.ticks.set_xspacing(0.000999999)  # degrees
F.ticks.set_yspacing(0.0008)
F.tick_labels.set_font(size='small')
F.axis_labels.set_xtext('Right Ascension (J2000)')
F.axis_labels.set_ytext('Declination (J2000)')
F.save('J123.pdf')
'''
