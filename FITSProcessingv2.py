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
from wcsaxes import WCSAxes
from astropy import wcs
import csv
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.odr import Model, Data, ODR
from scipy.stats import linregress
import numpy as np
import matplotlib.cm as cm
from scipy import constants
import pandas as pd
from matplotlib import rc
from matplotlib import rcParams

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams.update({'font.size': 22})
plt.ioff()

fig_size = plt.rcParams["figure.figsize"]
# Prints: [8.0, 6.0]
print "Current size:", fig_size

# Set figure width to 9 and height to 9
fig_size[1] = 9
fig_size[0] = 9
plt.rcParams["figure.figsize"] = fig_size

beaminfo = []
matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'

### Inputs ###
subimsize = 1000
units = 'uJy'
nlevs = 5
##############

if units == 'uJy':
    flux_scaler = 1E6
elif units == 'mJy':
    flux_scaler = 1E3
elif units == 'Jy':
    flux_scaler = 1
else:
    print 'No unit specified [Jy/mJy/uJy], assuming Jy'
    flux_scaler = 1

for file in os.listdir('./'):
    if file.endswith('casa.fits'):
        pixsiz = 40
        edge = 5
        if file.endswith('NA_IM_casa.fits'):
            pixsiz = 60
            edge = 10
        if file.endswith('PBCOR_NA_IM_large_casa.fits'):
            pixsiz = 60
            edge = 10
        print 'Plotting %s' % file
        hdu_list = fits.open(file)
        mywcs = wcs.WCS(hdu_list[0].header)
        naxis = hdu_list[0].header['NAXIS1']
        image_data = hdu_list[0].data[int((naxis-subimsize)/2.):int((naxis+subimsize)/2.),\
        int((naxis-subimsize)/2.):int((naxis+subimsize)/2.)]*flux_scaler
        fig = plt.figure(6)
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=mywcs)
        lon = ax.coords['ra']
        lat = ax.coords['dec']
        lon.set_major_formatter('hh:mm:ss.sss')
        if file.endswith('NA_PBCOR_IM_casa.fits'):
            lon.set_major_formatter('hh:mm:ss.ss')
        lat.set_major_formatter('dd:mm:ss.ss')
        lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
        lat.set_axislabel('Declination (J2000)', minpad=1)
        lon.set_ticks(number=3)
        maxpixcoord = np.argwhere(image_data==np.max(image_data))
        ax.set_xlim(maxpixcoord[0,1]-pixsiz,maxpixcoord[0,1]+pixsiz)
        ax.set_ylim(maxpixcoord[0,0]-pixsiz,maxpixcoord[0,0]+pixsiz)
        fig.add_axes(ax)
        if np.max(image_data) > 11:
            levs = [-1*np.std(image_data),np.std(image_data)]
            levs = np.append(levs,np.around(np.linspace(np.std(image_data),np.max(image_data),nlevs),decimals=0)[1:-1])
        else:
            levs = [-1*np.std(image_data),np.std(image_data)]
            levs = np.append(levs,np.around(np.linspace(np.std(image_data),np.max(image_data),nlevs),decimals=0)[1:-1])
        cont = ax.contour(image_data, levels=levs, cmap='gray_r', alpha=0.5)
        im = ax.imshow(image_data, origin='lower',cmap="magma",interpolation="bicubic")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="5%", pad=0.00,
                                  axes_class=matplotlib.axes.Axes)
        if np.max(image_data) > 11:
            tick = [np.min(image_data)]
            tick = np.append(tick,np.around(np.linspace(np.std(image_data),np.max(image_data),nlevs),decimals=0)[1:])
            tick = np.append(tick,np.max(image_data))
            tick= tick.astype(int)
            cb = plt.colorbar(orientation="horizontal",mappable=im, cax=cax, ticks=tick)
        else:
            tick = [np.min(image_data)]
            tick = np.append(tick,np.around(np.linspace(np.std(image_data),np.max(image_data),nlevs),decimals=0)[1:].astype(int))
            tick = np.append(tick,np.max(image_data))
            tick = tick.astype(int)
            cb = plt.colorbar(orientation="horizontal",mappable=im, cax=cax, ticks=tick)
        cb.ax.xaxis.set_ticks_position('top')
        #cb = plt.colorbar(mappable=im, cax=cax,  ticks=np.linspace(0,1,2).astype(int))
        cb.add_lines(cont)
        cax.set_xlabel("Flux Density ($\mathrm{\mu Jy\,beam^{-1}}$)", labelpad=-80)
        '''
        oldlabels = cb.ax.get_yticklabels()
        oldlabels[-1] ='900+'
        newlabels = oldlabels
        print(newlabels)
        cb.ax.set_yticklabels(newlabels)
        '''
        #cb.ax.xaxis.set_label_position('bottom')

        from matplotlib.patches import Ellipse
        bmaj = hdu_list[0].header['BMAJ']/hdu_list[0].header['CDELT2']
        bmin = hdu_list[0].header['BMIN']/hdu_list[0].header['CDELT2']
        bpa = hdu_list[0].header['BPA']
        if file.endswith('NA_PBCOR_IM_casa.fits'):
            with open('beamsizes_taper.csv', 'a') as myfile:
                beaminfo = [file, hdu_list[0].header['BMAJ']*1000*3600, hdu_list[0].header['BMIN']*1000*3600, bpa]
                wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
                wr.writerow(beaminfo)
        else:
            with open('beamsizes_cal.csv', 'a') as myfile:
                beaminfo = [file, hdu_list[0].header['BMAJ']*1000*3600, hdu_list[0].header['BMIN']*1000*3600, bpa]
                wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
                wr.writerow(beaminfo)
        circ=Ellipse((maxpixcoord[0,1]-(pixsiz-edge),maxpixcoord[0,0]-(pixsiz-edge)), width=bmin, height=bmaj, angle=bpa, fill=False, color='w', hatch='xxxxx')
        ax.add_patch(circ)
        ax.set_zorder(20)
        plt.savefig(file[:-5]+'_plot.pdf', bbox_inches='tight', clobber=True)
        plt.close(fig)
