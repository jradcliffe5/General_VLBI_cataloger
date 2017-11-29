import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib import rcParams
from astropy import wcs
from astropy.io import fits
import os.path
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import csv
import numpy as np
import pandas as pd


#### Matplotlib specifications ###
# All plots are 9x9 and use serif fonts
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
######################################

beaminfo = []
matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'

### Inputs ###
### Only take subsection of data (in pixels) to increase comp. speed
subimsize = 400
units = 'uJy'
nlevs = 5
### Set image size in pixels
pixsiz = 90
### Set edge no of pixels to not consider for rms calc
edge = 15
##############


### Scaler for the images ###
if units == 'uJy':
    flux_scaler = 1E6
elif units == 'mJy':
    flux_scaler = 1E3
elif units == 'Jy':
    flux_scaler = 1
else:
    print 'No unit specified [Jy/mJy/uJy], assuming Jy'
    flux_scaler = 1
#############################
### Looks in the current directory for fits files ending in .casa.pyfits

for file in os.listdir('./'):
    if file.endswith('casa.fits'):
        ## VLBI EG078B specifics
        '''
        pixsiz = 40
        edge = 5
        if file.endswith('NA_IM_casa.fits'):
            pixsiz = 60
            edge = 10
        if file.endswith('PBCOR_NA_IM_large_casa.fits'):
            pixsiz = 60
            edge = 10
        '''
        print 'Plotting %s' % file
        ### Pull fits file and info
        hdu_list = fits.open(file)
        mywcs = wcs.WCS(hdu_list[0].header)
        ### Get info on image size in pixels
        naxis = hdu_list[0].header['NAXIS1']
        ### Get data and set size of image (defined by subimsize)
        image_data = hdu_list[0].data[int((naxis-subimsize)/2.):int((naxis+subimsize)/2.),\
        int((naxis-subimsize)/2.):int((naxis+subimsize)/2.)]*flux_scaler
        ### Initialise figure and wcs axes
        ax = plt.subplot(projection=mywcs)
        ### Set coordinate formats
        lon = ax.coords['ra']
        lat = ax.coords['dec']
        lon.set_major_formatter('hh:mm:ss.sss')
        ### Some VLBI specific shit
        if file.endswith('NA_PBCOR_IM_casa.fits'):
            lon.set_major_formatter('hh:mm:ss.ss')
        lat.set_major_formatter('dd:mm:ss.ss')
        ### Set labels for axes
        lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
        lat.set_axislabel('Declination (J2000)', minpad=1)
        lon.set_ticks(number=3)
        ### Find maximum brightness
        maxpixcoord = np.argwhere(image_data==np.max(image_data))
        ax.set_xlim(maxpixcoord[0,1]-pixsiz,maxpixcoord[0,1]+pixsiz)
        ax.set_ylim(maxpixcoord[0,0]-pixsiz,maxpixcoord[0,0]+pixsiz)
        ### Add axes to figure
        ### Changes colorscale depending on dynamic range
        if np.max(image_data) > 11:
            levs = [-1*np.std(image_data),np.std(image_data)]
            levs = np.append(levs,np.around(np.linspace(np.std(image_data),np.max(image_data),nlevs),decimals=0)[1:-1])
        else:
            levs = [-1*np.std(image_data),np.std(image_data)]
            levs = np.append(levs,np.around(np.linspace(np.std(image_data),np.max(image_data),nlevs),decimals=0)[1:-1])
        ## Also contours same data
        cont = ax.contour(image_data, levels=levs, cmap='gray_r', alpha=0.5)
        ## Show image
        im = ax.imshow(image_data, origin='lower',cmap="magma",interpolation="bicubic")
        ### Makes colorbar on top
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="5%", pad=0.00,
                                  axes_class=matplotlib.axes.Axes)
        ## Dynamic range checks for tick labels so no overcrowding
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
        ### Set tick position
        cb.ax.xaxis.set_ticks_position('top')
        #cb = plt.colorbar(mappable=im, cax=cax,  ticks=np.linspace(0,1,2).astype(int))
        ### Add color bar and label
        cb.add_lines(cont)
        cax.set_xlabel("Flux Density ($\mathrm{\mu Jy\,beam^{-1}}$)", labelpad=-80)

        ### Pull beam information from fits file
        from matplotlib.patches import Ellipse
        ### Put bmaj, bmin in terms of pixel size!
        bmaj = hdu_list[0].header['BMAJ']/hdu_list[0].header['CDELT2']
        bmin = hdu_list[0].header['BMIN']/hdu_list[0].header['CDELT2']
        bpa = hdu_list[0].header['BPA']

        ### VLBI specific stuff (comment out if not working)
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
        ### Add beam information
        circ=Ellipse((maxpixcoord[0,1]-(pixsiz-edge),maxpixcoord[0,0]-(pixsiz-edge)), width=bmin, height=bmaj, angle=bpa, fill=False, color='w', hatch='xxxxx')
        ax.add_patch(circ)
        ax.set_zorder(20)
        ### Save figures
        plt.savefig(file[:-5]+'_plot.pdf', bbox_inches='tight', clobber=True)
        plt.close()
