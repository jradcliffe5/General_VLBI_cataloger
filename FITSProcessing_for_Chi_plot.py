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

matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'

### Inputs ###
subimsize = 4000
units = 'uJy'
nlevs = 5
df = pd.read_csv('/Users/jackradcliffe/PhD/GOODSN_Catalogues/eg078_EVN_catalogues/eg078b/VLBI_Catalogue_fix_v14.csv')
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

beaminfo = []
matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
for file in os.listdir('./'):
    if file[:8].upper() in df.Catalog_name.tolist():
        name =  df.NAME_VLA_1[df.Catalog_name == file[:8].upper()].tolist()[0]
    if file.endswith('casa.fits'):
        pixsiz = 400
        edge = 50
        if file.endswith('PBCOR_NA_IM_casa.fits'):
            pixsiz = 100
            edge = 10
        if file.endswith('PBCOR_NA_IM_large_casa.fits'):
            pixsiz = 600
            edge = 10
        print 'Plotting for Chi plot %s' % file
        hdu_list = fits.open(file)
        mywcs = wcs.WCS(hdu_list[0].header)
        naxis = hdu_list[0].header['NAXIS1']
        image_data = hdu_list[0].data[int((naxis-subimsize)/2.):int((naxis+subimsize)/2.),\
        int((naxis-subimsize)/2.):int((naxis+subimsize)/2.)]*flux_scaler
        fig = plt.figure(6)
        ### Initialise figure and wcs axes
        ax = plt.subplot(projection=mywcs)
        lon = ax.coords['ra']
        lat = ax.coords['dec']
        lon.set_ticks_visible(False)
        lon.set_ticklabel_visible(False)
        lat.set_ticks_visible(False)
        lat.set_ticklabel_visible(False)
        lon.set_axislabel('', minpad=1.5)
        lat.set_axislabel('', minpad=1)
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
        namesize = 82
        text = ax.text(0.98, 0.9,r'\textbf{J%s}' % (name.split('+')[0]),color='w', horizontalalignment='right', verticalalignment='center',transform = ax.transAxes,size=namesize)
        text2 = ax.text(0.98, 0.76,r'\textbf{+%s}' % (name.split('+')[1]),color='w', horizontalalignment='right', verticalalignment='center',transform = ax.transAxes,size=namesize)


        from matplotlib.patches import Ellipse
        bmaj = hdu_list[0].header['BMAJ']/hdu_list[0].header['CDELT2']
        bmin = hdu_list[0].header['BMIN']/hdu_list[0].header['CDELT2']
        bpa = hdu_list[0].header['BPA']
        beam_text = ax.text(0.98, 0.05,r'\textbf{%.1f$\times$%.1f mas}' % (hdu_list[0].header['BMAJ']*1000*3600,hdu_list[0].header['BMIN']*1000*3600),color='w', horizontalalignment='right', verticalalignment='center',transform = ax.transAxes,size=60)
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
        plt.savefig(file[:-5]+'_Chi_plot.pdf', bbox_inches='tight', clobber=True)
        plt.close(fig)
