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
import pandas as pd
df = pd.read_csv('/Users/jackradcliffe/PhD/GOODSN_Catalogues/eg078_EVN_catalogues/VLBI_Catalogue_v10.csv')
print df.columns
beaminfo = []
matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
for file in os.listdir('./'):
    if file[:8].lower() in df.Catalog_names.tolist():
        name =  df.NAME_VLA_1[df.Catalog_names == file[:8].lower()].tolist()[0]
        print name
    if file.endswith('casa.fits'):
        pixsiz = 40
        edge = 5
        if file.endswith('PBCOR_NA_IM_casa.fits'):
            pixsiz = 60
            edge = 10
        if file.endswith('PBCOR_NA_IM_large_casa.fits'):
            pixsiz = 60
            edge = 10
        hdu_list = fits.open(file)
        mywcs = wcs.WCS(hdu_list[0].header)
        image_data = hdu_list[0].data*1E6
        fig = plt.figure(6)
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=mywcs)
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
            levs = np.append(levs,np.around(np.linspace(np.std(image_data),np.max(image_data),7),decimals=0)[1:-1])
        else:
            levs = [-1*np.std(image_data),np.std(image_data)]
            levs = np.append(levs,np.around(np.linspace(np.std(image_data),np.max(image_data),7),decimals=0)[1:-1])
        print levs
        cont = ax.contour(image_data, levels=levs, cmap='gray_r', alpha=0.5)
        im = ax.imshow(image_data, origin='lower',cmap="magma",interpolation="bicubic")
        text = ax.text(0.5, 0.9,'J'+str(name),color='w', horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,size=30)
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("top", size="5%", pad=0.00,
        #                          axes_class=matplotlib.axes.Axes)
        #if np.max(image_data) > 11:
        #    tick = [np.min(image_data)]
        #    tick = np.append(tick,np.around(np.linspace(np.std(image_data),np.max(image_data),7),decimals=0)[1:])
            #tick = np.append(tick,np.max(image_data))
        #    tick= tick.astype(int)
        #    cb = plt.colorbar(orientation="horizontal",mappable=im, cax=cax, ticks=tick)
        #else:
        #    tick = [np.min(image_data)]
        #    tick = np.append(tick,np.around(np.linspace(np.std(image_data),np.max(image_data),7),decimals=0)[1:].astype(int))
            #tick = np.append(tick,np.max(image_data))
        #    tick = tick.astype(int)
        #    cb = plt.colorbar(orientation="horizontal",mappable=im, cax=cax, ticks=tick)
        #cb.ax.xaxis.set_ticks_position('top')
        #cb = plt.colorbar(mappable=im, cax=cax,  ticks=np.linspace(0,1,2).astype(int))
        #cb.add_lines(cont)
        #cax.set_xlabel("Flux Density ($\mathrm{\mu Jy\/bm^{-1}}$)", labelpad=-50)
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
        beam_text = ax.text(0.98, 0.05,r'%.1f$\times$%.1f mas' % (bmaj,bmin),color='w', horizontalalignment='right', verticalalignment='center',transform = ax.transAxes,size=24)
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
