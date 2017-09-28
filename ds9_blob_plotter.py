import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import fits
import pyregion
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import numpy as np
import os

def ds9_plotter(file,reg_name,colormap,log_scaler):
    plt.clf()
    f_xray = fits.open(file)
    try:
        from astropy.wcs import WCS
        from astropy.visualization.wcsaxes import WCSAxes

        wcs = WCS(f_xray[0].header)
        fig = plt.figure(figsize=(90, 90))
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
        fig.add_axes(ax)
    except ImportError:
        ax = plt.subplot(111)

    image_data = f_xray[0].data*1E6
    rms = fits.open(file[:-5]+'_rms.fits')[0].data*1E6
    lon = ax.coords['ra']
    lat = ax.coords['dec']
    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm:ss')
    lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
    lat.set_axislabel('Declination (J2000)', minpad=1)
    rms = np.std(rms)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.00,axes_class=matplotlib.axes.Axes)
    im = ax.imshow(image_data, origin='lower',cmap=colormap, vmin=rms/2., vmax=10*rms)
    cb = plt.colorbar(orientation="horizontal",mappable=im, cax=cax, ticks=np.around(np.linspace(rms/2.,10*rms,4),decimals=-1),extend='max')
    cax.set_xlabel("Flux Density ($\mathrm{\mu Jy\/bm^{-1}}$)", labelpad=-60)
    cb.ax.xaxis.set_ticks_position('top')
    #cb = plt.colorbar(mappable=im, cax=cax,  ticks=np.linspace(0,1,2).astype(int))

    r = pyregion.open(reg_name).as_imagecoord(header=f_xray[0].header)

    from pyregion.mpl_helper import properties_func_default


    # Use custom function for patch attribute
    def fixed_color(shape, saved_attrs):
        attr_list, attr_dict = saved_attrs
        attr_dict["color"] = "green"
        kwargs = properties_func_default(shape, (attr_list, attr_dict))

        return kwargs

    r1 = pyregion.ShapeList(r)
    patch_list1, artist_list1 = r1.get_mpl_patches_texts(fixed_color)

    for p in patch_list1:
        ax.add_patch(p)
    for t in artist_list1:
        ax.add_artist(t)
    plt.savefig(fitsfile[:-5]+'.pdf',bbox_inches='tight',overwrite=True)
    f_xray.close()
    plt.close()
    del f_xray, image_data

# read in the image
for file in os.listdir('./'):
    if file.endswith('ds9.reg'):
        fitsfile = file[:-8]+'.fits'
        reg_name = file
        ds9_plotter(fitsfile,reg_name,colormap='magma',log_scaler=250)
