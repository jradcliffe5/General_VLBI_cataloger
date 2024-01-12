from astropy import wcs
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from matplotlib import *
from matplotlib import rc
from matplotlib import rcParams
import pandas as pd
from matplotlib import gridspec

rc('font', **{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=False)
rcParams['mathtext.default'] = 'regular'
figsize = plt.rcParams["figure.figsize"]
figsize[1]=9
figsize[0]=9
plt.rcParams["figure.figsize"]=figsize
matplotlib.rcParams.update({'font.size': 22})



def grid(x, y, z, resX=1000j, resY=1000j):
    "Convert 3 column data t matplotlib grid"
    xi, yi = np.mgrid[min(x):max(x):resX, min(y):max(y):resY]
    Z = griddata(points=(x,y), values=z, xi=(xi,yi),method='linear')
    #Z = griddata(x, y, z, xi, yi,interp='linear')
    #Z2 = ndimage.gaussian_filter(Z, sigma=1.0, order=1)
    #Z = interp2d(x, y, z, kind='cubic')
    #Z2 = Z(xi,yi)
    X, Y = xi, yi
    return X, Y, Z

def generate_central_wcs(crval, cdelt, crpix):
	# Create a new WCS object.  The number of axes must be set
	# from the start
	w = wcs.WCS(naxis=2)

	# Set up an "Airy's zenithal" projection
	# Vector properties may be set with Python lists, or Numpy arrays
	#CTYPE1  = projection
	#CRVAL1  = central position in degrees
	#CDELT1  = pixel demarcation
	#CRPIX1  = reference pixel
	#CUNIT1  = values of angle objects
	w.wcs.crpix = np.array(crpix).astype(int)
	w.wcs.cdelt = np.array(cdelt).astype(float)
	w.wcs.crval = np.array(crval).astype(float)
	w.wcs.ctype = ["RA---SIN", "DEC--SIN"]

	# Some pixel coordinates of interest.
	pixcrd = np.array([[-10, -10], [24, 38], [45, 98]], np.float_)

	# Convert pixel coordinates to world coordinates
	world = w.wcs_pix2world(pixcrd, 1)
	#print(world)

	# Convert the same coordinates back to pixel coordinates.
	pixcrd2 = w.wcs_world2pix(world, 1)
	#print(pixcrd2)

	# These should be the same as the original pixel coordinates, modulo
	# some floating-point error.
	assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6

	return w

## Create array of ra, dec, & intensity
df = pd.read_csv('pbcor_rms.csv')
RA= df['ra']
DEC = df['dec']
rms = df['rms']*1e6
#print(min(rms),max(rms))

w = generate_central_wcs(crval=[np.average(RA),np.average(DEC)],cdelt=[1./60.,1./60.],crpix=[0,0])

c = SkyCoord(RA,DEC,unit='deg',frame='icrs')

### Create axes & scatter plot
gs = gridspec.GridSpec(ncols=1,nrows=2,height_ratios=[1/20.,1],wspace=0,hspace=0)
fig = plt.figure()
ax = fig.add_subplot(gs[1], projection=w)
lon = ax.coords['ra']
lat = ax.coords['dec']
lon.set_major_formatter('hh:mm:ss')
lat.set_major_formatter('dd:mm:ss')
lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
lat.set_axislabel('Declination (J2000)', minpad=1)
X, Y, Z = grid(c.ra.degree, c.dec.degree, rms)
Z_gauss = np.nan_to_num(np.array(gaussian_filter(Z,1.5)))
#Z_gauss = gaussian_filter(Z,2)
v = [min(rms),250]
Z_gauss[Z_gauss==0] = 250
im = ax.pcolormesh(X,Y,Z_gauss,transform=ax.get_transform('icrs'),cmap='magma',vmin=v[0], vmax=v[1], alpha=1)
ax.scatter(c.ra.degree, c.dec.degree, transform=ax.get_transform('icrs'), marker='+',color='g',alpha=1)
xlim = ax.get_xlim()
ylim = ax.get_ylim()
CS = ax.contour(X,Y,gaussian_filter(Z,10),levels=np.linspace(v[0]*2,v[1]*0.9,5),colors='w',transform=ax.get_transform('world'))
ax.set_xlim(xlim[0]*0.9,xlim[1]*1.1)
ax.set_ylim(ylim[0]*0.9,ylim[1]*1.1)
cax = fig.add_subplot(gs[0])
cb = plt.colorbar(orientation="horizontal",mappable=im, cax=cax,format='%d',ticks=np.linspace(v[0]*2,v[1],6),extend='max')
cb.add_lines(CS)
ticks = np.append(np.array([v[0]]),np.linspace(v[0]*2,v[1]*0.9,5))
ticks = np.append(ticks,v[1])
print(ticks)
cb.set_ticks(ticks)
cax.xaxis.set_ticks_position('top')
cax.set_xlabel('rms sensitivity (uJy/beam)',labelpad=-75)
fig.savefig('rms_pbcor.png',bbox_inches='tight',dpi=fig.dpi,format="png")
