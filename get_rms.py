from astropy.io import fits
import numpy as np
import os, sys
import pandas as pd

ra = []
dec = []
rms = []
name = []
for i in os.listdir('./'):
	if i.endswith('_casa.fits'):
		print('Getting rms for file %s'%i)
		hdu = fits.open(i)
		name.append('%s'%i)
		ra.append(hdu['PRIMARY'].header['CRVAL1'])
		dec.append(hdu['PRIMARY'].header['CRVAL2'])
		rms.append(np.std(hdu['PRIMARY'].data.squeeze()[100:1000,100:1000]))
		hdu.close()

pd.DataFrame({'src':name,'ra':ra,'dec':dec,'rms':rms}).to_csv('pbcor_rms.csv')
