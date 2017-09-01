import os

for file in os.listdir('./'):
    if file.endswith('.fits'):
        importfits(fitsimage=file, imagename=file[:-5]+'_casa.image')
        exportfits(imagename=file[:-5]+'_casa.image',fitsimage=file[:-5]+'_casa.fits',dropdeg=True,dropstokes=True)
        os.system('rm -r '+file[:-5]+'_casa.image')
	os.system('rm casa-*.log ipython-*.log')
