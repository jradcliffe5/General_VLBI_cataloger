import os

for file in os.listdir('./'):
    if (os.path.isdir(file) == True) and ((file.endswith('.image')) or (file.endswith('.pbcor')) or (file.endswith('.pbcor.image.tt0'))):
        exportfits(imagename=file,fitsimage=file.split('.fits')[0]+'_casa.fits',\
        dropdeg=True,dropstokes=True)
    if file.endswith('.fits'):
        importfits(fitsimage=file, imagename=file[:-5]+'_casa.image')
        exportfits(imagename=file[:-5]+'_casa.image',fitsimage=file[:-5]+'_casa.fits',\
        dropdeg=True,dropstokes=True)
        os.system('rm -r '+file[:-5]+'_casa.image')
os.system('rm casa-*.log ipython-*.log')
