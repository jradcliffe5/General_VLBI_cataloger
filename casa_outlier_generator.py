import pandas as pd
import os

def outlier_writer(RA,Dec,names,postfix,imsize,cell,gridder,deconvolver,nterms):
    #options = ['#outlier','imagename', 'imsize','cell', 'phasecenter', 'specmode', 'nchan', 'nterms', 'gridder', 'deconvolver', 'wprojplanes']
    #gridder=standard, wproject, widefield, mosaic, awproject
    #deconvolver=hogbom, clark, clarkstokes,mtmfs, mem
    os.system('touch catalogue_%s_CASA_outliers.txt' % postfix)
    text_file = open('catalogue_%s_CASA_outliers.txt' % postfix,'a')
    for i in range(len(names)):
        text_file.write('#outlier %s\n' % names[i])
        text_file.write('imagename=%s\n' % names[i])
        text_file.write('imsize=[%s,%s]\n' % (imsize,imsize))
        text_file.write('cell=[%s,%s]\n' % (cell,cell))
        text_file.write('phasecenter=J2000 %sdeg %sdeg\n' % (RA[i],Dec[i]))
        text_file.write('gridder=%s\n' % gridder)
        text_file.write('deconvolver=%s\n' % deconvolver)
        text_file.write('nterms=%s\n' % nterms)
        text_file.write('#\n')

for file in os.listdir('./'):
    if file.startswith('catalogue') and file.endswith('.csv'):
        postfix = file.split('_')[1].split('.')[0]
        df = pd.read_csv(file,delimiter=',')
        if RA < 0:
            RA = 360+df['RA_c_%s' % postfix]
        else:
            RA = df['RA_c_%s' % postfix]
        Dec= df['Dec_c_%s' % postfix]
        names = df['Catalog_name']

outlier_writer(RA=RA,Dec=Dec,names=names,postfix=postfix,imsize=1024,cell='0.2arcsec',\
gridder='standard',deconvolver='mtmfs',nterms=2)
