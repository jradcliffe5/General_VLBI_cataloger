import os,sys
### Inputs ###
detection_threshold = 6
shorthand = True
postfix='test'
def write_catalog_pybdsf(input_image,detection_threshold,shorthand):
        if shorthand == True:
            name = input_image.split('_MSSC')[0]
        else:
            name = input_image
        img = bdsf.process_image(input_image, adaptive_rms_box=True,thresh_pix=detection_threshold,\
        thresh='hard')
        # Write the source list catalog. File is named automatically.
        img.write_catalog(format='csv', catalog_type='srl',clobber=True)
        # Write the residual image. File is name automatically.
        img.export_image(outfile=name+'gaus.residual.fits',img_type='gaus_resid',clobber=True)
        # Write the model image. File name is specified.
        img.export_image(outfile=name+'gaus.model.fits',img_type='gaus_model',clobber=True)

#write_catalog_pybdsf('HDFA0002_MSSC_FG_NA_IM.fits',detection_threshold)

def combine_pybdsf(shorthand,postfix):
    os.system('rm catalogue_%s.csv' % postfix)
    if os.path.isfile('catalogue_%s.csv' % postfix) == False:
        s = 'Name, Source_id, Isl_id, RA, E_RA, DEC, E_DEC, Total_flux, E_Total_flux, Peak_flux, E_Peak_flux, RA_max, E_RA_max, DEC_max, E_DEC_max, Maj, E_Maj, Min, E_Min, PA, E_PA, Maj_img_plane, E_Maj_img_plane, Min_img_plane, E_Min_img_plane, PA_img_plane, E_PA_img_plane, DC_Maj, E_DC_Maj, DC_Min, E_DC_Min, DC_PA, E_DC_PA, DC_Maj_img_plane, E_DC_Maj_img_plane, DC_Min_img_plane, E_DC_Min_img_plane, DC_PA_img_plane, E_DC_PA_img_plane, Isl_Total_flux, E_Isl_Total_flux, Isl_rms, Isl_mean, Resid_Isl_rms, Resid_Isl_mean, S_Code\n'.format(postfix)
        print s
        os.system('touch catalogue_%s.csv' % postfix)
        text_file = open('catalogue_%s.csv' % postfix,'a')
        text_file.write(s)
    for file in os.listdir('./'):
        if file.endswith('.srl'):
            lines = open('%s' % file).readlines()
            print lines
            print len(lines)
            if shorthand == True:
                names = file[:8]
            else:
                names = file
            if len(lines) > 6:
                #detections = detections + [file]
                print ''.join(lines[6:])
                text_file.write(names+','+''.join(lines[6:]))
            #os.system('rm %s.fitout' % file)

for i in os.listdir('./'):
    if i.endswith('IM.fits'):
        write_catalog_pybdsf(i,detection_threshold,shorthand)
combine_pybdsf(shorthand=True,postfix=postfix)
