import os,sys,re

def headless(inputfile):
    ''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
    INPUTFILE = open(inputfile, "r")
    control = {}
    # a few useful regular expressions
    newline = re.compile(r'\n')
    space = re.compile(r'\s')
    char = re.compile(r'\w')
    comment = re.compile(r'#.*')
    # parse the input file assuming '=' is used to separate names from values
    for line in INPUTFILE:
        if char.match(line):
            line = comment.sub(r'', line)
            line = line.replace("'", '')
            (param, value) = line.split('=')
            param = newline.sub(r'', param)
            param = param.strip()
            param = space.sub(r'', param)
            value = newline.sub(r'', value)
            value = value.replace(' ','').strip()
            valuelist = value.split(',')
            if len(valuelist) == 1:
                if valuelist[0] == '0' or valuelist[0]=='1' or valuelist[0]=='2':
                    control[param] = int(valuelist[0])
                else:
                    control[param] = str(valuelist[0])
            else:
                control[param] = ','.join(valuelist)
    return control

inputs = headless('catalog_inputs.txt')
detection_threshold = float(inputs['S_N_ratio'])
postfix = str(inputs['postfix'])
shorthand = inputs['shorthand']
ds9 = inputs['ds9']

def write_catalog_pybdsf(input_image,detection_threshold,shorthand):
        if shorthand == 'True':
            name = input_image.split('_MSSC')[0]
        else:
            name = input_image
        img = bdsf.process_image(input_image, adaptive_rms_box=True,thresh_pix=detection_threshold,\
        thresh='hard')
        # Write the source list catalog. File is named automatically.
        img.write_catalog(format='csv', catalog_type='srl',clobber=True)
        if ds9 == 'True':
            img.write_catalog(format='ds9', catalog_type='srl',clobber=True)
        # Write the residual image. File is name automatically.
        img.export_image(outfile=name+'_gaus.residual.fits',img_type='gaus_resid',clobber=True)
        # Write the model image. File name is specified.
        img.export_image(outfile=name+'_gaus.model.fits',img_type='gaus_model',clobber=True)

#write_catalog_pybdsf('HDFA0002_MSSC_FG_NA_IM.fits',detection_threshold)

def combine_pybdsf(shorthand,postfix):
    os.system('rm catalogue_pybdsf_%s.csv' % postfix)
    if os.path.isfile('catalogue_pybdsf_%s.csv' % postfix) == False:
        s = 'Name, Source_id, Isl_id, RA, E_RA, DEC, E_DEC, Total_flux, E_Total_flux, Peak_flux, E_Peak_flux, RA_max, E_RA_max, DEC_max, E_DEC_max, Maj, E_Maj, Min, E_Min, PA, E_PA, Maj_img_plane, E_Maj_img_plane, Min_img_plane, E_Min_img_plane, PA_img_plane, E_PA_img_plane, DC_Maj, E_DC_Maj, DC_Min, E_DC_Min, DC_PA, E_DC_PA, DC_Maj_img_plane, E_DC_Maj_img_plane, DC_Min_img_plane, E_DC_Min_img_plane, DC_PA_img_plane, E_DC_PA_img_plane, Isl_Total_flux, E_Isl_Total_flux, Isl_rms, Isl_mean, Resid_Isl_rms, Resid_Isl_mean, S_Code\n'.format(postfix)
        os.system('touch catalogue_pybdsf_%s.csv' % postfix)
        text_file = open('catalogue_pybdsf_%s.csv' % postfix,'a')
        text_file.write(s)
    for file in os.listdir('./'):
        if file.endswith('.srl'):
            lines = open('%s' % file).readlines()
            if shorthand == 'True':
                names = file[:8]+','
            else:
                names = file+','
            if len(lines) > 6:
                #detections = detections + [file]
                #print names+names.join(lines[6:])
                text_file.write(names+names.join(lines[6:]))
            os.system('rm %s' % file)

for i in os.listdir('./'):
    if i.endswith('IM.fits'):
        write_catalog_pybdsf(i,detection_threshold,shorthand)
combine_pybdsf(shorthand=True,postfix=postfix)
