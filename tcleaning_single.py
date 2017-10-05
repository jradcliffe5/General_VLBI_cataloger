##########################
######Inputs##############
outlierfile = 'catalogue_Taper_CASA_outliers.txt'
measurement_sets = ['JVLA1_peeled_uvsub_mtmfsub_eMERGE_DR1_averaged_LLRR_uvfix.ms', 'JVLA2_peeled_uvsub_mtmfsub_eMERGE_DR1_averaged_LLRR_uvfix.ms', 'JVLA3_peeled_uvsub_mtmfsub_eMERGE_DR1_averaged_LLRR_uvfix.ms', 'JVLA4_peeled_uvsub_mtmfsub_eMERGE_DR1_averaged_LLRR_uvfix.ms', 'JVLA5_peeled_uvsub_mtmfsub_eMERGE_DR1_averaged_LLRR_uvfix.ms', 'JVLA6_peeled_uvsub_mtmfsub_eMERGE_DR1_averaged_LLRR_uvfix.ms', 'JVLA7_peeled_uvsub_mtmfsub_eMERGE_DR1_averaged_LLRR_uvfix.ms', 'JVLA8_peeled_uvsub_mtmfsub_eMERGE_DR1_averaged_LLRR_uvfix.ms']
datacolumn = 'data'
imsize = 1024
cell = '0.2arcsec'
niter = 1000
gridder = 'standard'
deconvolver = 'mtmfs'
weighting = 'briggs'
robust = 2
parallel = True


lines = [line.rstrip('\n') for line in open(outlierfile)]

imagenames = []
phasecenters = []
for line in lines:
    if line.startswith('imagename='):
        imagenames = imagenames + [line.strip('imagename=')]
    if line.startswith('phasecenter='):
        phasecenters = phasecenters + [line.strip('phasecenter=')]


for i in range(len(phasecenters)):
    tclean(vis=measurement_sets,selectdata=True,field="",spw="",timerange="",uvrange="",antenna="",scan="",observation="",intent="",datacolumn=datacolumn,imagename=imagenames[i],imsize=imsize,cell=cell,phasecenter=phasecenters[i],\
    stokes="I",projection="SIN",startmodel="",specmode="mfs",reffreq="",nchan=-1,start="",width="",outframe="LSRK",veltype="radio",restfreq=[],interpolation="linear",gridder=gridder,facets=1,chanchunks=1,wprojplanes=1,vptable="",\
    aterm=True,psterm=False,wbawp=True,conjbeams=True,cfcache="",computepastep=360.0,rotatepastep=360.0,pblimit=0.2,normtype="flatnoise",deconvolver=deconvolver,scales=[],nterms=2,smallscalebias=0.6,restoration=True,restoringbeam=[],\
    pbcor=False,outlierfile="",weighting=weighting,robust=robust,npixels=0,uvtaper=[],niter=niter,gain=0.1,threshold=8e-06,cycleniter=-1,cyclefactor=1.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,usemask="user",mask="",pbmask=0.0,maskthreshold="",maskresolution="",nmask=0,sidelobethreshold=3.0,noisethreshold=5.0,lownoisethreshold=1.5,smoothfactor=1.0,minbeamfrac=0.3,cutthreshold=0.01,restart=True,savemodel="none",calcres=True,calcpsf=True,parallel=parallel)
