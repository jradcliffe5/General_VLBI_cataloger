# General_VLBI_cataloger

1. Use BANE to make rms map of the field. This is part of the AEGEAN package and can be installed easily using pip install AegeanTools  and then can be run using BANE <image.fits>
2. Use blobcat using the following inputs python blobcat.py --ppe=0.01 --pasbe=0.2 --dSNR=5 --fSNR=3 —rmsmap=<image_rms.fits> --edgemin=%d --ds9 --write <image.fits>'
3. Use the pyregion package to take the ds9 regions and output a pdf of the fits file to show the locations of the catalogue entries.

It is probably more useful to use the scripts I already have made to instantly make a catalogue and the region plots + blobcat inputs. You can download this using git clone https://github.com/jradcliffe5/General_VLBI_cataloger.git

To use the scripts you need to do the following steps:

1. Copy the scripts into the current directory with the files.

2. Use casa_convert_fits.py to remove degenerate axes from the images in the current directory (I would write this into a python only but casa works well). You run this with casa -c casa_convert_fits.py. It will convert all files with suffixes .pbcor .image .fits .image .tt0.

3. Edit catalog_inputs.txt to edit inputs

4. Then run the cataloger.py using ParselTongue cataloger.py (or comment out the import AIPS lines (8-11) etc. and use python cataloger.py).

5. If cataloger.py has ds9 = True, run python ds9_blob_plotter.py to create pdf images of the fields.
