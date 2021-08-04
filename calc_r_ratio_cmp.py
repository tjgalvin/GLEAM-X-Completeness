#!/usr/bin/env python

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import math as m
from optparse import OptionParser

# Read input parameters
usage="Usage: %prog [options] <srcname> <map> <map_psf> <output_file>\n"
parser = OptionParser(usage=usage)
parser.add_option('--z',type="float", dest="z", default=-26.7, help="Dec at zenith, in deg [default=%default]")
parser.add_option('--flux',type="float", dest="flux", default=1.0, help="Flux density of simulated sources, in Jy [default=%default]")
(options, args) = parser.parse_args()
srcname=args[0]; map=args[1]; map_psf=args[2]; output_file=args[3]
z=options.z
flux=options.flux

# Define output file
f=open(output_file,'w')

print('# RA Dec S bmaj bmin bpa R', file=f)

# -------------------------------------------

# Read file with list of source positions in deg
ra_list=[]
dec_list=[]
for line in open(srcname):
  columns = line.split()
  if columns[0] != '#':
    ra_list.append(float(columns[0]))
    dec_list.append(float(columns[1]))

# Read bmaj and bmin from fits header of surface brightness map
img = fits.open(map)
bmaj=img[0].header['BMAJ']
bmin=img[0].header['BMIN']

# Convert bmaj and bmin from deg to arcsec, and print to screen
bmaj=bmaj*3600.0
bmin=bmin*3600.0
print('bmaj / arcsec =',bmaj)
print('bmin / arcsec =',bmin)

# Read PSF map
img_psf = fits.open(map_psf)
mappix = (np.squeeze(img_psf[0].data))
img_psf.close()
header_psf = fits.getheader(map_psf)
w=WCS(header_psf)

for (ra,dec) in np.nditer((ra_list,dec_list)):  
# Convert (RA, Dec) position from degrees to pixel
  wpos = np.array([ra, dec, 0])
  wpos.shape = (1, 3)
  px=w.wcs_world2pix(wpos, 0)[0,0]
  py=w.wcs_world2pix(wpos, 0)[0,1]
  x=int(round(px))
  y=int(round(py))
  a=mappix[0,y,x]
  b=mappix[1,y,x]
  c=mappix[2,y,x]
  proj=1/m.cos(m.radians(z-dec))
  if a == a and b == b and c == c:
# Convert a and b from deg to arcsec
    a=a*3600.0
    b=b*3600.0
  else:
    a=bmaj*proj
    b=bmin
    c=0.0
# Calculate R ratio for each source - R = int_flux/peak_flux  
  r=a*b/(bmaj*bmin*proj)
  if r < 1.0:
    r=1.0
  print(ra,dec,flux,a,b,c,r, file=f)

    
