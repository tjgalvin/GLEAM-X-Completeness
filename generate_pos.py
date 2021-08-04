#!/usr/bin/env python

# Generates a list of random RA and Dec positions in specified region of sky

import numpy as np
from optparse import OptionParser
from astropy.coordinates import SkyCoord

### gcd() was written by Paul Hancock
# The following functions are explained at http://www.movable-type.co.uk/scripts/latlong.html
# phi ~ lat ~ Dec
# lambda ~ lon ~ RA
def gcd(ra1, dec1, ra2, dec2):
  """
  Great circle distance as calculated by the haversine formula
  ra/dec in degrees
  returns:
  sep in degrees
  """
  dlon = ra2 - ra1
  dlat = dec2 - dec1
  a = np.sin(np.radians(dlat) / 2) ** 2
  a += np.cos(np.radians(dec1)) * np.cos(np.radians(dec2)) * np.sin(np.radians(dlon) / 2) ** 2
  sep = np.degrees(2 * np.arcsin(min(1, np.sqrt(a))))
  return sep

# Read input parameters
usage="Usage: %prog [options] <output_file>"
parser = OptionParser(usage=usage)
parser.add_option('--nsrc',type="int", dest="nsrc", default=1000, help="Number of random source positions to simulate [default=%default]")
parser.add_option('--region',type="string", dest="region", default="0,360,-90,90", help="Region of sky. Enter ra_min, ra_max, dec_min, "\
"dec_max, in deg. For example, for 60 < RA < 300 enter: 60,300,-90,90. For RA < 60 or RA > 300, and -40 < Dec < -10, enter: "\
"300,60,-40,-10. [default=%default]")
parser.add_option('--sep_min',type="float", dest="sep_min", default=0.0, help="Minimum separation between simulated sources, in arcmin "\
"[default=%default]")
(options, args) = parser.parse_args()
output_file=args[0]
nsrc=options.nsrc
region=options.region
sep_min=options.sep_min

# Read region parameter
ra_min=float(region.split(',')[0])
ra_max=float(region.split(',')[1])
dec_min=float(region.split(',')[2])
dec_max=float(region.split(',')[3])

# Check RA and Dec ranges
if dec_min < -90.0 or dec_min > 90.0 or dec_max < -90.0 or dec_max > 90.0 or dec_min >= dec_max:
  print("Error: Dec range entered incorrectly. Aborting.")
  exit()
if ra_min < 0.0 or ra_min > 360.0 or ra_max < 0.0 or ra_max > 360.0 or ra_min == ra_max:
  print("Error: RA range entered incorrectly. Aborting.")
  exit()
if ra_max < ra_min:
  ra_max=ra_max+360.0

# Convert dec_min and dec_max from deg to rad
dec_min=np.radians(dec_min)
dec_max=np.radians(dec_max)

# Convert sep_min from arcmin to deg
sep_min=sep_min/60.0

i=0; n=0
ra_store=[]; dec_store=[]
step=1.0
f=step
while i < nsrc:
  n=n+1
# Draw random position inside specified region
  ra=np.random.uniform(low=ra_min, high=ra_max)
  if ra >= 360.0:
    ra=ra-360.0
  r=np.random.uniform(low=0.0, high=1.0)
  dec=np.arcsin((np.sin(dec_max)-np.sin(dec_min))*r+np.sin(dec_min))
  dec=np.degrees(dec) # Convert from rad to deg
# Reject position if it lies within sep_min deg from any other drawn position
  accept=True
  if sep_min > 0:
    for j in range(0,i):
      if np.absolute(dec-dec_store[j]) < sep_min:
        sep=gcd(ra, dec, ra_store[j], dec_store[j])
        if sep < sep_min:
          accept=False
          break
  if accept:
    ra_store.append(ra)
    dec_store.append(dec)
    i=i+1
  if i/nsrc >= f/100:
    print('%.0f' % f,'% complete')
    f=f+step
  if n>10*nsrc:
    print('Too many sources are being rejected as a result of the minimum separation constraint.', \
    '\nReduce sep_min parameter. Aborting.')
    exit()

# Print simulated source positions to file
f=open(output_file,'w')
print('# ra dec', file=f)
for i in range(0,nsrc):
  print('%.6f' % ra_store[i], '%.6f' % dec_store[i], file=f)

  

