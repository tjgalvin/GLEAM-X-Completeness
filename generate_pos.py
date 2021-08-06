#!/usr/bin/env python

# Generates a list of random RA and Dec positions in specified region of sky

import numpy as np
from argparse import ArgumentParser
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u

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
    a += (
        np.cos(np.radians(dec1))
        * np.cos(np.radians(dec2))
        * np.sin(np.radians(dlon) / 2) ** 2
    )
    sep = np.degrees(2 * np.arcsin(min(1, np.sqrt(a))))
    return sep


# Read input parameters
description = "Generates a set os source positions with a minimum separation"
parser = ArgumentParser(description=description)
parser.add_argument(
    "outfile", type=str, help="Output path to write the source positions to"
)
parser.add_argument(
    "--nsrc",
    type=int,
    default=1000,
    help="Number of random source positions to simulate [default=%default]",
)
parser.add_argument(
    "--region",
    type=str,
    default="0,360,-90,90",
    help="Region of sky. Enter ra_min, ra_max, dec_min, "
    "dec_max, in deg. For example, for 60 < RA < 300 enter: 60,300,-90,90. For RA < 60 or RA > 300, and -40 < Dec < -10, enter: "
    "300,60,-40,-10. [default=%default]",
)
parser.add_argument(
    "--sep-min",
    type=float,
    default=0.0,
    help="Minimum separation between simulated sources, in arcmin "
    "[default=%default]",
)
parser.add_argument(
    "--max-attempts",
    type=int,
    default=50,
    help="Maximum number of passes to make over generating valid source positions before failing. ",
)

args = parser.parse_args()
output_file = args.outfile
nsrc = args.nsrc
region = args.region
sep_min = args.sep_min
max_attempts = args.max_attempts

# Read region parameter
ra_min = float(region.split(",")[0])
ra_max = float(region.split(",")[1])
dec_min = float(region.split(",")[2])
dec_max = float(region.split(",")[3])

# Check RA and Dec ranges
if (
    dec_min < -90.0
    or dec_min > 90.0
    or dec_max < -90.0
    or dec_max > 90.0
    or dec_min >= dec_max
):
    print("Error: Dec range entered incorrectly. Aborting.")
    exit()
if ra_min < 0.0 or ra_min > 360.0 or ra_max < 0.0 or ra_max > 360.0 or ra_min == ra_max:
    print("Error: RA range entered incorrectly. Aborting.")
    exit()
if ra_max < ra_min:
    ra_max = ra_max + 360.0

# Convert dec_min and dec_max from deg to rad
dec_min = np.radians(dec_min)
dec_max = np.radians(dec_max)

# Convert sep_min from arcmin to deg
sep_min = sep_min / 60.0

ras = np.zeros(nsrc)
decs = np.zeros(nsrc)
mask = ras == 0
min_sep = args.sep_min * u.arcmin

c = 0
success = False
while (success is False) and (c < max_attempts):
    ras[mask] = np.random.uniform(low=ra_min, high=ra_max, size=np.sum(mask))
    ras = np.where(ras >= 360, ras - 360, ras)

    r = np.random.uniform(low=0.0, high=1.0, size=np.sum(mask))
    decs[mask] = np.rad2deg(
        np.arcsin((np.sin(dec_max) - np.sin(dec_min)) * r + np.sin(dec_min))
    )

    pos = SkyCoord(ras * u.deg, decs * u.deg)
    seps = match_coordinates_sky(pos, pos, nthneighbor=2)

    mask = seps[1] < min_sep
    print(f"Pass {c+1} : Accepted {np.sum(~mask)}")
    if np.all(~mask):
        success = True

    c += 1


if not success:
    print("Injection did not converge. Was --sep-min to large?")
    exit(1)

# Print simulated source positions to file
with open(output_file, "w") as f:
    print("# ra dec", file=f)
    for i in range(0, nsrc):
        print("%.6f" % ras[i], "%.6f" % decs[i], file=f)

