#!/usr/bin/env python

# (1) Calculate overall completeness as a function of flux density
# (2) Generate completeness map at each flux density (optional)

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from optparse import OptionParser
import glob
import matplotlib.pyplot as plt
import math

# Define function to read catalogue with detected simulated sources
def read_catalogue(input_catalogue):
    ra = []
    dec = []
    table = fits.open(input_catalogue)[1].data
    ra = table.field("ra")
    dec = table.field("dec")
    return ra, dec


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


# Define function to count number of sources within a given radius from a reference position
# Input coordinates and search radius should be in deg
def count_sources(ra_ref, dec_ref, ra_list, dec_list, rad):
    n = 0
    for i in range(0, len(ra_list)):
        if np.absolute(dec_ref - dec_list[i]) <= rad:
            sep = gcd(ra_ref, dec_ref, ra_list[i], dec_list[i])
            if sep <= rad:
                n = n + 1
    return n


# Read input parameters
usage = (
    "Usage: %prog [options] <injected_sources> <detected_sources> \n"
    "injected_sources: file containing the positions of the injected sources \n"
    "detected_sources: directory containing catalogues of detected simulated sources"
)
parser = OptionParser(usage=usage)
parser.add_option(
    "--flux",
    type="string",
    dest="flux",
    help="Fluxes at which to calculate completeness. "
    "Enter flux_min, flux_max, flux_interval, in mJy (e.g. 5,10,1). This parameter must be given. No default.",
)
parser.add_option(
    "--template_map",
    type="string",
    dest="template_map",
    help="Template FITS map for "
    "generating completeness maps. "
    "This can be set to the PSF map. If not given, no completeness maps are generated. No default.",
)
parser.add_option(
    "--region",
    type="string",
    dest="region",
    default="0,360,-90,90",
    help="Sky coverage of "
    "completeness maps. This should normally correspond to the area of sky in which the simulated sources "
    "were injected. Enter ra_min, ra_max, dec_min, dec_max, in deg. For example, for 60 < RA < 300 enter: "
    "60,300,-90,90. For RA < 60 or RA > 300, and -40 < Dec < -10, enter: 300,60,-40,-10. [default=%default]",
)
parser.add_option(
    "--rad",
    type="float",
    dest="rad",
    default=6.0,
    help="Radius of circle centred on each "
    "pixel in which to calculate the completeness, in deg [default=%default]",
)
parser.add_option(
    "--output_root",
    type="string",
    dest="output_root",
    default="completeness",
    help="Output " "root filename [default=%default]",
)
(options, args) = parser.parse_args()
injected_sources = args[0]
detected_sources = args[1]
if not options.flux:
    parser.error("Flux parameter not given")

# Read region parameter
ra_min = float(options.region.split(",")[0])
ra_max = float(options.region.split(",")[1])
dec_min = float(options.region.split(",")[2])
dec_max = float(options.region.split(",")[3])

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

# Read flux parameter and create array with flux levels
s_min = float(options.flux.split(",")[0])
s_max = float(options.flux.split(",")[1])
s_step = float(options.flux.split(",")[2])
s = np.arange(s_min, s_max + 0.0001, s_step)
sdim = len(s)
slin = 10 ** (s + 3.0)  # convert to mJy in linear space

# Read input file with injected sources
ra_inj = []
dec_inj = []
for line in open(injected_sources):
    columns = line.split()
    if columns[0] != "#":
        ra_inj.append(float(columns[0]))
        dec_inj.append(float(columns[1]))

ra_inj = np.array(ra_inj)
dec_inj = np.array(dec_inj)
ninj = len(ra_inj)

# Read FITS catalogues with detected simulated sources
ra_det = []
dec_det = []
for i in range(0, sdim):
    f = detected_sources + "/flux*/det_source_list_flux" + str("%.4f" % s[i]) + ".fits"
    input_catalogue = glob.glob(f)
    if len(input_catalogue) < 1:
        print('Error: no match found for FITS catalogue "' + f + '". Aborting.')
        exit()
    elif len(input_catalogue) > 1:
        print('Error: multiple matches found for FITS catalogue "' + f + '". Aborting.')
        exit()
    ra, dec = read_catalogue(input_catalogue[0])
    ra_det.append(ra)
    dec_det.append(dec)

ra_det = np.array(ra_det)
dec_det = np.array(dec_det)

# Calculate completeness as a function of flux density and print results to file
cmp = []
cmp_err = []
f = open(options.output_root + ".txt", "w")
print("# log10(S_Jy) cmp_perc cmp_err_perc", file=f)
for i in range(0, sdim):
    ndet = len(ra_det[i])
    cmp.append(ndet / ninj * 100.0)
    cmp_err.append(np.sqrt(ndet) / ninj * 100.0)

    print("%.4f" % s[i], "%.2f" % cmp[i], "%.2f" % cmp_err[i], file=f)

# Plot completeness as a function of flux density
fontsize = 12
plt.xlabel("Flux density (mJy)", fontsize=fontsize)
plt.ylabel("Completeness (%)", fontsize=fontsize)
plt.grid(True)
plt.plot(slin, cmp, color="k", linewidth=0.7)
plt.scatter(slin, cmp, s=2.0, color="k")
plt.ylim(0, 100)
plt.xscale("log")
plt.tight_layout()
plt.savefig(options.output_root + ".pdf")
plt.close()

if options.template_map:
    print("Calculating completeness cube")
    # Update header of template FITS map
    template = fits.open(options.template_map)
    hdr = template[0].header
    hdr.set("CDELT3", s_step, before="CTYPE3")
    hdr.set("CRPIX3", 1.0, before="CTYPE3")
    hdr.set("CRVAL3", s_min, before="CTYPE3")
    hdr.set("CTYPE3", "log10[ S/(Jy/beam) ]", "")
    hdr.set("CUNIT3", "Jy/beam")
    hdr.set("BUNIT", "Completeness (per cent)")
    # If CDELT1 (the RA cell size) is positive, multiply by -1 to ensure that RA increases from right to left
    cdelt1 = hdr["CDELT1"]
    if cdelt1 > 0:
        hdr.set("CDELT1", -cdelt1)
    w = WCS(hdr)
    # Count total number of (RA,Dec) pixels
    ra1, dec1, dec2 = math.ceil(ra_min), math.ceil(dec_min), math.floor(dec_max)
    if ra_max < ra_min:
        ra2 = math.floor(ra_max + 360.0)
    else:
        ra2 = math.floor(ra_max)
    npix = (ra2 - ra1 + 1) * (dec2 - dec1 + 1)
    # Calculate completeness map at each flux level, blanking out pixels outside the specified region
    data = template[0].data
    shape = data.shape
    ydim = shape[1]
    xdim = shape[2]
    shape = (sdim, ydim, xdim)
    cmp_cube = np.empty(shape)
    cmp_cube[:] = np.nan
    step = 1.0
    f = step
    n = 0
    for j in range(0, ydim):
        for k in range(0, xdim):
            px = np.array([[k, j, 0]])
            t = w.wcs_pix2world(px, 0)
            ra, dec = t[0, 0], t[0, 1]
            if dec >= dec_min and dec <= dec_max:
                if (ra_max > ra_min and (ra >= ra_min and ra <= ra_max)) or (
                    ra_max < ra_min and (ra >= ra_min or ra <= ra_max)
                ):
                    ninj_rad = count_sources(ra, dec, ra_inj, dec_inj, options.rad)
                    for i in range(0, sdim):
                        ndet_rad = count_sources(
                            ra, dec, ra_det[i], dec_det[i], options.rad
                        )
                        cmp_cube[i, j, k] = ndet_rad / ninj_rad * 100.0
                    n = n + 1
                    p = (n / npix) * 100
                    if p >= f:
                        print("%.0f" % p, "% complete")
                        f = f + step
    out = template
    out[0].header = hdr
    out[0].data = cmp_cube
    out.writeto(options.output_root + ".fits", clobber=True)

