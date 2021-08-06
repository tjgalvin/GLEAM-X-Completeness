#!/bin/bash
# Inject artificial sources into wide-band image used for source detection

start_time=$(date +%s)

# Read input parameters
if [[ $1 ]] && [[ $2 ]] && [[ $3 ]] && [[ $4 ]] && [[ $5 ]] && [[ $6 ]] ; then
    # Directory containing GLEAM-X mosaics (wide-band image used for source detection + rms, bkg and PSF maps)
    # Naming should be XG_wideband.fits, XG_wideband_rms.fits, XG_wideband_bkg.fits and XG_wideband_projpsf_psf.fits
    input_map_dir=$(echo $1 | awk -F"=" '{print $NF}')
    # File specifying the positions of the simulated sources
    input_sources=$(echo $2 | awk -F"=" '{print $NF}')
    # Directory containing simulated source flux files
    flux_dir=$(echo $3 | awk -F"=" '{print $NF}')
    # Sigma detection threshold used for source finding
    sigma=$(echo $4 | awk -F"=" '{print $NF}')
    # Output directory
    output_dir=$(echo $5 | awk -F"=" '{print $NF}')
    # Base file names for images
    imageset_name=$(echo $6 | awk -F"=" '{print $NF}')
else
    echo "Give me: input_map_dir input_sources flux_dir sigma output_dir imageset_name"
    exit 1
fi

# File containing list of simulated source fluxes
flux_list=$flux_dir/flux_list${SLURM_ARRAY_TASK_ID}.txt

# Create main output directory
if [ ! -e "${output_dir}" ]; then
    mkdir "${output_dir}"
fi

# Create flux output directory
if [ -e "${output_dir}/flux${SLURM_ARRAY_TASK_ID}" ]; then
    echo "Error: Output directory $output_dir/flux${SLURM_ARRAY_TASK_ID} already exists. Aborting."
    exit 1
else
    mkdir "${output_dir}/flux${SLURM_ARRAY_TASK_ID}"
    cd "${output_dir}/flux${SLURM_ARRAY_TASK_ID}" || exit
fi

# Set other input parameters
z=-26.7 # Dec at zenith in deg
search_rad=75 # Search radius to use when cross-matching sources to measure completeness, in arcsec
ncpus=48 # number of cores for running Aegean

# Write input parameters to file for record
cat >> input_parameters_inject_sources.txt <<EOPAR
input_map_dir = $input_map_dir
input_sources = $input_sources
flux_list = $flux_list
sigma = $sigma
z = $z
search_rad = $search_rad
ncpus = $ncpus
output_dir = $output_dir/flux${SLURM_ARRAY_TASK_ID}
imageset_name = $imageset_name
EOPAR

# Read source fluxes to simulate
i=0
while read line; do
    number=$(echo "$line" | grep -c "#")
    if [ "$number" != 1 ]; then
        (( i++ ))
        flux["$i"]=$(echo "$line")
    fi
done < $flux_list
nflux=$i

input_map="${input_map_dir}/${imageset_name}.fits" # Potentially this may miss ddmod.fits
input_map_rms="${input_map_dir}/${imageset_name}_rms.fits"
input_map_bkg="${input_map_dir}/${imageset_name}_bkg.fits"
input_map_psf="${input_map_dir}/${imageset_name}_projpsf_psf.fits"
for file in "${input_map}" "${input_map_rms}" "${input_map_bkg}" "${input_map_psf}" "${input_sources}"; do
    if [ ! -e "${file}" ]; then
        echo "Error: $file does not exist. Aborting."
        exit 1
    fi
done

# Define power function
pow(){
    echo "e($2 * l($1))" | bc -l
}

# ----------------------------------------------------------------
# Create real_and_sim_list.txt file for the image. This file lists the RA and Dec position of each real
# and simulated source in the image. A 3rd col gives the source type (1 for real, 0 for simulated).
# The list of real sources is obtained by running Aegean on the image.

# Run Aegean on real image
singularity exec \
"$CONTAINER" \
aegean \
--cores=$ncpus \
--out=aegean_list.txt \
--table=aegean_list.vot \
--noise="$input_map_rms" \
--background="$input_map_bkg" \
--seedclip="$sigma" \
--floodclip=4 \
--maxsummits=5 \
--psf="$input_map_psf" \
"$input_map"

rm -f aegean_list.txt

# Select RA and Dec columns in Aegean list of real sources; add type=1 col to indicate that these are real sources
singularity exec \
"$CONTAINER" \
stilts tpipe \
ifmt=votable \
in=aegean_list_comp.vot \
ofmt=ascii \
omode=out \
out=t \
cmd='addcol "type" 1' \
cmd='keepcols "ra dec type"'

rm -f aegean_list_comp.vot

# Select RA and Dec columns in list of simulated sources; add type=0 col to indicate these are simulated sources
singularity exec \
"$CONTAINER" \
stilts tpipe \
ifmt=ascii \
in="$input_sources" \
ofmt=ascii \
omode=out \
out=t2 \
cmd='addcol "type" 0' \
cmd='keepcols "ra dec type"'

# Concatenate real and simulated source lists
singularity exec \
"$CONTAINER" \
stilts tcat \
ifmt=ascii \
in=t \
in=t2 \
out=real_and_sim_list.txt

rm -f t t2

# ----------------------------------------------------------------

# Loop over fluxes in flux_list
for ((i=1; i<=($nflux); i++ )); do
    
    s=${flux[$i]}
    s_lin=$( pow 10 $s ) # convert flux to linear space
    
    # Get PSF size and blurring factor at the location of each simulated source
    singularity exec \
    "$CONTAINER" \
    "$MYCODE/calc_r_ratio_cmp.py" \
    --z=$z \
    --flux="$s_lin" \
    "$input_sources" \
    "$input_map" \
    "$input_map_psf" \
    "source_pos_flux${s}.txt"
    
    # Convert source position file to Aegean format (votable)
    # Since the map in which we will inject the point sources has already been rescaled to account for ionospheric smearing,
    # the peak fluxes of the injected sources should NOT be suppressed by the blurring factor
    # (i.e. the peak fluxes should be equal to the integrated fluxes)
    singularity exec \
    "$CONTAINER" \
    stilts tpipe \
    ifmt=ascii \
    ofmt=votable \
    omode=out \
    in=source_pos_flux${s}.txt \
    out=aegean_source_list.vot \
    cmd='addcol "island" $0' \
    cmd='addcol "source" "0"' \
    cmd='addcol "background" "0.0"' \
    cmd='addcol "local_rms" "0.0"' \
    cmd='addcol "ra_str" degreesToHms(RA,2)' \
    cmd='addcol "dec_str" degreesToDms(Dec,2)' \
    cmd='addcol "ra" RA' \
    cmd='addcol "err_ra" "0.0"' \
    cmd='addcol "dec" Dec' \
    cmd='addcol "err_dec" "0.0"' \
    cmd='addcol "peak_flux" -S' \
    cmd='addcol "err_peak_flux" "0.0"' \
    cmd='addcol "int_flux" -S' \
    cmd='addcol "err_int_flux" "0.0"' \
    cmd='addcol "a" bmaj' \
    cmd='addcol "err_a" "0.0"' \
    cmd='addcol "b" bmin' \
    cmd='addcol "err_b" "0.0"' \
    cmd='addcol "pa" bpa' \
    cmd='addcol "err_pa" "0.0"' \
    cmd='addcol "flags" "0"' \
    cmd='addcol "residual_mean" "0.0"' \
    cmd='addcol "residual_std" "0.0"' \
    cmd='addcol "uuid" "0"' \
    cmd='addcol "psf_a" "0.0"' \
    cmd='addcol "psf_b" "0.0"' \
    cmd='addcol "psf_pa" "0.0"' \
    cmd='delcols "RA Dec S bmaj bmin bpa R"'
    
    # Add simulated sources to real map
    singularity exec \
    "$CONTAINER" \
    AeRes \
    -c aegean_source_list.vot \
    -f "$input_map" \
    -r "sim_and_real_map_flux${s}.fits" \
    -m sim_map.fits
    
    rm -f sim_map.fits aegean_source_list.vot
    
    # Run Aegean on sim_and_real_map.fits (this is the real image + simulated sources); use existing rms and background images
    singularity exec \
    "$CONTAINER" \
    aegean \
    --cores=$ncpus \
    --out=aegean_SIM_list.txt \
    --table=aegean_SIM_list.vot \
    --noise="$input_map_rms" \
    --background="$input_map_bkg" \
    --seedclip="$sigma" \
    --floodclip=4 \
    --maxsummits=5 \
    --psf="$input_map_psf" "sim_and_real_map_flux${s}.fits"
    
    rm -f aegean_SIM_list.txt
    rm -f sim_and_real_map_flux${s}.fits
    
    # Match sources detected in the simulated image with the list of real & simulated sources for the image
    singularity exec \
    "$CONTAINER" \
    stilts tskymatch2 \
    in1=aegean_SIM_list_comp.vot \
    in2=real_and_sim_list.txt \
    out=match_list.txt \
    error=$search_rad \
    ra1='ra' \
    dec1='dec' \
    ra2='ra' \
    dec2='dec' \
    ifmt1=votable \
    ifmt2=ascii \
    find=best \
    join=1and2
    
    rm -f aegean_SIM_list_comp.vot
    
    # Select sources in match_list.txt that have type=0
    singularity exec \
    "$CONTAINER" \
    stilts tpipe \
    ifmt=ascii \
    in=match_list.txt \
    ofmt=fits \
    omode=out \
    out=det_source_list_flux${s}.fits \
    cmd='select "type == 0"' \
    cmd='addcol -after "ra_1" "ra" ra_1' \
    cmd='addcol -after "dec_1" "dec" dec_1' \
    cmd='addcol -after "ra_2" "ra_sim" ra_2' \
    cmd='addcol -after "dec_2" "dec_sim" dec_2' \
    cmd='addcol "flux_sim" '$s_lin'' \
    cmd='delcols "ra_1 dec_1 ra_2 dec_2 type"'
    
    rm -f match_list.txt
    
done

# -------------------------------------------------------------------------------------

end_time=$(date +%s)
duration=$(echo "$end_time-$start_time" | bc -l)
echo "Total runtime = $duration sec"

exit 0
