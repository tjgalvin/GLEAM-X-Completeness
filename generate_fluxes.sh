#!/bin/bash -l

# Prepare flux files for inject_sources.sh

if [[ -z ${MYCODE} ]]
then
    echo "Error: The MYCODE variable is missing. Exiting."
    exit 1
fi

# module load singularity
echo "${SINGULARITY_BINDPATH}"

# Read input parameters
if [[ $1 ]] && [[ $2 ]] && [[ $3 ]] && [[ $4 ]] && [[ $5 ]] && [[ $6 ]]; then
    
    # Number of random source positions to simulate
    nsrc=$(echo $1 | awk -F"=" '{print $NF}')
    # Region of sky - enter ra_min, ra_max, dec_min, dec_max, in deg.
    # For example, for 60 < RA < 300 enter: 60,300,-90,90
    # For RA < 60 or RA > 300, and -40 < Dec < -10, enter: 300,60,-40,-10
    region=$(echo $2 | awk -F"=" '{print $NF}')
    # Minimum separation between simulated sources, in arcmin
    sep_min=$(echo $3 | awk -F"=" '{print $NF}')
    
    # Fluxes at which to measure completeness. Enter flux_min, flux_max, flux_interval, in Jy and in log space.
    # (e.g. for 21 equally-spaced fluxes in log space between 10 mJy and 1 Jy, enter -2,0,0.1)
    flux=$(echo $4 | awk -F"=" '{print $NF}')
    # Number of files in which to divide the fluxes
    nfiles=$(echo $5 | awk -F"=" '{print $NF}')
    # Output directory
    output_dir=$(echo $6 | awk -F"=" '{print $NF}')
else
    echo "Give me: nsrc region sep_min flux nfiles output_dir"
    exit 1
fi

# Read flux parameter and perform some basic checks
flux_min=$(echo ${flux//,/ } | awk '{print $1}')
flux_max=$(echo ${flux//,/ } | awk '{print $2}')
flux_step=$(echo ${flux//,/ } | awk '{print $3}')
if [ $(echo "$flux_min > $flux_max"|bc) -eq 1 ]; then
    echo "Error: flux_min > flux_max. Aborting."
    exit 1
fi
if [ $(echo "$flux_step <= 0"|bc) -eq 1 ]; then
    echo "Error: flux_step <= 0. Aborting."
    exit 1
fi

# Calculate flux levels
i=1; s[$i]=$flux_min
while [ $(echo "${s[$i]} < $flux_max" | bc -l) -eq 1 ]; do
    (( i++ ))
    s[$i]=$(echo "${s[$i-1]}+$flux_step" | bc -l)
done
nflux=$i
if [ $(echo "$nfiles > $nflux"|bc) -eq 1 ]; then
    echo "Error: nfiles > number of fluxes. Aborting."
    exit 1
fi

pos_outdir="${output_dir}/source_pos"
flux_outdir="${output_dir}/fluxes"

if [[ -e "${pos_outdir}" ]]; then
    echo "Error: Output directory ${pos_outdir} already exists. Aborting."
    exit 1
fi

if [[ -e "${flux_outdir}" ]]; then
    echo "Error: Output directory ${flux_outdir} already exists. Aborting."
    exit 1
fi

if [[ ! -e "${output_dir}" ]]
then
    echo "Creating ${output_dir}"
    mkdir "${output_dir}"
fi

basedir=$(pwd)

# Write input parameters to file for record
cat >> "${output_dir}"/input_parameters_generate_sources.txt <<EOPAR
nsrc = $nsrc
region = $region
sep_min = $sep_min
flux = $flux
nfiles = $nfiles
output_dir = $output_dir
pos_outputdir = $pos_outdir
flux_outputdir = $flux_outdir
EOPAR

mkdir "${pos_outdir}"
cd "${pos_outdir}" || exit 1

# Run Python script to generate RA and Dec positions
singularity exec \
-B "$MYCODE" \
"$CONTAINER" \
"$MYCODE/generate_pos.py" \
--nsrc="$nsrc" \
--region="$region" \
--sep-min="$sep_min" \
source_pos.txt

cd "${basedir}" || exit 1
mkdir "${flux_outdir}"
cd "${flux_outdir}" || exit 1

# Generate flux files
quotient=$((nflux/nfiles))
remainder=$((nflux%nfiles))
k=0
for ((i=1; i<=($nfiles); i++ )); do
    echo "# log10(flux/Jy)" > flux_list${i}.txt
    if [ $(echo "$i <= $remainder"|bc) -eq 1 ]; then
        n=$(echo "$quotient+1" | bc -l)
    else
        n=$quotient
    fi
    echo "Writing flux_list${i}.txt"
    for ((j=1; j<=($n); j++ )); do
        (( k++ ))
        printf "%0.4f\n" ${s[$k]} >> flux_list${i}.txt
    done
done

exit 0
