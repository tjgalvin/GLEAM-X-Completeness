#!/bin/bash -l

# Prepare flux files for inject_sources.sh

# Read input parameters
if [[ $1 ]] && [[ $2 ]] && [[ $3 ]]; then
    # Fluxes at which to measure completeness. Enter flux_min, flux_max, flux_interval, in Jy and in log space.
    # (e.g. for 21 equally-spaced fluxes in log space between 10 mJy and 1 Jy, enter -2,0,0.1)
    flux=$(echo $1 | awk -F"=" '{print $NF}')
    # Number of files in which to divide the fluxes
    nfiles=$(echo $2 | awk -F"=" '{print $NF}')
    # Output directory
    output_dir=$(echo $3 | awk -F"=" '{print $NF}')
else
    echo "Give me: flux nfiles output_dir"
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

# Create output directory
if [ -e $output_dir ]; then
    echo "Error: Output directory $output_dir already exists. Aborting."
    exit 1
else
    mkdir $output_dir
    cd $output_dir || exit 1
fi

# Write input parameters to file for record
cat >> input_parameters_generate_fluxes.txt <<EOPAR
flux = $flux
nfiles = $nfiles
output_dir = $output_dir
EOPAR

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
    for ((j=1; j<=($n); j++ )); do
        (( k++ ))
        printf "%0.4f\n" ${s[$k]} >> flux_list${i}.txt
    done
done

exit 0
