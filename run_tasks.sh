#!/bin/env bash

# A simple helper script to string together the completeness tasks that need to be run.

cluster=replaceme
nsrc=30000
region=75,210,-43,-13
flux=-3,-0.5,0.1
nfiles=6
sep_min=5
outdir=/replace/me
imageset='XG_170-231MHz'

imageset_dir=/replace/me/with/images

export GLEAMX="${outdir}"
export MYCODE=
export NCPUS=38
export CONTAINER=/replace/me


if [[ -z ${MYCODE} ]]
then
    echo "Error. Completeness code directory not set. Exiting. "
    return 1
fi

if [[ ! -d $outdir ]]
then
    echo "Making directory ${outdir}"
    mkdir -p "${outdir}"
fi

#TODO: See how well this works with symlinks. Need to be sure the container can follow them.
for suffix in "" "_bkg" "_rms" "_projpsf_psf"
do
    if [[ -e "${imageset_dir}/${imageset}${suffix}.fits" ]]
    then
        cp -v "${imageset_dir}/${imageset}${suffix}.fits" "${GLEAMX}/inputimages"
    else
        echo "Could not find ${imageset_dir}/${imageset}${suffix}.fits. Exiting. "
        return 1
    fi
done

"$MYCODE"/generate_fluxes.sh \
$nsrc \
$region \
$sep_min \
$flux \
$nfiles \
$outdir

# We will be blocking until we are finished
srun \
--array=1-$nfiles \
--time:24:00:00 \
--ntasks-per-node=$NCPUS \
--export=ALL \
-o="${outdir}/inject_source.o%A" \
-e="${outdir}/inject_source.e%A" \
"$MYCODE/inject_sources.sh" \
input_map_dir="${GLEAMX}/input_images" \
input_sources="${GLEAMX}/source_pos/source_pos.txt" \
flux_dir="${GLEAMX}/fluxes" \
sigma=4.0 \
output_dir="${GLEAMX}/inject" \
imageset_name="${imageset}"

"$MYCODE"/make_cmp_map.sh \
injected_sources="${GLEAMX}/source_pos/source_pos.txt" \
detected_sources="${GLEAMX}/inject" \
flux="$flux" \
template_map="${GLEAMX}/input_images/XG_wideband_projpsf_psf.fits" \
region="${region}" \
rad=6 \
output_dir="${GLEAMX}/results"


