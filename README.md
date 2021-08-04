This is an set of codes to perform completness simulations on the GLEAM and GLEAM-X data. This code was originally developed by Thomas Franzen.

The basic procedure is described below. 

There are 4 steps to follow:





1) Generate a list of random RA and Dec positions where the simulated point sources will be injected.





Command line example:



sbatch --time=06:00:00 generate_pos.sh nsrc=30000 region=75,195,-40,-13 sep_min=5 output_dir=$GLEAMX/source_pos


The input parameters are explained at the top of the script. 10
 simulated sources / deg^2 should be sufficient, I think. Your mosaic covers about 3000 deg^2, so I set the number of simulated sources to 30,000. I
 set the minimum separation between the simulated sources to 5 arcmin in order to avoid introducing
 an artifical factor of confusion. The output directory will be created by the script. ‘$GLEAMX' is just the location of the project (in my case '/astro/mwasci/tfranzen/GLEAMX’,
 obviously you will want to use a different location).




2) Prepare file(s) listing the fluxes at which to measure the completeness


Command line example:


generate_fluxes.sh flux=-2.3,-0.5,0.1 nfiles=5 output_dir=$GLEAMX/fluxes


The input parameters are explained at the top of the script. The fluxes can be divided into multiple
 files to speed up the next step (then you can run the next step in parallel on each flux file). In
 this example, the fluxes are divided into 5 separate files. The script should only take a few seconds to run, so it can be run on the login node.




3) Inject point sources into the wideband image and
 run source finding on the simulated sources


Command line example:


sbatch --array=1-5 --time=12:00:00 inject_sources.sh input_map_dir=$GLEAMX/input_images input_sources=$GLEAMX/source_pos/source_pos.txt
 flux_dir=$GLEAMX/fluxes sigma=4.0 output_dir=$GLEAMX/inject


The input parameters are explained at the top of the script. There is one realisation per flux level (i.e.
 30,000 sources of the same flux are injected into the wideband mosaic). The script will loop over the fluxes listed in the flux file. The array parameter indicates which flux file to use. So 'array=1-5’ allows you to run the script in parallel on the 5 flux
 files generated in the previous step.


The script expects to find the following maps in the input map directory:
- XG_wideband.fits
- XG_wideband_rms.fits
- XG_wideband_bkg.fits
- XG_wideband_projpsf_psf.fits
Change this inside the script if needed.


Important: the wide-band image should be rescaled to account for ionospheric smearing, otherwise
 the completeness results will be too optimistic.


Please check the parameters used to run Aegean inside the script. They should be exactly the same as the parameters used to generate
 the final source catalogue. In this example, I set the source detection limit to 4 sigma, but you should obviously adjust that as required.




4) Calculate the overall completeness as a function of flux and generate completeness
 maps (one per flux level)


Command line example:


sbatch --time=12:00:00 make_cmp_map.sh injected_sources=$GLEAMX/source_pos/source_pos.txt detected_sources=$GLEAMX/inject flux=-2.3,-0.5,0.1 template_map=$GLEAMX/input_images/XG_wideband_projpsf_psf.fits
 region=75,195,-40,-13 rad=6 output_dir=$GLEAMX/results



The input parameters are explained at the top of the script. I set rad=6 but you can increase this
 if you want to apply more smoothing to the completeness maps. It may take several hours to create the completeness maps as it involves calculating lots of distances on a sphere.

