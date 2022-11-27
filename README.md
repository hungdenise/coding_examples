## Searching for the Yarkovsky Effect in Asteroid Orbits
[make_syn_obs.py](https://github.com/hungdenise/coding_examples/blob/main/yarkovsky/make_syn_obs.py) is a Python script to create options files for generating ephemerides in OrbFit assuming a user-defined a priori A2 value 
and then use the ephemerides to create synthetic observations of varying arc length.

[963.oop](https://github.com/hungdenise/coding_examples/blob/main/yarkovsky/963.oop) is the example template options file, [13936_apriori_A1.oep](https://github.com/hungdenise/coding_examples/blob/main/yarkovsky/13936_apriori_A1.oep) is an example ephemerides file, and [13936_A1_as0.03_0_arc25.rwo](https://github.com/hungdenise/coding_examples/blob/main/yarkovsky/13936_A1_as0.03_0_arc25.rwo) is an 
example output synthetic observations file.

-----------------------------------------------------------------------------------------------------------------------------------------
## Finding Galaxy Clusters in Voronoi tessellation Monte-Carlo Mapping
[regmatch_vuds.pro](https://github.com/hungdenise/coding_examples/blob/main/vmc/regmatch_vuds.pro) is an IDL script for linking Source Extractor detections across neighboring redshift slices, where each redshift slice is a separate file. [zlinks.py](https://github.com/hungdenise/coding_examples/blob/main/vmc/zlinks.py) is a Python script that fits Gaussians to the redshifts and fluxes of the linked chains to obtain a list of candidate clusters.

[COSMOS_d3_a100_db01_567.cat](https://github.com/hungdenise/coding_examples/blob/main/vmc/COSMOS_d3_a100_db01_567.cat) is an example Source Extractor detections file and [regmatch_COSMOS_d3_a100_db01_2d.txt](https://github.com/hungdenise/coding_examples/blob/main/vmc/regmatch_COSMOS_d3_a100_db01_2d.txt) is an example full linked output file for the COSMOS field using a DETECT_THRESH sigma of 3, a DETECT_MINAREA of 100 pixels, a DEBLEND_MINCONT of 0.1, and a linking radius of 2 Mpc.

-----------------------------------------------------------------------------------------------------------------------------------------
