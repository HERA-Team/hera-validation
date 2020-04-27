"""Script for modifying and chunking a simulation.

This script serves as the first step in the hera_validation pipeline for 
the end-to-end validation test in support of the H1C power spectrum upper 
limits paper. This script takes a simulation file, a set of reference 
observation files, and modifies the simulation file contents such that:

    * The antenna array has been downselected and translated such that 
      its intersection with the H1C array yields the greatest number of 
      unique baselines. The simulation antenna positions have been 
      updated to exactly match the corresponding H1C antennas, and their 
      numbers (as well as all relevant metadata) have been updated 
      accordingly.

    * The simulation data has been downselected in time such that the 
      simulation covers the same LST range as the obsevation files. The 
      simulation data has also been rephased to match the observed LSTs.

After the simulation file is modified, it is chunked into files in a way 
that matches the chunking of the H1C files and follows the same naming 
convention with a minor modification.
"""
import glob
import os
import re
import sys
import yaml

import numpy as np
from . import sim_prep

a = sim_prep.sim_prep_argparser()
# This will need to be updated for future validation tests,
# as naming conventions have changed since H1C!
obsfiles = sorted(glob.glob(os.path.join(a.obsdir, "*.HH.uvh5")))
if a.config is not None:
    with open(a.config, 'r') as cfg:
        systematics_params = yaml.load(cfg.read(), Loader=yaml.FullLoader)

# Load in some things necessary for doing fileprep in chunks (which is 
# necessary for large files)
jd_pattern = re.compile('[0-9]{7}.[0-9]{5}')
start_jd, end_jd = (
    float(jd_pattern.findall(f)[0]) for f in (obsfiles[0], obsfiles[-1])
)
jd = int(start_jd)
sky_cmp = sim_prep._parse_filename_for_cmp(a.simfile)
new_config = os.path.join(a.savedir, f'{jd}.config.{sky_cmp}.yaml')
chunk_len = int(np.ceil(len(obsfiles) / a.Nchunks))
if 'gains' in systematics_params.keys():
    time_vary_params = systematics_params['gains'].get('time_vary_params', None)
else:
    time_vary_params = None
if a.Nchunks > 1 and time_vary_params is not None:
    file_duration = end_jd - start_jd
    center_jd = 0.5 * (start_jd + end_jd)
    for vary_mode, vary_params in time_vary_params.items():
        if vary_params.get('variation_ref_times', None) is None:
            vary_params['variation_ref_times'] = center_jd
        if vary_params.get('variation_timescales', None) is None:
            vary_params['variation_timescales'] = file_duration

for N in range(a.Nchunks):
    obsfile_chunk = obsfiles[N * chunk_len : (N+1) * chunk_len]
    sim_prep.prepare_sim_files(
        a.simfile, obsfile_chunk, a.savedir, systematics_params=systematics_params, 
        save_truth=not a.skip_truth, clobber=a.clobber, verbose=a.verbose
    )
    with open(new_config, 'r') as cfg:
        systematics_params = yaml.load(cfg.read(), Loader=yaml.FullLoader)

