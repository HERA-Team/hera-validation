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
from pyuvdata import UVData
from . import sim_prep

a = sim_prep.sim_prep_argparser()
# This will need to be updated for future validation tests,
# as naming conventions have changed since H1C!
obsfiles = sorted(glob.glob(os.path.join(a.obsdir, "*.HH.uvh5")))
if a.lst_min != 0 or a.lst_max != 24:
    obsfiles = sim_prep._apply_lst_cut(
        obsfiles, a.lst_min, a.lst_max, uniformly_sampled=a.uniformly_sampled
    )
if a.config is not None:
    with open(a.config, 'r') as cfg:
        systematics_params = yaml.load(cfg.read(), Loader=yaml.FullLoader)
else:
    systematics_params = {}

# Load in some things necessary for doing fileprep in chunks (which is 
# necessary for large files)
jd_pattern = re.compile('[0-9]{7}.[0-9]{5}')
start_jd = float(jd_pattern.findall(obsfiles[0])[0])
uvd = UVData()
uvd.read(obsfiles[-1], read_data=False)
end_jd = uvd.time_array.max()
sky_cmp = sim_prep._parse_filename_for_cmp(a.simfile)
chunk_len = int(np.ceil(len(obsfiles) / a.Nchunks))
for param in ("gains", "reflections", "reflection_spectrum"):
    time_vary_params = systematics_params.get(param, {}).get('time_vary_params', None)
    if a.Nchunks > 1 and time_vary_params is not None:
        file_duration = float(end_jd - start_jd)
        center_jd = 0.5 * float(start_jd + end_jd)
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
    jd = jd_pattern.findall(obsfile_chunk[0])[0]
    if a.config is not None:
        new_config = os.path.join(a.savedir, f'{jd}.config.{sky_cmp}.yaml')
        with open(new_config, 'r') as cfg:
            systematics_params = yaml.load(cfg.read(), Loader=yaml.FullLoader)
        # Ensure that the noise realization is unique.
        if 'noise' in systematics_params.keys():
            systematics_params['noise']['seed'] = 'random'
