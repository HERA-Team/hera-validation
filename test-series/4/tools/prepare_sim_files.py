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
import sys
import yaml

from . import sim_prep

a = sim_prep.sim_prep_argparser()
# This will need to be updated for future validation tests,
# as naming conventions have changed since H1C!
obsfiles = sorted(glob.glob(os.path.join(a.obsdir, "*.HH.uvh5")))
if a.config is not None:
    with open(a.config, 'r') as cfg:
        systematics_params = yaml.load(cfg.read(), Loader=yaml.FullLoader)

sim_prep.prepare_sim_files(
    a.simfile, obsfiles, a.savedir, systematics_params=systematics_params, 
    save_truth=not a.skip_truth, clobber=a.clobber, verbose=a.verbose
)
