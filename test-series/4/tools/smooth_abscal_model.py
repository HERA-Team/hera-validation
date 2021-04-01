"""
Smooth reference files for absolute calibration.
"""
import numpy as np
from astropy import units
from pyuvdata import UVData
import glob
import hera_cal as hc
import os

from . import sim_prep

a = sim_prep.smooth_abscal_model_argparser()
clean_params = vars(a)

# get datafiles
dfiles = sorted(glob.glob(os.path.join(a.datadir, a.datafile)))

# iterate over dfiles
for df in dfiles:
    outfile = df.replace(".uvh5", ".M.uvh5")

    # load filter
    DF = hc.delay_filter.DelayFilter(df, filetype='uvh5')
    DF.read(polarizations=['ee', 'nn'])

    # run filter
    DF.run_filter(min_dly=a.min_dly, window=a.window, alpha=a.alpha, tol=a.tol, gain=a.gain,
                  maxiter=a.maxiter, skip_wgt=a.skip_wgt, horizon=a.horizon, standoff=a.standoff,
                  edgecut_low=a.edgecut_low, edgecut_hi=a.edgecut_hi)

    history = "Delay CLEANed with: "
    for p in ['tol', 'maxiter', 'window', 'alpha', 'skip_wgt', 'min_dly', 'horizon', 'standoff',
              'edgecut_low', 'edgecut_hi', 'gain']:
        history += "{}={}; ".format(p, clean_params[p])

    # write clean model
    print("...Writing {}".format(outfile))
    DF.write_data(DF.clean_model, outfile, add_to_history=history, overwrite=True)

