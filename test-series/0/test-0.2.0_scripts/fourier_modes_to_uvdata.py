import os
import numpy as np

import time

from RIMEz import management, utils

realizations = [str(i) for i in range(10,20)]

script_dir = os.path.dirname(os.path.realpath(__file__))

data_dir = 'hera_hex37_100-200MHz_HERA_dipole_beam_pow_neg2_pspec'

base_name = os.path.join(script_dir, data_dir, 'realization_{}.rimezh5')
fourier_mode_files = [base_name.format(rel) for rel in realizations]

start_time = time.time()

for ii in range(len(fourier_mode_files)):

    t1 = time.time()

    print "Generating uvh5 for realization_{0}...".format(realizations[ii])
    file_path = fourier_mode_files[ii]
    VC = management.VisibilityCalculation(restore_file_path=file_path)

    delta_era_axis = np.linspace(0., 2*np.pi, 8052, endpoint=False)
    jd0 = VC.initial_time_sample_jd
    era_axis = delta_era_axis + utils.JD2era_tot(jd0)
    time_sample_jds = np.array(map(lambda era: utils.era2JD(era, jd0), era_axis))
    integration_time = 10.73

    VC.compute_time_series(time_sample_jds=time_sample_jds, integration_time=integration_time)

    uvdata_file_name = 'realization_{}.uvh5'.format(realizations[ii])
    uvdata_file_path = os.path.join(script_dir, data_dir, uvdata_file_name)

    history = "time series generated from {}".format(file_path)
    uvdata_extra_args = {
    'instrument': 'RIMEz calculation',
    'telescope_name': 'mock-HERA',
    'history': history,
    'object_name': "Gaussian cosmological field with P(k) = A0 k^-2, A0 = 1e1*(1./0.2)**-2 (Kelvin), realization_{0}".format(realizations[ii]),
    }

    VC.write_uvdata_time_series(uvdata_file_path, clobber=True, **uvdata_extra_args)

    t2 = time.time()
    print "Took", (t2-t1)/60.,"minutes."

end_time = time.time()
print "Elapsed time:", (end_time - start_time)/60., "minutes."
