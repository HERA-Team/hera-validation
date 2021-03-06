import os
import inspect
import time

import numpy as np
import h5py

from scipy import linalg
from astropy import cosmology

from redshifted_gaussian_fields import generator

from RIMEz import sky_models, beam_models, rime_funcs, utils, management

data_dir = 'hera_hex37_100-200MHz_HERA_dipole_beam_pow_neg2_pspec'
if os.path.exists(data_dir) is False:
    os.mkdir(data_dir)

TEST_MODE = False

#############################################################
########## Set indices of realizations to generate ##########
#############################################################

realizations_use = [str(ii) for ii in range(35,40)]

#############################################################

if TEST_MODE:
    realizations_use = [str(ii) for ii in range(900,902)]

array_latitude = utils.HERA_LAT
array_longitude = utils.HERA_LON
array_height = utils.HERA_HEIGHT

JD_INIT = 2458529.500000 # no particular signifigance
era_init = 2*np.pi - array_longitude + 2*2*np.pi/86164. # an ERA such that the LST very close to zero, and > 0

initial_time_sample_jd = utils.era2JD(era_init, JD_INIT)

integration_time = 1e-9

frequency_samples_hz = 1e6*np.linspace(100.,200.,1024, endpoint=True)

if TEST_MODE:
    frequency_samples_hz = frequency_samples_hz[::1024/3]

if TEST_MODE:
    antenna_positions_meters = utils.generate_hex_positions(u_lim=2, v_lim=2, w_lim=1)
else:
    hex_lim = 4 # 37-element hex array
    antenna_positions_meters = utils.generate_hex_positions(u_lim=hex_lim, v_lim=hex_lim, w_lim=hex_lim)

antenna_pairs_used, u2a, a2u = utils.get_minimal_antenna_set(antenna_positions_meters, precision=3)

antenna_beam_function_map = np.zeros(antenna_positions_meters.shape[0], dtype=np.int64)

integral_kernel_cutoff = utils.kernel_cutoff_estimate(
                            np.amax(u2a.keys()),
                            np.amax(frequency_samples_hz),
                            width_estimate=125)

parameters = {
    'array_latitude': array_latitude,
    'array_longitude': array_longitude,
    'array_height': array_height,
    'initial_time_sample_jd': initial_time_sample_jd,
    'integration_time': integration_time,
    'frequency_samples_hz': frequency_samples_hz,
    'antenna_positions_meters': antenna_positions_meters,
    'antenna_pairs_used': antenna_pairs_used,
    'antenna_beam_function_map': antenna_beam_function_map,
    'integral_kernel_cutoff': integral_kernel_cutoff,
}

beam_data_file_path = '/users/zmartino/zmartino/beam_data_from_nicolas/HERA_dipole_beams_spin1_model.h5'

beam_funcs = beam_models.model_data_to_spline_beam_func(
                beam_data_file_path,
                frequency_samples_hz,
                L_synth=180,
                indexed=True)

beam_function_source = 'spin1_model_data'
spin1_model_data = 'HERA_dipole_beams_spin1_model'
parameters['beam_function_source'] = np.string_(beam_function_source)
parameters['spin1_model_data'] = np.string_(spin1_model_data)

print("Beam function computed.")

realizations_file = 'power_law_neg2_100_200MHz_realizations.h5'
realizations_file_path = os.path.join(os.getcwd(), realizations_file)

def get_realization(index_str, L_cut):
    seed = np.random.randint(1, 100000)
    np.random.seed(seed)

    nu_hz = frequency_samples_hz
    del_nu = np.diff(nu_hz)[0]

    cosmo = cosmology.Planck15

    A0 = 1e1 * (1./0.2)**-2.

    ell_max = L_cut

    barC = generator.pow_neg2_CAPS(cosmo, 1e-6*nu_hz, 1e-6*del_nu, A0, ell_max)

    eig_vals = np.zeros(barC.shape[:2])
    eig_vecs = np.zeros_like(barC)

    for ll in range(barC.shape[0]):
        eig_vals[ll], eig_vecs[ll] = linalg.eigh(barC[ll])

    alm = generator.sh_realization_healpy_index(eig_vals, eig_vecs, seed)

    Ilm = sky_models.hp2ssht_index(alm.T)

    Ilm = Ilm.reshape(Ilm.shape + (1,))

    git_hash = generator.get_commit_hash()
    with h5py.File(realizations_file_path, 'w') as h5f:
        grp = h5f.create_group(index_str)
        grp.create_dataset('numpy_random_seed', data=seed)
        grp.create_dataset('git_hash', data=git_hash)

    return Ilm

sky_data_source = realizations_file_path
parameters['sky_data_source'] = np.string_(sky_data_source)

# this is probably not needed, but why not? just in case...
this_file = inspect.getmodule(inspect.currentframe())
script_code = np.string_(inspect.getsource(this_file))

parameters['script_code'] = script_code

start_time = time.time()

for ii in range(len(realizations_use)):
    print "Working on realization", realizations_use[ii], "..."

    time1 = time.time()

    index_str = realizations_use[ii]
    Slm = get_realization(index_str, integral_kernel_cutoff+1)

    VC = management.VisibilityCalculation(parameters, beam_funcs, Slm)

    VC.compute_fourier_modes()

    output_file = os.path.join(data_dir, 'realization_{}.rimezh5'.format(index_str))
    VC.write_visibility_fourier_modes(output_file, overwrite=True)

    time2 = time.time()
    print "Finished realization {0}. Elapsed time:".format(index_str), (time2-time1)/60.

end_time = time.time()

print "Total elapsed time:", (end_time - start_time)/60.
print "Done."
