import numpy as np
import datetime
import argparse
import os
import yaml
import ast
import astropy.units as u
from pyuvdata import UVData
from RIMEz import sky_models, beam_models, rime_funcs, utils, management


# --------------------------------------------------------
# Parse command line arguments

parser = argparse.ArgumentParser(description='Run RIMEz from a pyuvsim reference ' +
                                 'simulation config file.')
parser.add_argument('--sim', default='1.1', help='Reference simulation number ' +
                    '(default: 1.1)', required=False)
parser.add_argument('--beam', default='uniform', help='Beam type (default: ' +
                    'uniform)', required=False)
parser.add_argument('--lmax', help='Maximum L for the spherical harmonics ' +
                    '(default: calculated from frequencies and baselines)',
                    required=False)
parser.add_argument('--config_path', default='/lustre/aoc/projects/hera/Validation/' +
                    'test-neg1.1.0/config_files', help='Path to configuration files ' +
                    '(default: /lustre/aoc/projects/hera/Validation/test-neg1.1.0/' +
                    'config_files/)', required=False)
parser.add_argument('--ref_path', default='/lustre/aoc/projects/hera/Validation/' +
                    'test-neg1.1.0/ref_sims', help='Path to reference simulations ' +
                    '(default: /lustre/aoc/projects/hera/Validation/test-neg1.1.0/' +
                    'ref_sims/)', required=False)
parser.add_argument('--rimez_path', default='/lustre/aoc/projects/hera/Validation/' +
                    'test-neg1.1.0/rimez_sims', help='Path to save RIMEz simulations ' +
                    '(default: /lustre/aoc/projects/hera/Validation/test-neg1.1.0/' +
                    'rimez_sims/)', required=False)
args = parser.parse_args()

sim = args.sim
beam = args.beam.lower()
sim_name = '{}_{}'.format(sim, beam)
if sim == '1.2':
    sim_name = '{}_zenith_1freq_{}'.format(sim, beam)
if beam == 'hera':
    sim_name = sim_name.replace('hera', 'uvbeam')

ref_path = args.ref_path
rimez_path = args.rimez_path
config_path = args.config_path

# --------------------------------------------------------
# Observation parameters

param_path = ''
if sim == '1.1':
    param_path = os.path.join(config_path, 'obsparam_ref_' +
                              '{}_downsampled_mwa.yaml'.format(sim))
elif sim in ['1.2', '1.3']:
    param_path = os.path.join(config_path, 'obsparam_ref_' +
                              '{}.yaml'.format(sim_name))
print(param_path)
assert os.path.exists(param_path), ('Configuration file for ' +
                                    '{} not found'.format(sim_name))
print('Reading simulation parameters from {}'.format(param_path))

ref_data_path = os.path.join(ref_path, 'ref_{}.uvh5'.format(sim_name))
assert os.path.exists(ref_data_path), ('Reference simulation ' +
                                      '{} not found'.format(sim_name))
print('Loading reference simulation from {}'.format(ref_data_path))

ref_uvd = UVData()
ref_uvd.read(ref_data_path, run_check_acceptability=False)

with open(param_path) as param_file:
    param_dict = yaml.load(param_file, Loader=yaml.FullLoader)
freq_params = param_dict['freq']
source_params = param_dict['sources']
tel_params = param_dict['telescope']
time_params = param_dict['time']

# Times
integration_time = time_params['integration_time']
time_sample_jds = np.unique(ref_uvd.time_array)

# Frequencies
freq_samples_hz = ref_uvd.freq_array[0]  # First spectral window
# Channel width (for saving the RIMEz simulation later)
channel_width = 'derived'
if len(freq_samples_hz) == 1:
    channel_width = 0.

# Location
tel_param_path = os.path.join(config_path,
                              tel_params['telescope_config_name'])
with open(tel_param_path) as tel_param_file:
    tel_params_dict = yaml.load(tel_param_file, Loader=yaml.FullLoader)
location = ast.literal_eval(tel_params_dict['telescope_location'])
array_latitude = np.deg2rad(location[0])
array_longitude = np.deg2rad(location[1])
array_height = location[2]

# Antennas
bl_param_path = os.path.join(config_path, tel_params['array_layout'])
ant_pos_meters = np.genfromtxt(bl_param_path, skip_header=1,
                               usecols=(3, 4, 5))
ant_pairs_used, u2a, a2u = utils.get_minimal_antenna_set(ant_pos_meters,
                                                         precision=3)
ant_beam_function_map = np.zeros(ant_pos_meters.shape[0], dtype=np.int64)

if beam == 'uniform':
    beam_funcs = beam_models.uniform
elif beam == 'gauss':
    sigma = tel_params_dict['sigma']
    beam_funcs = beam_models.make_gaussian_dipole(sigma, pol=False)
elif beam == 'airy':
    diameter = tel_params_dict['diameter']
    # Config files give the diameter, the RIMEz function takes
    # the radius (we think...)
    beam_funcs = beam_models.make_airy_dipole(diameter / 2., pol=False)
elif beam == 'hera':
    # Assuming H1C HERA dipole feed (Vivaldi beam files do exist)
    beam_path = ('/lustre/aoc/projects/hera/Validation/' +
                 'HERA_beams/HERA_dipole_beams_spin1_model.h5')
    beam_funcs = beam_models.model_data_to_spline_beam_func(
                     beam_path,
                     freq_samples_hz,
                     L_synth=180,
                     indexed=True)

# This is something related to the Fourier series cutoff
integral_kernel_cutoff = utils.kernel_cutoff_estimate(
                            np.amax(list(u2a.keys())),
                            np.amax(freq_samples_hz),
                            width_estimate=100)

parameters = {
    'array_latitude': array_latitude,
    'array_longitude': array_longitude,
    'array_height': array_height,
    'initial_time_sample_jd': time_sample_jds[0],
    'integration_time': integration_time,
    'frequency_samples_hz': freq_samples_hz,
    'antenna_positions_meters': ant_pos_meters,
    'antenna_pairs_used': ant_pairs_used,
    'antenna_beam_function_map': ant_beam_function_map,
    'integral_kernel_cutoff': integral_kernel_cutoff,
    'time_sample_jds': time_sample_jds
}

# --------------------------------------------------------
# Sky parameters

# RA/dec
source_param_path = os.path.join(config_path,
                                 source_params['catalog'])
cols_to_use = (1, 2, 3, 4)  # RA, dec, flux, freq
if sim == '1.2':
    cols_to_use = (1, 2, 3)  # RA, dec, flux
sources = np.genfromtxt(source_param_path, skip_header=1,
                        usecols=cols_to_use)
intensities = sources[:, 2:3].T  # Get the shape right for RIMEz
# Broadcast to the number of frequencies we're simulating
# (assuming a uniform source)
intensities = np.broadcast_to(intensities, (len(freq_samples_hz),
                              len(sources)))

# ICRS
ra, dec = sources[:, 0] * u.deg, sources[:, 1] * u.deg
ra, dec = ra.to(u.rad), dec.to(u.rad)

# --------------------------------------------------------
# Simulate the sky and calculate visibilities

# Set maximum L value
if args.lmax:
    lmax = int(args.lmax)
else:
    lmax = integral_kernel_cutoff

Slm = sky_models.threaded_point_sources_harmonics(intensities, ra, dec, L=lmax)
Slm = Slm.reshape(Slm.shape + (1,))
VC = management.VisibilityCalculation(parameters, beam_funcs, Slm)
VC.compute_time_series()

if not os.path.exists(rimez_path):
    os.makedirs(rimez_path)
save_path = os.path.join(rimez_path, 'rimez_{}_lmax{}.uvh5'.format(sim_name, lmax))
print('Saving RIMEz simulation to {}'.format(save_path))

beam_name = beam
if beam == 'gauss':
    beam_name = 'Gaussian'
elif beam == 'airy':
    beam_name = 'Airy'
elif beam == 'hera':
    beam_name = 'HERA'
time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
history = ('Created ' + time + ' with RIMEz to match first reference ' +
           'simulation ' + sim + ' with ' + beam_name + ' beam') 
VC.write_uvdata_time_series(save_path, channel_width=channel_width,
                            clobber=True, instrument='RIMEz',
                            telescope_name='HERA', history=history)
