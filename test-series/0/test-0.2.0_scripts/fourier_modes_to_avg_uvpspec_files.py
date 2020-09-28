import os
import numpy as np

import hera_pspec as hps
import pyuvdata

import h5py
import time

from RIMEz import management, utils

from scipy.special import i0 as bessel_i0

########################################
######### Realizations to do ###########

realizations = [str(i) for i in range(0,50)]

########################################

script_dir = os.path.dirname(os.path.realpath(__file__))

data_dir = 'hera_hex37_100-200MHz_HERA_dipole_beam_pow_neg2_pspec'

base_name = os.path.join(script_dir, data_dir, 'realization_{}.rimezh5')
fourier_mode_files = [base_name.format(rel) for rel in realizations]


zdir = '/lustre/aoc/projects/hera/zmartino'

vis_data_dir = zdir + '/hera_validation/test-0.1/hera_hex37_100-200MHz_HERA_dipole_beam_pow_neg2_pspec'

def vis_data_path(realization):
    path = os.path.join(vis_data_dir, 'realization_{0}.uvh5'.format(realization))
    return path

# Load beam data
omegas_data_path = zdir + '/hera_validation/test-0.1/HERA_dipole_Omegas.h5'

with h5py.File(omegas_data_path, 'r') as h5f:
    Omega = h5f['Omega'].value
    Omegapp = h5f['Omegapp'].value

# Add a new window function to hera_pspec's aipy instance
kaiser6 = lambda x, L: bessel_i0(np.pi * 6 * np.sqrt(1-(2*x/(L-1) - 1)**2)) / bessel_i0(np.pi * 6)
hps.pspecdata.aipy.dsp.WINDOW_FUNC['kaiser6'] = kaiser6

def astropyPlanck15_for_hera_pspec():
    H0 = 67.74
    h = H0/100.

    Om_b = 0.02230/h**2.
    Om_c = 0.1188/h**2.
    Om_L = 0.6911
    Om_k = 1. - (Om_b + Om_c + Om_L)

    hps_cosmo = hps.conversions.Cosmo_Conversions(Om_L=Om_L,
                                                Om_b=Om_b,
                                                Om_c=Om_c,
                                                H0=H0,)
    return hps_cosmo

cosmo = astropyPlanck15_for_hera_pspec()

def get_VI_data(vis_data_path):
    uvd = pyuvdata.UVData()
    uvd.read_uvh5(vis_data_path)

    # one of these days...
    xx_integer = pyuvdata.utils.polstr2num('xx')
    yy_integer = pyuvdata.utils.polstr2num('yy')

    xx_ind = np.argwhere(uvd.polarization_array == xx_integer)[0][0]
    yy_ind = np.argwhere(uvd.polarization_array == yy_integer)[0][0]

    VI_data = uvd.data_array[:,:,:,xx_ind] + uvd.data_array[:,:,:,yy_ind]

    uvd.select(polarizations=(-5))
    uvd.polarization_array[0] = 1
    uvd.data_array = VI_data.reshape(VI_data.shape + (1,))

    return uvd

def write_realization_uvp(realization):
    uvdI = get_VI_data(vis_data_path(realization))

    hpsb = hps.pspecbeam.PSpecBeamFromArray(Omega, Omegapp, uvdI.freq_array[0], cosmo=cosmo)

    nu_e = 1420405751.7667 # Hz
    nu_hz = uvdI.freq_array[0]

    Jy_to_mK = hpsb.Jy_to_mK(nu_hz, pol='pI')

    Jy_to_mK_src = (nu_e/nu_hz) * Jy_to_mK

    uvdI.data_array *= Jy_to_mK_src[None,None,:,None]

    ds = hps.PSpecData(dsets=[uvdI, uvdI], wgts=[None,None], beam=hpsb)

    ds.dsets[0].vis_units = 'mK'
    ds.dsets[1].vis_units = 'mK'

    ant_pairs = [ant_pair for ant_pair in uvdI.get_antpairs() if ant_pair[0] != ant_pair[1]]

    edge_inds = [np.argmin(np.abs(nu_hz - nu_i)) for nu_i in 1e6*np.linspace(100.,200.,8, endpoint=True)]

    spw_ranges = zip(edge_inds, edge_inds[1:])

    uvp = ds.pspec(ant_pairs, ant_pairs,
                   dsets=(0,1),
                   pols=('pI', 'pI'),
                   spw_ranges=spw_ranges,
                   input_data_weight='identity',
                   norm='I',
                   taper='blackman-harris',
                   verbose=False,
                   little_h=False)
    blpair_group = [sorted(np.unique(uvp.blpair_array))]

    uvp_avg = uvp.average_spectra(blpair_groups=blpair_group, time_avg=True, inplace=False)

    uvp_avg.fold_spectra()

    uvp_save_path = os.path.join(vis_data_dir, 'uvp_avg_realization_{0}.h5'.format(realization))
    uvp_avg.write_hdf5(uvp_save_path, overwrite=True)

    return

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

    print "Computing pspec for realization {0}...".format(realizations[ii])

    write_realization_uvp(realizations[ii])

    if ii > 0:
        os.remove(uvdata_file_path)

    t2 = time.time()
    print "Took ", (t2 - t1)/60., "minutes."

end_time = time.time()

print "Elapsed time:", (end_time - start_time)/60., "minutes."
print "Done."
