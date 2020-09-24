import numpy as np

import h5py

from RIMEz import beam_models, utils

nu_hz = 1e6*np.linspace(100., 200., 1024, endpoint=True)

beam_data_file_path = '/users/zmartino/zmartino/beam_data_from_nicolas/HERA_dipole_beams_spin1_model.h5'

beam_funcs = beam_models.model_data_to_spline_beam_func(
                beam_data_file_path,
                nu_hz,
                L_synth=180,
                indexed=True)

print "Beam function computed."

Omega, Omegapp = utils.beam_func_to_Omegas_pyssht(nu_hz, beam_funcs, L_use=200)

with h5py.File('HERA_dipole_Omegas.h5', 'w') as h5f:
    h5f.create_dataset('Omega', data=Omega)
    h5f.create_dataset('Omegapp', data=Omegapp)
    h5f.create_dataset('frequency_samples_hz', data=nu_hz)
