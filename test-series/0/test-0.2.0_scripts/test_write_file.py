import os
import h5py

zdir = '/lustre/aoc/projects/hera/zmartino'

vis_data_dir = zdir + '/hera_validation/test-0.1/hera_hex37_100-200MHz_HERA_dipole_beam_pow_neg2_pspec'
realization = '11'
uvp_save_path = os.path.join(vis_data_dir, 'uvp_realization_{0}.h5'.format(realization))

with h5py.File(uvp_save_path, 'w') as h5f:
    h5f.create_dataset('test', data=5.)


with h5py.File(uvp_save_path, 'r') as h5f:
    print h5f['test'].value
