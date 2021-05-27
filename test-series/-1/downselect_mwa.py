import numpy as np
from astropy.io import ascii
import os


config_path = ('/lustre/aoc/projects/hera/Validation/test-neg1.1.0/' +
               'config_files/telescope_config')

# Read in the full MWA layout
mwa_layout_path = os.path.join(config_path, 'mwa_nocore_layout.csv')
ants = ascii.read(mwa_layout_path)
ants_e, ants_n = ants['E'], ants['N']

# Calculate the distance of all antennas from the center of the array
distance = np.sqrt(ants_e**2 + ants_n**2)

# Restrict baselines to HERA's maximum baseline length, including outriggers
max_bl_length = 876  # DeBoer+2017
print('Maximum baseline length: {} m'.format(max_bl_length))

max_radius = max_bl_length / 2
ants_inside = np.where(distance <= max_radius)[0]
ants_outside = np.where(distance > max_radius)[0]

# Write out the downsampled array layout
new_ants = {}
save_path = os.path.join(config_path, 'mwa_nocore_downsampled_layout.csv')
for key in ants.keys():
    new_ants[key] = ants[key][ants_inside]
ascii.write([new_ants[key] for key in new_ants.keys()],
            save_path, names=list(new_ants.keys()), overwrite=True)
