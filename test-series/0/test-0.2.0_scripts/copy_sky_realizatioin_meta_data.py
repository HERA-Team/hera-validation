import h5py

realizations = [str(i) for i in range(20)]

with h5py.File('power_law_neg2_100_200MHz_realizations.h5', 'a') as h5f_in:
    with h5py.File('power_law_neg2_100_200MHz_realizations_meta_data.h5', 'w') as h5f_out:

        for key in h5f_in:
            if not (key in realizations):
                h5f_out.create_dataset(key, data=h5f_in[key].value)

        for realization in realizations:
            dset = h5f_out.create_group(realization)
            dset.create_dataset('numpy_random_seed', data=h5f_in[realization]['numpy_random_seed'].value)
            
