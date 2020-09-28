import h5py

realizations = [str(i) for i in range(20)]

with h5py.File('power_law_neg2_100_200MHz_realizations.h5', 'a') as h5f:
    for realization in realizations:
        del h5f[realization]['Ilm']
        print "Deleted realization",realization
