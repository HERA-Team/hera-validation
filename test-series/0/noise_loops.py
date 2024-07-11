# Load the things
from pyuvdata import UVData
import hera_pspec as hp
import numpy as np
import matplotlib.pyplot as plt
import copy, os, itertools, inspect
from hera_pspec.data import DATA_PATH

def Powerg_m(n, mu, sigma):
    """
    Parameters

    n: number of noise realizations
    mu: mean of Gaussian noise
    sigma: standard deviation of Gaussian noise. The true values of sigma shown in the plots are scaled by 0.5

    """

    # P(k) with Gaussian noise for all 491 time samples, all 100 delays, and all 3 baselines
    Powers_g = []
    # Same noise is added each time
    np.random.seed(1)
    for i in range(n):
        # select the data file to load
        dfile = os.path.join('eorsky_3.00hours_Nside128_sigma0.03_fwhm12.13_uv.uvh5')
        # Load into UVData objects
        uvd = UVData()
        uvd.read_uvh5(dfile)
        # Add random Gaussian noise, repeat code
        uvdg = UVData()
        uvdg.read_uvh5(dfile)
        
        # (From the PSpecBeam Tutorial)
        # Each beam is defined over a frequency interval:
        beam_freqs = uvd.freq_array[0]

        # Create a new Gaussian beam object with full-width at half-max. of 12.13 degrees, converted to radians
        beam_gauss = hp.PSpecBeamGauss(fwhm=12.13*np.pi/180, beam_freqs=beam_freqs)
        
        # Apply unit conversion factor to UVData
        # The expression [None, None, :, None] reshapes the conversion factor into the same shape as as the data_array
        uvd.data_array *= beam_gauss.Jy_to_mK(np.unique(uvd.freq_array))[None, None, :, None]
        uvdg.data_array *= beam_gauss.Jy_to_mK(np.unique(uvdg.freq_array))[None, None, :, None]
        
        # Add the noise manually
        # (until this point, uvdg and uvd are the same)
        d = uvd.data_array
        uvdg.data_array += 0.5*(np.random.normal(mu, sigma, size=d.shape) + 1.j*np.random.normal(mu, sigma, size=d.shape))
        
        # Instantiate a Cosmo Conversions object
        # we will need this cosmology to put the power spectra into cosmological units
        cosmo = beam_gauss.cosmo
        print(cosmo)
        
        # Gaussian beam object
        uvb = beam_gauss
        
        # slide the time axis of uvd by one integration
        uvdg1 = uvdg.select(times=np.unique(uvdg.time_array)[:-1:2], inplace=False)
        uvdg2 = uvdg.select(times=np.unique(uvdg.time_array)[1::2], inplace=False)
        
        # Create a new PSpecData object, and don't forget to feed the beam object
        dsg = hp.PSpecData(dsets=[uvdg1, uvdg2], wgts=[None, None], beam=uvb)
        
        # Because the LST integrations are offset by more than ~15 seconds we will get a warning
        # but this is okay b/c it is still **significantly** less than the beam-crossing time and we are using short
        # baselines...
        # here we phase all datasets in dsets to the zeroth dataset
        dsg.rephase_to_dset(0)
        
        # change units of UVData objects
        dsg.dsets[0].vis_units = 'mK'
        dsg.dsets[1].vis_units = 'mK'
        
        # Specify which baselines to include
        # For the simulated data, baselines are (0,11), (0,12), (11,12)
        baselines = [(0, 11), (0, 12), (11, 12)]
        
        # we will use the baselines list to produce 491 power spectra 
        # whose data will be drawn from the dsets[0] and dsets[1]
        # across two spectral windows with identity weighting 
        # and a blackman-harris taper.
        
        # Take the full spw range: (0, 384)
        uvpg = dsg.pspec(baselines, baselines, (0, 1), [('pI', 'pI')], spw_ranges=[(0, 384)], input_data_weight='identity',
                         norm='I', taper='blackman-harris', verbose=True)
        
        # Spectral window
        spw = 0
        # keys for 3 baselines
        key_bl = []
        # Power with Gaussian noise
        powerg_bl = []
        for j in range(len(baselines)):
            key_bl.append((spw, (baselines[j], baselines[j]), 'pI'))
        for j in range(len(baselines)):
            powerg_bl.append((np.real(uvpg.get_data(key_bl[j]))))
        # delays
        dlys = uvpg.get_dlys(spw) * 1e9
        
        # P(k) with Gaussian noise
        Powers_g.append(powerg_bl)
        
        # Mean power over all delays and time samples
        Powers_g_dm = np.mean(np.mean(Powers_g, axis=2), axis=2)
        
    return Powers_g, Powers_g_dm, dlys

# Get arrays of mean powers and delays for different sigmas.
# Set the number of loops (noise realizations) n here

n_l = 500
# sigma = 0.1, mu = 0
P_01 = Powerg_m(n_l, 0, 0.1)
# sigma = 0.5, mu = 0
P_05 = Powerg_m(n_l, 0, 0.5)
# sigma = 0.25, mu = 0
P_025 = Powerg_m(n_l, 0, 0.25)
# sigma = 0.1, mu = 0
P_1 = Powerg_m(n_l, 0, 0.1)

# Export the mean power arrays 
hi = np.savez("Parr_500loops_f.npz", power_1 = P_1[1], power_05 = P_05[1], power_025 = P_025[1], power_01 = P_01[1])
