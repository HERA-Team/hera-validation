import numpy as np
from hera_sim import noise, Simulator
from scipy import interpolate


def add_noise(sim, seed=1010):
    """
    Add thermal noise to a UVData object, based on its auto-correlations.

    Parameters
    ----------
    sim : :class:`pyuvdata.UVData` or :class:`hera_sim.Simulator` instance
        The simulation data to add noise to
    seed : int, optional
        The random seed.
    """
    if not isinstance(sim, Simulator):
        sim = Simulator(data=sim)

    lsts = np.unique(sim.data.lst_array)
    freqs = np.unique(sim.data.freq_array)

    # Set up to use the autos to set the noise level
    freqs_GHz = freqs / 1e9  # GHz
    omega_p = noise.bm_poly_to_omega_p(freqs_GHz)
    Jy_to_K = noise.jy2T(freqs_GHz, omega_p) / 1000
    autos = sim.data.get_data((0, 0, 'xx')) * Jy_to_K[None, :]
    autos_interp = interpolate.RectBivariateSpline(lsts, freqs_GHz, autos.real)

    np.random.seed(seed)
    noise = sim.add_noise('thermal_noise', Tsky_mdl=autos_interp, ret_vis=True)
    return sim, noise