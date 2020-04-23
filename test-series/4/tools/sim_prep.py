import argparse
import copy
import os
import re
import time
import yaml

import numpy as np
from scipy import interpolate, stats
from warnings import warn

from .data import DATA_PATH
from pyuvdata import UVData
from pyuvdata.utils import polstr2num
from hera_cal.io import to_HERAData
from hera_cal.abscal import rephase_vis, get_d2m_time_map
from hera_sim import Simulator

import hera_sim
if hera_sim.version.version.startswith('0'):
    from hera_sim.rfi import _listify
else:
    from hera_sim.utils import _listify

# ------- Functions for adding systematics ------- #

def add_noise(sim, Trx=100, seed=None, ret_cmp=True):
    """
    Add thermal noise to a UVData object, based on its auto-correlations.

    Parameters
    ----------
    sim : :class:`pyuvdata.UVData` or :class:`hera_sim.Simulator` instance
        The simulation data to add noise to
    Trx : float, optional
        The receiver temperature, in Kelvin.
    seed : int, optional
        The random seed. Not used if not specified. Use value 'random' 
        to have a seed automatically generated.
    """
    if not isinstance(sim, Simulator):
        sim = Simulator(data=sim)

    lsts = np.unique(sim.data.lst_array)
    freqs = np.unique(sim.data.freq_array)
    seed = _gen_seed(seed)

    # Find out which antenna has the autocorrelation data, in case we
    # add noise before inflating.
    antpair = next(ants for ants in sim.data.get_antpairs() if ants[0] == ants[1])

    # Set up to use the autos to set the noise level
    freqs_GHz = freqs / 1e9  # GHz
    beam_poly = np.load(os.path.join(DATA_PATH, "RIMEz_beam_poly.npy"))
    omega_p = np.polyval(beam_poly, freqs_GHz)
    Jy_to_K = hera_sim.noise.jy2T(freqs_GHz, omega_p) / 1000
    xx_autos = sim.data.get_data(*antpair, 'xx') * Jy_to_K[None, :]
    xx_autos_interp = interpolate.RectBivariateSpline(lsts, freqs_GHz, xx_autos.real)
    yy_autos = sim.data.get_data(*antpair, 'yy') * Jy_to_K[None, :]
    yy_autos_interp = interpolate.RectBivariateSpline(lsts, freqs_GHz, yy_autos.real)

    # Generate noise for each polarization separately.
    blts, _, xx_pol_ind = sim.data._key2inds('xx')
    _, _, yy_pol_ind = sim.data._key2inds('yy')
    xx_slice = (blts, 0, slice(None), xx_pol_ind[0])
    yy_slice = (blts, 0, slice(None), yy_pol_ind[0])
    if seed is not None:
        np.random.seed(seed)
    if ret_cmp:
        noise = np.zeros_like(sim.data.data_array, dtype=np.complex)
        noise[xx_slice] += sim.add_noise(
            'thermal_noise', Tsky_mdl=xx_autos_interp, Trx=Trx, omega_p=omega_p,
            ret_vis=True, add_vis=False
        )[xx_slice]
        noise[yy_slice] += sim.add_noise(
            'thermal_noise', Tsky_mdl=yy_autos_interp, Trx=Trx, omega_p=omega_p,
            ret_vis=True, add_vis=False
        )[yy_slice]
        sim.data.data_array += noise
        return sim, noise
    else:
        sim.data.data_array[xx_slice] += sim.add_noise(
            'thermal_noise', Tsky_mdl=xx_autos_interp, Trx=Trx, omega_p=omega_p,
            ret_vis=True, add_vis=False
        )[xx_slice]
        sim.data.data_array[yy_slice] += sim.add_noise(
            'thermal_noise', Tsky_mdl=yy_autos_interp, Trx=Trx, omega_p=omega_p,
            ret_vis=True, add_vis=False
        )[yy_slice]
        return sim


def add_gains(
    sim, 
    seed=None, 
    time_vary_params=None,
    ret_cmp=True,
    **gain_params
):
    """
    Add per-antenna bandpass gains to the simulation.

    Parameters
    ----------
    sim : :class:`hera_sim.Simulator` or :class:`pyuvdata.UVData`
        The object containing the simulation data and metadata.
    seed : int, optional
        The random seed. Not used if not specified. Use value 'random' 
        to have a seed automatically generated.
    time_vary_params : dict, optional
        Parameters for adding time variation to the gains. Keys should 
        be any of ('amp', 'amplitude', 'phs', 'phase'). Values should 
        be the set of variation parameters (parameters prefixed by 
        'variation' in the ``vary_gains_in_time`` function) to be used.
        Default behavior is to leave the gains constant in time.
        See ``vary_gains_in_time`` function for details on parameters.
    **gain_params
        The parameters for simulating bandpass gains.

    Returns
    -------
    corrupted_sim : :class:`pyuvdata.UVData`
        Simulation with bandpass gains applied. 
    gains : dict
        Dictionary mapping antenna numbers to bandpass gains as a 
        function of frequency (and potentially time).
    """
    # Setup
    sim = _sim_to_uvd(sim)
    freqs_GHz = np.unique(sim.freq_array) / 1e9
    ants = sim.antenna_numbers
    seed = _gen_seed(seed)

    # Simulate and apply the gains
    if seed is not None:
        np.random.seed(seed)
    gains = hera_sim.sigchain.gen_gains(freqs_GHz, ants, **gain_params)
    apply_gains(sim, gains, time_vary_params)

    return sim if not ret_cmp else sim, gains

def add_reflections(
    sim, 
    seed=None, 
    dly=1200, 
    dly_spread=0, 
    amp=1e-3, 
    amp_scale=0, 
    time_vary_params=None,
    ret_cmp=True
):
    """
    Add per-antenna reflection gains to the simulation.

    Parameters
    ----------
    sim : :class:`hera_sim.Simulator` or :class:`pyuvdata.UVData`
        The object containing the simulation data and metadata.
    seed : int, optional
        The random seed. Not used if not specified. Use value 'random' 
        to have a seed automatically generated.
    dly : float, optional
        Delay at which the reflection appears, in nanoseconds. 
        Default is 1200 ns.
    dly_spread : float, optional
        Absolute amount the delays vary between antennas, in ns. 
        Default is 0 ns.
    amp : float, optional
        Amplitude of the reflection. Default is 1e-3.
    amp_scale : float, optional
        Fractional variation in the reflection amplitude between antennas.
        Default is 0. 
    time_vary_params : dict, optional
        Parameters for adding time variation to the gains. Keys should 
        be any of ('amp', 'amplitude', 'phs', 'phase'). Values should 
        be the set of variation parameters (parameters prefixed by 
        'variation' in the ``vary_gains_in_time`` function) to be used.
        Default behavior is to leave the gains constant in time.
        See ``vary_gains_in_time`` function for details on parameters.

    Returns
    -------
    corrupted_sim : :class:`pyuvdata.UVData`
        Simulation with reflection gains applied. 
    gains : dict
        Dictionary mapping antenna numbers to reflection gains as a 
        function of frequency (and potentially time).
    """
    # Setup
    sim = _sim_to_uvd(sim)
    freqs_GHz = np.unique(sim.freq_array) / 1e9
    ants = sim.antenna_numbers
    Nants = len(ants)
    seed = _gen_seed(seed)

    # Randomize parameters in a realistic way.
    if seed is not None:
        np.random.seed(seed)
    delays_ns = stats.norm.rvs(dly, dly_spread, Nants)
    phases = stats.uniform.rvs(0, 2*np.pi, Nants)
    amps = amp * stats.norm.rvs(1, amp_scale, Nants)

    # Simulate and apply the gains.
    gains = hera_sim.sigchain.gen_reflection_gains(
        freqs_GHz, ants, amp=amps, dly=delays_ns, phs=phases
    )
    apply_gains(sim, gains, time_vary_params)

    return sim if not ret_cmp else sim, gains

def add_xtalk(
    sim, 
    seed=None, 
    Ncopies=10, 
    amp_range=(-4,-6), 
    dly_rng=(900,1300),
    ret_cmp=True
):
    """
    Add cross-coupling crosstalk to the simulation's cross-correlations.

    Parameters
    ----------
    sim : :class:`hera_sim.Simulator` or :class:`pyuvdata.UVData`
        The object containing the simulation data and metadata.
    seed : int, optional
        The random seed. Not used if not specified. Use value 'random' 
        to have a seed automatically generated.
    Ncopies : int, optional
        Number of delays at which to generate crosstalk. Default is 
        to use 10 copies.
    amp_range : tuple of float, optional
        Base-10 logarithm of the amplitude of the crosstalk at the 
        first and last delays where crosstalk is injected. Default 
        is to use -4 and -6.
    dly_rng : tuple of float, optional
        Minmum and maximum delays, in ns, at which to inject 
        crosstalk. Crosstalk is inserted at ``Ncopies`` delays spaced 
        linearly between the minimum and maximum delays. 
        Default is 900-1300 ns.

    Returns
    -------
    corrupted_sim : :class:`pyuvdata.UVData`
        Simulation with crosstalk applied.
    xtalk : np.ndarray
        Data for the crosstalk visibilities.
    """
    # Setup
    sim = _sim_to_uvd(sim)
    freqs_GHz = np.unique(sim.freq_array) / 1e9
    amps = np.logspace(*amp_range, Ncopies)
    dlys = np.linspace(*dly_rng, Ncopies)
    seed = _gen_seed(seed)

    if ret_cmp:
        xtalk = np.zeros_like(sim.data_array, dtype=np.complex)
    if seed is not None:
        np.random.seed(seed)
    for antpairpol in sim.get_antpairpols():
        ai, aj, pol = antpairpol
        if ai == aj:
            continue
        blt_inds, _, pol_inds = sim._key2inds(antpairpol)
        this_slice = (blt_inds, 0, slice(None), pol_inds[0])

        # Add some jitter to the delays and amplitudes.
        ddlys = stats.norm.rvs(1, 1e-2, Ncopies)
        damps = stats.norm.rvs(1, 1e-4, Ncopies)

        # Generate random phases and pull the autocorrelation.
        phs = stats.uniform.rvs(0, 2*np.pi, Ncopies)
        autovis = sim.get_data(ai, ai, pol)
        if ret_cmp:
            xtalk[this_slice] = gen_xtalk(
                autovis, freqs_GHz, amps * damps, dlys + ddlys, phs
            )
        else:
            sim.data_array[this_slice] = gen_xtalk(
                autovis, freqs_GHz, amps * damps, dlys + ddlys, phs
            )

    if ret_cmp:
        sim.data_array += xtalk
        return sim, xtalk
    else:
        return sim

def apply_gains(sim, gains, time_vary_params=None):
    """Apply per-antenna gains to a simulation."""
    # Support gain variation in time.
    times = np.unique(sim.time_array)
    time_vary_params = time_vary_params or {}
    for mode, vary_params in time_vary_params.items():
        gains = vary_gains_in_time(gains, times, mode, **vary_params)
    for antpairpol in sim.get_antpairpols():
        blt_inds, _, pol_inds = sim._key2inds(antpairpol)
        this_slice = (blt_inds, 0, slice(None), pol_inds[0])
        vis = sim.get_data(antpairpol)
        sim.data_array[this_slice] = hera_sim.sigchain.apply_gains(
            vis, gains, antpairpol[:2]
        )
    return

def vary_gains_in_time(
    gains, 
    times, 
    mode='amp',
    variation_ref_times=None, 
    variation_timescales=None,
    variation_amps=(0.05,),
    variation_modes=('linear',),
    clip_gains=True
):
    """Vary gain amplitudes or phases in time.
    
    Parameters
    ----------
    gains : dict
        Dictionary containing the per-antenna gains to make vary in time.
    times : array-like of float
        Times at which to simulate time variation. Should match the unique 
        times in the visibility metadata.
    mode : str, optional
        Whether to vary the amplitudes or the phases. Accepts 'amp', 
        'amplitude', 'phs', or 'phase' as a value. Default is 'amp'.
    variation_ref_times : float or array-like of float, optional
        Reference time(s) to use for generating time variation. If not 
        using the 'noiselike' option, then this is the time where the 
        gains are unchanged from their original values. Default is to 
        use the middle time. This should have the same units as the 
        `times` parameter.
    variation_timescales : float or array-like of float, optional
        Timescale on which the variation is supposed to take place. For 
        sinusoidal variation, this is the period of the variation. This 
        should be in the same units as the `times` parameter. Default 
        is to use the entire length of `times`.
    variation_amps : float or array-like of float, optional
        Amplitude of the variation, as a fraction of the input gains. 
        This is *not* the peak-to-peak amplitude! Default is 5%.
    variation_modes : str or array-like of str, optional
        Which type(s) of time variation to simulate. Currently supported 
        options are 'linear', 'sinusoidal', and 'noiselike'. Default is 
        to vary gains linearly over the given times.
    clip_gains : bool, optional
        Whether to clip the gains so that they do not vary by more than 
        100% of their original values. Default is True. Raises a warning 
        if gains are clipped.

    Returns
    -------
    time_varied_gains : dict
        Dictionary of gains with time variation applied.
    """
    # Setup for handling multiple modes of variation at once.
    if variation_ref_times is None:
        variation_ref_times = (np.median(times),)
    if variation_timescales is None:
        variation_timescales = (times[-1] - times[0],)
    variation_modes = _listify(variation_modes)
    variation_amps = _listify(variation_amps)
    variation_ref_times = _listify(variation_ref_times)
    variation_timescales = _listify(variation_timescales)
    variation_params = (
        variation_modes, variation_amps, variation_ref_times, variation_timescales
    )
    Nmodes = max(len(param) for param in variation_params)
    if any(Nmodes % param != 0 for param in variation_params):
        raise ValueError(
            "Input variation parameters are coprime in length. "
            "Please ensure variation parameters are either all "
            "the same length, or can be broadcast to the same "
            "length as the longest variation parameter."
        )
    iterator = zip(
        variation_ref_times * (Nmodes // len(variation_ref_times)),
        variation_timescales * (Nmodes // len(variation_timescales)),
        variation_modes * (Nmodes // len(variation_modes)),
        variation_amps * (Nmodes // len(variation_amps)),
    )

    # Generate an envelope for adding time variation.
    envelope = 1
    for ref_time, timescale, variation_amp, variation_mode in iterator:
        phases = (times - ref_time) / timescale
        if variation_mode == 'linear':
            envelope *= 1 + 2 * vary_amp * phases
        elif variation_mode == 'sinusoidal':
            envelope *= 1 + vary_amp * np.sin(2 * np.pi * phases)
        elif variation_mode == 'noiselike':
            envelope *= stats.norm.rvs(1, vary_amp, times.size)
    if not np.all(envelope == 1):
        warn("Time variation method not supported; returning.")
        return gains
    if clip_gains and not np.allclose(envelope, np.clip(envelope, 0, 2)):
        warn(
            "Gains vary in time by an amount more than their time-"
            "independent values. Clipping time-variation envelope."
            "Check your parameters if this is unexpected."
        )
        envelope = np.clip(envelope, 0, 2)
    envelope = np.outer(envelope, np.ones(gain_shape[-1]))
    if mode in ('phs', 'phase'):
        envelope = np.exp(1j*envelope)

    # Actually apply the time variation.
    gain_shape = list(gains.values())[0].shape
    if len(gain_shape) == 1:
        gains = {ant : gain[None,:] * envelope for ant, gain in gains.items()}
    else:
        gains = {ant : gain * envelope for ant, gain in gains.items()}
    return gains

def gen_xtalk(autovis, freqs, xamps, xdlys, xphs):
    """Generate a series of cross-coupling crosstalk visibilities."""
    xtalk = np.zeros_like(autovis, dtype=np.complex)
    _gen_xtalk = hera_sim.sigchain.gen_cross_coupling_xtalk
    for amp, dly, phs in zip(xamps, xdlys, xphs):
        xtalk += _gen_xtalk(freqs, autovis, amp, dly, phs)
        xtalk += _gen_xtalk(freqs, autovis, amp, -dly, phs)
    return xtalk

SYSTEMATICS_SIMULATORS = {
    'noise' : add_noise,
    'gains' : add_gains,
    'reflections' : add_reflections,
    'xtalk' : add_xtalk
}

def apply_systematics(
    sim, 
    seed=None,
    noise=None, 
    gains=None, 
    reflections=None,
    xtalk=None,
    return_systematics=True
):
    """One-stop shop for applying systematics to a simulation.

    # TODO: docstring
    """
    # If sim components are specified, initialize parameter dictionaries.
    # This allows users to pass an empty dictionary and use the same 
    # value for ``seed`` throughout. This is useful in case the user 
    # wants to have the random seeds automatically generated by using 
    # the 'random' setting for ``seed``, but all other settings default.
    if noise is not None:
        noise = noise or {'seed' : seed}
    if gains is not None:
        gains = gains or {'seed' : seed}
    if reflections is not None:
        reflections = reflections or {'seed' : seed}
    if xtalk is not None:
        xtalk = xtalk or {'seed' : seed}
    
    # Collect all of the parameters into a dictionary.
    parameters = {
        'noise' : noise,
        'gains' : gains,
        'reflections' : reflections,
        'xtalk' : xtalk
    }
    # Update random seeds in case user wants to "randomly" generate seeds.
    for params in parameters.values():
        if params is not None and 'seed' in params.keys():
            params['seed'] = _gen_seed(params['seed'])

    # Apply the systematics and track the results.
    if return_systematics:
        systematics = {}
    sim = _sim_to_uvd(sim)
    for systematic, params in parameters.items():
        # This is a bit of a hack, but I can't think of a better way...
        add_systematic = SYSTEMATICS_SIMULATORS[systematic]
        if return_systematics:
            sim, systematics[systematic] = add_systematic(sim, **params)
        else:
            sim = add_systematic(sim, ret_cmp=False, **params)

    # Simulating this and keeping all the systematics may be very 
    # memory-intensive, so sometimes we might not want to keep 
    # track of the intermediate products.
    if not return_systematics:
        systematics = None

    return sim, systematics, parameters

# ------- Functions for preparing files ------- #

def adjust_sim_to_data(sim_file, data_files, verbose=False):
    """
    Modify simulated data to be consistent with an observation's metadata.
    
    Parameters
    ----------
    sim_file : str or path-like object
        Path to simulation file to modify.
    data_files : array-like of str or path-like objects
        Collection of paths to observed data.

    Returns
    -------
    modified_sim_data : :class:`pyuvdata.UVData`
        Simulation data that has been modified to match the metadata 
        contents of ``data_files``. See notes for details.
        
    Notes
    -----
    The modified simulation data has its antennas relabeled and repositioned to match 
    the labels and positions of the observational data. Any antenna from the 
    simulation that does not correspond to an antenna used in observation is excluded 
    from the modified simulation data. See ``downselect_antennas`` to learn more 
    about how the antenna-related metadata is modified. Note that this process 
    requires that the data be fully inflated; in the interest of keeping memory 
    consumption low, this step is performed after LST rephasing.
    
    It is also worth noting that the modified simulation data has not had its flag or 
    nsamples arrays modified in any way, but its data array has been modified. The 
    modified data has been downselected in time, rephased to have its LSTs match the 
    observed LSTs, and its time and LST arrays have been updated to match the observed 
    data.
    
    The modified simulation data is chunked into files containing the same number of 
    integrations as the observation files (assuming a uniform number of integration 
    per file) and written to disk. 
    """
    # Ensure the data files are a list and load in their metadata.
    data_files = _listify(data_files)
    ref_uvd = UVData()
    ref_uvd.read(data_files, read_data=False)
    
    # Load in the simulation data, but only the linear vispols.
    use_pols = [polstr2num(pol) for pol in ('xx', 'yy')]
    sim_uvd = UVData()
    sim_uvd.read(sim_file, polarizations=use_pols)

    # Downselect in time and rephase LSTs to match the reference data.
    if verbose:
        print("Rephasing simulation to reference data...")
    sim_uvd = rephase_to_reference(sim_uvd, ref_uvd)

    # Inflate the simulation so antenna downselection can actually be done.
    if verbose:
        print("Inflating simulation by redundancy...")
    sim_uvd.inflate_by_redundancy()
    
    # Find and use the intersection of the RIMEz and H1C arrays.
    if verbose:
        print("Choosing subset of antennas to keep...")
    sim_uvd = downselect_antennas(sim_uvd, ref_uvd)
    
    # Make sure the data is conjugated properly so redcal doesn't break.
    if verbose:
        print("Conjugating simulation to ant1<ant2 convention...")
    sim_uvd.conjugate_bls('ant1<ant2')
    
    return sim_uvd

def prepare_sim_files(
    sim_file, 
    data_files, 
    save_dir, 
    sky_cmp=None,
    systematics_params=None, 
    save_truth=True,
    clobber=True,
    verbose=True
):
    """
    Modify simulation data to match reference metadata and apply systematics.
    
    Parameters
    ----------
    sim_file : str or path-like object
        Path to simulation file to modify.
    data_files : array-like of str or path-like objects
        Collection of paths to observed data.
    save_dir : str or path-like object
        Path to directory where modified simulation files should be saved.
    sky_cmp : str, optional
        Sky component simulated in ``sim_file``. If not provided, then the component 
        is inferred from the file name.
    systematics_params : dict, optional
        Dictionary mapping systematics names (i.e. 'noise', 'gains') to 
        a dictionary of parameters to use for generating the systematic.
        Default is to not simulate/apply systematics.
    save_truth : bool, optional
        Whether to save the true visibilities. Default is True.
    clobber : bool, optional
        Whether to overwrite files that may share the same name as the new files to 
        be saved in ``save_dir``. Default is to overwrite files.
    """
    if verbose:
        print("Beginning file preparation...")

    # Get the sky component if not specified.
    sky_cmp = sky_cmp or _parse_filename_for_cmp(sim_file)

    # Don't do anything if the files already exist and clobber is False.
    true_files_exist = _sim_files_exist(data_files, save_dir, sky_cmp, 'true')
    corrupt_files_exist = _sim_files_exist(data_files, save_dir, sky_cmp, 'corrupt')
    save_corrupt = bool(systematics_params)
    write_truth = save_truth and (clobber or not true_files_exist)
    write_corrupt = save_corrupt and (clobber or not corrupt_files_exist)
    if not (write_truth or write_corrupt):
        print("No new files to make; returning.")
        return

    # Modify the simulation data to match the reference metadata.
    sim_uvd = adjust_sim_to_data(sim_file, data_files, verbose=verbose)

    # Chunk the simulation data and write to disk.
    if write_truth:
        if verbose:
            print("Writing true visibilities to disk...")
        chunk_sim_and_save(
            sim_uvd, data_files, save_dir, 
            sky_cmp=sky_cmp, state='true', clobber=clobber
        )

    # Apply systematics if desired.
    if write_corrupt:
        # TODO: update this to handle being able to save the systematics
        # but maybe do a warning if the task may cause a MemoryError
        if verbose:
            print("Simulating and applying systematics...")
        sim_uvd, systematics, params = apply_systematics(
            sim_uvd, return_systematics=False, **systematics_params
        )
        if verbose:
            print("Writing corrupted visibilities to disk...")
        chunk_sim_and_save(
            sim_uvd, data_files, save_dir, 
            sky_cmp=sky_cmp, state='corrupt', clobber=clobber
        )
        if verbose:
            print("Writing sysytematics parameters to disk...")
        save_config(params, data_files[0], save_dir, sky_cmp, clobber)
    return

def downselect_antennas(sim_uvd, ref_uvd, tol=1.0):
    """
    Downselect simulated array to match unique baselines from reference.
    
    Parameters
    ----------
    sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the simulation data.
    ref_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the reference data. 
        Only needs to have the metadata loaded.
    tol : float, optional
        The maximum allowable discrepancy in antenna position components.
        Default is 1 meter.
        
    Returns
    -------
    new_sim_uvd : :class:`pyuvdata.UVData`
        Simulation data adjusted to contain a subset of the simulated array
        such that, after appropriately translating the antenna positions, 
        the maximum number of antennas/unique baselines remain in the 
        intersection of the simulation subarray and the reference array. 
        The relevant attributes of the simulation data are updated 
        accordingly.

    Notes
    -----
    See the memo "Antennas, Baselines, and Maximum Overlap" for details on 
    the antenna selection process. It is worth noting that the algorithm 
    implemented here is *not* appropriate for general arrays, but the 
    general case should reduce to this case for an array with appropriate 
    symmetry, as we have for the RIMEz simulation.
    """
    sim_uvd = _sim_to_uvd(sim_uvd)

    # Find the optimal intersection, get the map between antennas, and downselect the data.
    sim_ENU_antpos = _get_antpos(sim_uvd)
    ref_ENU_antpos = _get_antpos(ref_uvd)
    array_intersection = _get_array_intersection(sim_ENU_antpos, ref_ENU_antpos, tol)
    sim_to_ref_ant_map = _get_antenna_map(array_intersection, ref_ENU_antpos, tol)
    sim_uvd.select(antenna_nums=list(sim_to_ref_ant_map.keys()), keep_all_metadata=False)
    ref_uvd.select(antenna_nums=list(sim_to_ref_ant_map.values()), keep_all_metadata=False)
    
    # Prepare the new data array.
    new_sim_data = np.zeros_like(sim_uvd.data_array, dtype=np.complex)
    new_sim_flags = np.zeros_like(sim_uvd.flag_array, dtype=np.bool)
    new_sim_nsamples = np.zeros_like(sim_uvd.nsample_array, dtype=np.float)
    for antpairpol in sim_uvd.get_antpairpols():
        ai, aj, pol = antpairpol
        ref_antpairpol = (sim_to_ref_ant_map[ai], sim_to_ref_ant_map[aj], pol)
        blts, conj_blts, pol_inds = ref_uvd._key2inds(ref_antpairpol)
        sim_data = sim_uvd.get_data(antpairpol)
        # Correctly choose which slice to use depending on whether the 
        # reference baseline corresponding to (ai, aj) is conjugated.
        # (If it's conjugated in the reference, then blts is empty.)
        if len(blts) > 0:
            this_slice = (blts, 0, slice(None), pol_inds[0])
        else:
            this_slice = (conj_blts, 0, slice(None), pol_inds[1])
            sim_data = sim_data.conj()
        new_sim_data[this_slice] = sim_data
        new_sim_flags[this_slice] = sim_uvd.get_flags(antpairpol)
        new_sim_nsamples[this_slice] = sim_uvd.get_nsamples(antpairpol)
        
    # Update the data-like parameters.
    sim_uvd.data_array = new_sim_data
    sim_uvd.flag_array = new_sim_flags
    sim_uvd.nsample_array = new_sim_nsamples

    # Update the antenna-related simulation metadata to reflect the changes 
    # made to the data array; this ensures that the data-like parameters
    # will be accessed correctly.
    sim_uvd.ant_1_array = ref_uvd.ant_1_array
    sim_uvd.ant_2_array = ref_uvd.ant_2_array
    sim_uvd.antenna_numbers = ref_uvd.antenna_numbers
    sim_uvd.antenna_positions = ref_uvd.antenna_positions
    sim_uvd.telescope_location = ref_uvd.telescope_location
    sim_uvd.history += "\nAntennas adjusted to optimally match H1C antennas."
    
    return sim_uvd

def rephase_to_reference(sim_uvd, ref_uvd):
    """
    Rephase simulation LSTs to reference LSTs after downselection in time.
    
    After downselection and rephasing, this function overwrites the 
    simulation LST and time arrays to match the reference LSTs and times.
    
    Parameters
    ----------
    sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the simulation data.
    ref_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the reference data. 
        Only the metadata for this object needs to be read.

    Returns
    -------
    new_sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the modified 
        simulation data. The modified data has been rephased to match 
        the reference LSTs and has had its time and LST arrays adjusted 
        to contain only the times and LSTs that are present in the 
        reference data.
    """
    # Convert the simulation to a HERAData object. 
    sim_uvd = to_HERAData(copy.deepcopy(sim_uvd))
    hd_metas = sim_uvd.get_metadata_dict()

    # Load in useful metadata.
    sim_lsts = hd_metas['lsts']
    sim_times = hd_metas['times']
    ref_lsts = np.unique(ref_uvd.lst_array)
    ref_times = np.unique(ref_uvd.time_array)
    sim_to_ref_time_map = get_d2m_time_map(sim_times, sim_lsts, ref_times, ref_lsts)
    use_times = np.array([
        sim_time for sim_time, ref_time in sim_to_ref_time_map.items()
        if ref_time is not None
    ])
    use_lsts = np.array(
        [sim_lsts[np.argmin(np.abs(sim_times - time))] for time in use_times]
    )
    ref_lsts = np.array([
        ref_lsts[np.argmin(np.abs(ref_times - sim_to_ref_time_map[sim_time]))] 
        for sim_time in use_times
    ])
    ref_times = np.array([sim_to_ref_time_map[sim_time] for sim_time in use_times])

    # Downselect in time and load data.
    #ref_uvd.select(times=ref_times) # Figure out why the times don't match exactly.
    sim_uvd.select(times=use_times) # We should only need to do this one select.
    data, _, _ = sim_uvd.build_datacontainers()

    # Build the antpair -> ENU baseline vector dictionary.
    antpos, ants = sim_uvd.get_ENU_antpos()
    bls = {bl : None for bl in data.bls()}
    for bl in bls:
        ai, aj = bl[:2]
        i, j = ants.tolist().index(ai), ants.tolist().index(aj)
        bls[bl] = antpos[j] - antpos[i]

    # Rephase and update the data.
    new_data, new_flags = rephase_vis(
        data, use_lsts, ref_lsts, bls, hd_metas['freqs']
    )
    sim_uvd.update(data=new_data, flags=new_flags)
    new_sim_lsts = np.zeros_like(sim_uvd.lst_array)
    new_sim_times = np.zeros_like(sim_uvd.time_array)
    for use_time, ref_time, ref_lst in zip(use_times, ref_times, ref_lsts):
        blt_slice = np.argwhere(sim_uvd.time_array == use_time)
        new_sim_lsts[blt_slice] = ref_lst
        new_sim_times[blt_slice] = ref_time
    sim_uvd.time_array = new_sim_times
    sim_uvd.lst_array = new_sim_lsts

    # Ensure that we return a UVData object so that chunk_sim_and_save doesn't break.
    return_uvd = UVData()
    for _property in return_uvd:
        setattr(return_uvd, _property, getattr(sim_uvd, _property))
    return return_uvd

def chunk_sim_and_save(
    sim_uvd, 
    save_dir, 
    ref_files=None,
    Nint_per_file=None,
    sky_cmp=None, 
    state=None, 
    clobber=True
):
    """
    Chunk the simulation data to match the reference file and write to disk.
    
    Chunked files have the following naming convention:
    save_dir/zen.{jd_major}.{jd_minor}[.{sky_cmp}][.{state}].uvh5
    The entires in brackets are optional and may be omitted.

    Parameters
    ----------
    sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the simulation data 
        to chunk and write to disk.
    save_dir : str or path-like object
        Path to the directory where the chunked files will be saved.
    ref_files : iterable of str
        Iterable of filepaths to use for reference when chunking. This must 
        be specified if ``Nint_per_file`` is not specified. This determines 
        (and overrides, if also provided) ``Nint_per_file`` if provided.
    Nint_per_file : int, optional
        Number of integrations per chunked file. This must be specified 
        if ``ref_files`` is not specified.
    sky_cmp : str, optional
        String denoting which sky component has been simulated. Should 
        be one of the following: ('foregrounds', 'eor', 'sum').
    state : str, optional
        String denoting whether the file is the true sky or corrupted.
    clobber : bool, optional
        Whether to overwrite any existing files that share the new 
        filenames. Default is to overwrite files.
    """
    if ref_files is None and Nint_per_file is None:
        raise ValueError(
            "Either a glob of reference files or the number of integrations "
            "per file must be provided."
        )

    # Pull the number of integrations per file if needed.
    if ref_files is not None:
        uvd = UVData()
        uvd.read(ref_files[0], read_data=False)
        Nint_per_file = uvd.Ntimes
        jd_pattern = re.compile(r"\.(?P<major>[0-9]{7})\.(?P<minor>[0-9]{5}).")

    # Pull the simulation times, then start the chunking process.
    sim_times = np.unique(sim_uvd.time_array)
    Nfiles = int(sim_uvd.Ntimes / Nint_per_file)
    for Nfile in range(Nfiles):
        # Figure out filing and slicing information.
        if ref_files is not None:
            jd = re.search(jd_pattern, ref_files[Nfile]).groupdict()
            jd = float(f"{jd['major']}.{jd['minor']}")
            start_ind = np.argmin(np.abs(sim_times - jd))
        else:
            start_ind = Nfile * Nint_per_file
            jd = np.round(sim_times[start_ind], 5)
        this_slice = slice(start_ind, start_ind + Nint_per_file)
        filename = f"zen.{jd}.uvh5"
        if sky_cmp is not None:
            filename = filename.replace(".uvh5", f".{sky_cmp}.uvh5")
        if state is not None:
            filename = filename.replace(".uvh5", f".{state}.uvh5")
        save_path = os.path.join(save_dir, filename)

        # Chunk it and write to disk.
        times = sim_times[this_slice]
        this_uvd = sim_uvd.select(times=times, inplace=False)
        this_uvd.write_uvh5(save_path, clobber=clobber)
    return

def save_config(config, ref_file, save_dir, sky_cmp, clobber=True):
    """
    Write the configuration used for systematics simulation to disk.
    
    Parameters
    ----------
    config : dict
        Dictionary mapping names of systematics components (e.g. 'noise' 
        or 'gains') to parameters used to simulate the systematic.
    ref_file : str
        Path to a reference file, to be used for extracting the JD.
    save_dir : str or path-like object
        Path to the directory where the chunked files will be saved.
    sky_cmp : str
        String denoting which sky component has been simulated. Should 
        be one of the following: ('foregrounds', 'eor', 'sum').
    clobber : bool, optional
        Whether to overwrite any existing files that share the new 
        filenames. Default is to overwrite files.
    """
    jd_pattern = re.compile("[0-9]{7}")
    jd = jd_pattern.findall(ref_file)[0]
    filename = f"{jd}.config.{sky_cmp}.yaml"
    save_path = os.path.join(save_dir, filename)
    if not clobber and os.path.exists(save_path):
        print("Config file exists and clobber set to False; returning.")
        return
    with open(save_path, 'w') as f:
        f.write(yaml.dump(config))
    return

# ------- Helper Functions ------- #

def _sim_to_uvd(sim):
    """Update simulation object type."""
    if isinstance(sim, hera_sim.Simulator):
        sim = sim.data
    elif isinstance(sim, str):
        sim_ = UVData()
        sim_.read(sim)
        sim = sim_
    return sim

def _gen_seed(seed):
    """Generate a random seed pseudo-randomly if desired."""
    if seed is None or type(seed) is int:
        return seed
    elif seed == 'random':
        return time.time_ns() % 2**32

def _get_array_intersection(sim_antpos, ref_antpos, tol=1.0):
    """Find the optimal choice of simulation subarray and return it."""
    optimal_translation = _get_optimal_translation(sim_antpos, ref_antpos, tol)
    new_antpos = {ant : pos + optimal_translation for ant, pos in sim_antpos.items()}
    intersection = {
        ant : pos for ant, pos in new_antpos.items()
        if any(np.allclose(pos, ref_pos, atol=tol) for ref_pos in ref_antpos.values())
    }
    return intersection

def _get_antenna_map(sim_antpos, ref_antpos, tol=1.0):
    """Find a mapping from simulation antennas to reference antennas."""
    antenna_map = {}
    refants = list(ref_antpos.keys())
    for ant, pos in sim_antpos.items():
        refant_index = np.argwhere(
            [np.allclose(pos, ref_pos, atol=tol) for ref_pos in ref_antpos.values()]
        ).flatten()
        if refant_index.size == 0:
            continue
        antenna_map[ant] = refants[refant_index[0]]
    return antenna_map

def _get_antpos(uvd, ENU=True):
    """Retrieve the {ant : pos} dictionary from the data."""
    if ENU:
        pos, ant = uvd.get_ENU_antpos()
    else:
        ant = uvd.antenna_numbers
        pos = uvd.antenna_positions
    return dict(zip(ant, pos))

def _get_optimal_translation(sim_antpos, ref_antpos, tol=1.0):
    """Find the translation that maximizes the overlap between antenna arrays."""
    # Get a dictionary of translations; keys are sim_ant -> ref_ant.
    translations = _build_translations(sim_antpos, ref_antpos, tol)
    intersection_sizes = {}
    
    # Calculate the number of antennas in the intersection for each pair of arrays.
    for sim_ant_to_ref_ant, translation in translations.items():
        new_sim_antpos = {
            ant : pos + translation for ant, pos in sim_antpos.items()
        }
        Nintersections = sum(
            any(np.allclose(pos, ref_pos, atol=tol) 
                for ref_pos in ref_antpos.values()
            )
            for pos in new_sim_antpos.values()
        )
        intersection_sizes[sim_ant_to_ref_ant] = Nintersections
        
    # Choose the translation that has the most antennas in the intersection.
    sim_to_ref_keys = list(translations.keys())
    intersections_per_translation = np.array(list(intersection_sizes.values()))
    index = np.argmax(intersections_per_translation)
    return translations[sim_to_ref_keys[index]]

def _build_translations(sim_antpos, ref_antpos, tol=1.0):
    """Build all possible translations that overlap at least one antenna."""
    sim_ants = list(sim_antpos.keys())
    ref_ants = list(ref_antpos.keys())
    translations = {
        f"{sim_ant}->{ref_ant}" : ref_antpos[ref_ant] - sim_antpos[sim_ant]
        for sim_ant in sim_ants for ref_ant in ref_ants
    }
    unique_translations = {}
    for key, translation in translations.items():
        if not any(
            np.allclose(translation, unique_translation, atol=tol)
            for unique_translation in unique_translations.values()
        ):
            unique_translations[key] = translation
    return unique_translations

def _parse_filename_for_cmp(filename):
    """Infer the sky component from the provided filename."""
    cmp_re = re.compile("[a-z]+.uvh5")
    sky_cmp = cmp_re.findall(str(filename))
    if sky_cmp == []:
        raise ValueError(
            f"Simulation component could not be inferred from {filename}."
        )
    return sky_cmp[0][:-5]

def _sim_files_exist(data_files, save_dir, sky_cmp, state):
    """Check if the chunked simulation files already exist."""
    file_ext = os.path.splitext(data_files[0])[1]
    data_file_names = [os.path.split(dfile)[1] for dfile in data_files]
    sim_file_names = [data_file.replace("HH", f"{sky_cmp}") for data_file in data_file_names]
    # Catch cases where the data files don't have 'HH' in their name.
    if all(
        sim_file == data_file for sim_file, data_file in zip(sim_file_names, data_file_names)
    ):
        sim_file_names = [
            sim_file.replace(file_ext, f".{sky_cmp}" + file_ext)
            for sim_file in sim_file_names
        ]
    sim_files = [
        os.path.join(save_dir, sim_file.replace(file_ext, f".{state}" + file_ext))
        for sim_file in sim_file_names
    ]
    return all(os.path.exists(sim_file) for sim_file in sim_files)

# ------- Argparsers ------- #

def sim_prep_argparser():
    """Argparser for preparing simulation files and adding systematics."""
    desc = "Modify simulation files to match observation parameters "
    desc += "and optionally apply systematics."
    a = argparse.ArgumentParser(description=desc)
    a.add_argument("--verbose", default=False, action="store_true")
    file_opts = a.add_argument_group(title="File preparation parameters.")
    file_opts.add_argument("simfile", type=str, help="Simulation file to be modified.")
    file_opts.add_argument("obsdir", type=str, help="Directory containing observation files.")
    file_opts.add_argument("savedir", type=str, help="Destination to write modified files.")
    file_opts.add_argument(
        "--skip_truth", default=False, action="store_true",
        help="Skip writing the true visibilities to disk; only write corrupted."
    )
    file_opts.add_argument(
        "--clobber", default=False, action="store_true", 
        help="Overwrite existing modified simulation files."
    )
    file_opts.add_argument(
        "--Nchunks", default=1, type=int, 
        help="Number of chunks to break obsfiles glob into for file preparation."
    )
    systematics = a.add_argument_group(title="Options for simulating systematics.")
    systematics.add_argument("--seed", default=None, help="Random seed for all components.")
    systematics.add_argument(
        "--config", default=None, type=str, help="Configuration file for systematics."
    )
    args = a.parse_args()
    return args

def abscal_model_argparser():
    """Argparser for preparing the abscal model."""
    desc = "Prepare an absolute calibration model from a simulation file."
    a = argparse.ArgumentParser(description=desc)
    a.add_argument("simfile", type=str, help="Simulation file to use.")
    a.add_argument("reffile", type=str, help="Reference for array selection.")
    a.add_argument("savedir", type=str, help="Where to write files to disk.")
    a.add_argument("Nint_per_file", type=int, help="Number of integrations per file.")
    a.add_argument("lst_min", type=float, help="Lower bound of LST range, in hours.")
    a.add_argument("lst_max", type=float, help="Upper bound of LST range, in hours.")
    a.add_argument("--clobber", default=False, action="store_true", help="Overwrite files.")
    args = a.parse_args()
    return args
