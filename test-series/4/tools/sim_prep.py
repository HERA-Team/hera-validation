import copy
import os
import re
import numpy as np
from scipy import interpolate

from pyuvdata import UVData
from pyuvdata.utils import polstr2num
from hera_cal.io import to_HERAData
from hera_cal.utils import lst_rephase
from hera_sim import noise, Simulator
import hera_sim
if hera_sim.version.version.startswith('0'):
    from hera_sim.rfi import _listify
else:
    from hera_sim.utils import _listify


def add_noise(sim, Trx=100, seed=None):
    """
    Add thermal noise to a UVData object, based on its auto-correlations.

    Parameters
    ----------
    sim : :class:`pyuvdata.UVData` or :class:`hera_sim.Simulator` instance
        The simulation data to add noise to
    Trx : float, optional
        The receiver temperature, in Kelvin.
    seed : int, optional
        The random seed.
    """
    if not isinstance(sim, Simulator):
        sim = Simulator(data=sim)

    lsts = np.unique(sim.data.lst_array)
    freqs = np.unique(sim.data.freq_array)

    # Find out which antenna has the autocorrelation data, since we'll
    # add noise before inflating and adding other systematics.
    antpair = list(
        [ants for ants in sim.data.get_antpairs() if ants[0] == ants[1]]
    )[0]

    # Set up to use the autos to set the noise level
    freqs_GHz = freqs / 1e9  # GHz
    # XXX Has anyone checked that the h1c beam polynomial used in hera_sim
    # actually matches the beam used in the RIMEz simulation?
    omega_p = noise.bm_poly_to_omega_p(freqs_GHz)
    Jy_to_K = noise.jy2T(freqs_GHz, omega_p) / 1000
    
    # XXX The xx and yy polarizations are rotated 90 degrees relative to one
    # another, so does it make sense to use the xx polarization for adding 
    # noise to the data for both the xx and yy polarizations?
    xx_autos = sim.data.get_data(*antpair, 'xx') * Jy_to_K[None, :]
    xx_autos_interp = interpolate.RectBivariateSpline(lsts, freqs_GHz, xx_autos.real)

    yy_autos = sim.data.get_data(*antpair, 'yy') * Jy_to_K[None, :]
    yy_autos_interp = interpolate.RectBivariateSpline(lsts, freqs_GHz, yy_autos.real)

    if seed is not None:
        np.random.seed(seed)
    xx_noise = sim.add_noise(
        'thermal_noise', Tsky_mdl=xx_autos_interp, Trx=Trx, ret_vis=True, add_vis=False
    )
    yy_noise = sim.add_noise(
        'thermal_noise', Tsky_mdl=yy_autos_interp, Trx=Trx, ret_vis=True, add_vis=False
    )
    noise = np.zeros_like(xx_noise, dtype=np.complex)
    
    return sim, noise

def adjust_sim_to_data(sim_file, data_files, save_dir, sky_cmp=None, clobber=True):
    """
    Modify simulated data to be consistent with an observation's metadata.
    
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
        
    clobber : bool, optional
        Whether to overwrite files that may share the same name as the new files to 
        be saved in ``save_dir``. Default is to overwrite files.
        
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
    
    # Get the sky component if not specified.
    sky_cmp = sky_cmp or _parse_filename_for_cmp(sim_file)
    
    # Load in the simulation data, but only the linear vispols.
    use_pols = [polstr2num(pol) for pol in ('xx', 'yy')]
    sim_uvd = UVData()
    sim_uvd.read(sim_file, polarizations=use_pols)

    # Downselect in time and rephase LSTs to match the reference data.
    sim_uvd = rephase_to_reference(sim_uvd, ref_uvd)

    # Inflate the simulation so antenna downselection can actually be done.
    sim_uvd.inflate_by_redundancy()
    
    # Find and use the intersection of the RIMEz and H1C arrays.
    downselect_antennas(sim_uvd, ref_uvd, inplace=True)
    
    # Make sure the data is conjugated properly so redcal doesn't break.
    sim_uvd.conjugate_bls('ant1<ant2')
    
    # Chunk the simulation data and write to disk.
    chunk_sim_and_save(sim_uvd, data_files[0], save_dir, sky_cmp, clobber)
    
    return

def downselect_antennas(sim_uvd, ref_uvd, tol=1.0, inplace=True):
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
        
    inplace : bool, optional
        Whether to perform the downselection in-place. Default is True.
                
    Returns
    -------
    new_sim_uvd : :class:`pyuvdata.UVData`
        Simulation data adjusted to contain a subset of the simulated array
        such that, after appropriately translating the antenna positions, 
        the maximum number of antennas/unique baselines remain in the 
        intersection of the simulation subarray and the reference array. 
        The relevant attributes of the simulation data are updated 
        accordingly. This is only returned if downselection is not performed 
        in-place.

    Notes
    -----
    See the memo "Antennas, Baselines, and Maximum Overlap" for details on 
    the antenna selection process.
    """
    if not inplace:
        sim_uvd = copy.deepcopy(sim_uvd)

    # Find the optimal intersection and get the map between antenna numbers.
    sim_ENU_antpos = _get_antpos(sim_uvd)
    ref_ENU_antpos = _get_antpos(ref_uvd)
    array_intersection = _get_array_intersection(sim_ENU_antpos, ref_ENU_antpos, tol)
    sim_to_ref_ant_map = _get_antenna_map(array_intersection, ref_ENU_antpos, tol)
    
    # Perform a downselection and update the relevant metadata accordingly.
    sim_uvd.select(antenna_nums=list(array_intersection.keys()), keep_all_metadata=False)
    ref_antpos = _get_antpos(ref_uvd, ENU=False)
    sim_antpos = sim_uvd.antenna_positions
    sim_ants = sim_uvd.antenna_numbers
    
    # Prepare new metadata attributes.
    new_sim_ants = np.empty_like(sim_ants)
    new_sim_antpos = np.empty_like(sim_antpos)
    new_ant_1_array = np.empty_like(sim_uvd.ant_1_array)
    new_ant_2_array = np.empty_like(sim_uvd.ant_2_array)
    for sim_ant, ref_ant in sim_to_ref_ant_map.items():
        new_ant_1_array[sim_uvd.ant_1_array == sim_ant] = ref_ant
        new_ant_2_array[sim_uvd.ant_2_array == sim_ant] = ref_ant
        sim_ant_index = sim_ants.tolist().index(sim_ant)
        new_sim_ants[sim_ant_index] = ref_ant
        new_sim_antpos[sim_ant_index] = ref_antpos[ref_ant]
        
    # Actually update the metadata.
    sim_uvd.ant_1_array = new_ant_1_array
    sim_uvd.ant_2_array = new_ant_2_array
    sim_uvd.antenna_numbers = new_sim_ants
    sim_uvd.antenna_positions = new_sim_antpos
    sim_uvd.history += "\nAntennas adjusted to optimally match H1C antennas."
    
    if not inplace:
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
    
    # Find out which times to use and how much to rephase in LST.
    start_lst = ref_lsts[0]
    start_ind = np.argmin(np.abs(sim_lsts - start_lst))
    lst_slice = slice(start_ind, start_ind + ref_uvd.Ntimes)
    use_times = sim_times[lst_slice]
    use_lsts = sim_lsts[lst_slice]
    dlst = ref_lsts - use_lsts

    # Downselect in time and load data.
    sim_uvd.select(times=use_times)
    data, _, _ = sim_uvd.build_datacontainers()

    # Build the antpair -> ENU baseline vector dictionary.
    antpos, ants = sim_uvd.get_ENU_antpos()
    bls = {bl : None for bl in data.bls()}
    ants = list(ants)
    for bl in bls:
        ai, aj = bl[:2]
        i, j = ants.index(ai), ants.index(aj)
        bls[bl] = antpos[j] - antpos[i]

    # Rephase and update the data.
    lst_rephase(data, bls, hd_metas['freqs'], dlst)
    sim_uvd.update(data=data)
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

def chunk_sim_and_save(sim_uvd, ref_file, save_dir, sky_cmp, clobber=True):
    """
    Chunk the simulation data to match the reference file and write to disk.
    
    Parameters
    ----------
    sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the simulation data 
        to chunk and write to disk.
        
    ref_file : str
        Path to a file to use for reference when chunking. It is assumed 
        that all files for the same JD have the same number of 
        integrations as the reference file.
        
    save_dir : str or path-like object
        Path to the directory where the chunked files will be saved.
        
    sky_cmp : str
        String denoting which sky component has been simulated. Should 
        be one of the following: ('foregrounds', 'eor', 'sum').
        
    clobber : bool, optional
        Whether to overwrite any existing files that share the new 
        filenames. Default is to overwrite files.
    """
    # Get some useful metadata.
    uvd = UVData()
    uvd.read(ref_file, read_data=False)
    integrations_per_file = uvd.Ntimes
    times = np.unique(sim_uvd.time_array)
    jd_major = int(np.floor(times[0]))
    Nfiles = int(times.size / integrations_per_file)
    
    # Actually chunk and save the data.
    for Nfile in range(Nfiles):
        start = Nfile * integrations_per_file
        stop = start + integrations_per_file
        use_times = times[start:stop]
        jd_minor = str(times[start] - jd_major)[2:7]
        filename = f"zen.{jd_major}.{jd_minor}.{sky_cmp}.uvh5"
        save_path = os.path.join(save_dir, filename)
        if os.path.exists(save_path):
            if not clobber:
                print(f"File {save_path} exists; skipping.")
                continue
            if clobber:
                print(f"File {save_path} exists; clobbering.")
        this_uvd = sim_uvd.select(times=use_times, inplace=False)
        this_uvd.write_uvh5(save_path, clobber=clobber)
    return

def _get_array_intersection(sim_antpos, ref_antpos, tol=1.0):
    """Find the optimal choice of simulation subarray and return it."""
    optimal_translation = _get_optimal_translation(sim_antpos, ref_antpos, tol)
    new_antpos = {ant : pos + optimal_translation for ant, pos in sim_antpos.items()}
    intersection = {
        ant : pos for ant, pos in new_antpos.items()
        if any([np.allclose(pos, ref_pos, atol=tol) for ref_pos in ref_antpos.values()])
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
        Nintersections = 0
        for pos in new_sim_antpos.values():
            pos_in_ref = any(
                [np.allclose(pos, ref_pos, atol=tol)
                 for ref_pos in ref_antpos.values()]
            )
            Nintersections += pos_in_ref
        intersection_sizes[sim_ant_to_ref_ant] = Nintersections
        
    # Choose the translation that has the most antennas in the intersection.
    sim_to_ref_keys = list(translations.keys())
    intersections_per_translation = np.array(list(intersection_sizes.values()))
    index = np.argwhere(
        intersections_per_translation == intersections_per_translation.max()
    ).flatten()[0]
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
        if not unique_translations:
            unique_translations[key] = translation
            continue
        keep_translation = not any(
            [np.allclose(translation, unique_translation, atol=tol)
             for unique_translation in unique_translations.values()]
        )
        if keep_translation:
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
