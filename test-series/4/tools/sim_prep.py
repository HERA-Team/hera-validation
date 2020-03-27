import re
import numpy as np
from scipy import interpolate

from pyuvdata import UVData
from pyuvdata.utils import polstr2num
from hera_cal.io import to_HERAData
from hera_cal.utils import lst_rephase
from hera_sim import noise, Simulator
if hera_sim.version.startswith('0'):
    from hera_sim.rfi import _listify
else:
    from hera_sim.utils import _listify


def add_noise(sim, Trx=100, seed=1010):
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
    autos = sim.data.get_data(*antpair, 'xx') * Jy_to_K[None, :]
    autos_interp = interpolate.RectBivariateSpline(lsts, freqs_GHz, autos.real)

    np.random.seed(seed)
    noise = sim.add_noise(
        'thermal_noise', Tsky_mdl=autos_interp, Trx=Trx, ret_vis=True
    )
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
        
    Returns
    -------
    ref_ants : list of int
        List of reference antennas remaining after antenna downselection.
        
    sim_ants : list of int
        List of simulation antennas remaining after antenna downselection.
        
    baseline_map : dict
        Dictionary mapping unique reference baselines to unique simulation baselines.
    """
    # Ensure the data files are a list and load in their metadata.
    data_files = _listify(data_files)
    ref_uvd = UVData()
    ref_uvd.read(data_files, read_data=False)
    
    # Get the sky component if not specified
    sky_cmp = sky_cmp or _parse_filename_for_cmp(sim_file)
    
    # Load in the simulation data, but only the linear vispols.
    use_pols = [polstr2num(pol) for pol in ('xx', 'yy')]
    sim_uvd = UVData()
    sim_uvd.read(sim_file, polarizations=use_pols)
    
    # Find and use the intersection of the RIMEz and H1C arrays.
    ref_ants, sim_ants, baseline_map = downselect_antennas(sim_uvd, ref_uvd)
    
    # Downselect in time and rephase LSTs to match the reference data.
    rephase_to_reference(sim_uvd, ref_uvd)
    
    # Chunk the simulation data and write to disk.
    chunk_sim_and_save(sim_uvd, data_files[0], save_dir, sky_cmp, clobber)
    
    return ref_ants, sim_ants, baseline_map

def downselect_antennas(sim_uvd, ref_uvd, inplace=True):
    """
    Downselect simulated array to match unique baselines from reference.
    
    Parameters
    ----------
    sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the simulation data.
        
    ref_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the reference data. 
        Only needs to have the metadata loaded.
                
    Returns
    -------
    new_sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object whose antennas have been relabeled, 
        and whose data array has been rearranged to reflect the change in 
        antenna labels. If the number of antennas to keep is less than 
        those in the simulated data, then a downselection is performed on 
        the simulated data.
        
    ref_ants : list of int
        List of antenna numbers that remain in the reference data.
        
    sim_ants : list of int
        List of antenna numbers that remain in the simulated data.
        
    baseline_map : dict
        Dictionary mapping unique baselines between data sets. The keys 
        of this dictionary are the antenna pairs from the reference data, 
        and the values are the antenna pairs from the simulated data.
    """
    # TODO: actually write the code
    ref_ants, sim_ants, baseline_map = [], [], {}

    if not inplace:
        return ref_ants, sim_ants, baseline_map

    return sim_uvd, ref_ants, sim_ants, baseline_map

def rephase_to_reference(sim_uvd, ref_uvd, inplace=True):
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
        
    inplace : bool, optional
        Whether to perform the rephasing and downselection on the 
        simulated data or on a copy of the data. Default is to perform 
        the downselection and rephasing in-place.

    Returns
    -------
    new_sim_uvd : :class:`hera_cal.io.HERAData`
        :class:`hera_cal.io.HERAData` object containing the modified 
        simulation data. The modified data has been rephased to match 
        the reference LSTs and has had its time and LST arrays adjusted 
        to contain only the times and LSTs that are present in the 
        reference data. Only returned if ``inplace`` is set to False.
    """
    # Convert the simulation to a HERAData object. 
    use_uvd = sim_uvd if inplace else copy.deepcopy(sim_uvd)
    hd = to_HERAData(use_uvd)
    hd_metas = hd.get_metadata_dict()

    # Load in useful metadata.
    sim_lsts = hd_metas['lsts']
    sim_times = hd_metas['times']
    ref_lsts = np.unique(ref_uvd.lst_array)
    ref_times = np.unique(ref_uvd.time_array)
    
    # Find out which times to use and how much to rephase in LST
    start_lst = ref_lsts[0]
    dist_from_start = np.abs(sim_lsts - start_lst)
    start_ind = np.argwhere(
        dist_from_start == dist_from_start.min()
    ).flatten()[0]
    lst_slice = slice(start_ind, start_ind + ref_uvd.Ntimes)
    use_times = sim_times[lst_slice]
    use_lsts = sim_lsts[lst_slice]
    dlst = ref_lsts - use_lsts

    # Downselect in time and load data.
    hd.select(times=use_times, keep_all_metadata=False)
    data, _, _ = hd.build_datacontainers()

    # Build the baseline -> ENU baseline vector dictionary.
    antpos, ants = hd.get_ENU_antpos()
    bls = {bl : None for bl in data.bls()}
    ants = list(ants)
    for bl in bls:
        ai, aj = bl[:2]
        i, j = ants.index(ai), ants.index(aj)
        bls[bl] = antpos[i] - antpos[j]

    # Rephase and update the data.
    hera_cal.utils.lst_rephase(data, bls, hd_metas['freqs'], dlst)
    hd.update(data=data)
    new_sim_lsts = np.zeros_like(hd.lst_array)
    new_sim_times = np.zeros_like(hd.time_array)
    for use_time, ref_time, ref_lst in zip(use_times, ref_times, ref_lsts):
        blt_slice = np.argwhere(sim_times == use_time)
        new_sim_lsts[blt_slice] = ref_lst
        new_sim_times[blt_slice] = ref_time
    hd.time_array = new_sim_times
    hd.lst_array = new_sim_lsts

    return hd if not inplace else None

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
    # get some useful metadata
    uvd = UVData()
    uvd.read(ref_file, read_data=False)
    integrations_per_file = uvd.Ntimes
    times = np.unique(sim_uvd.time_array)
    jd_major = int(np.floor(times[0]))
    Nfiles = int(times.size / integrations_per_file)
    
    # actually chunk and save the data
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
        this_uvd.write_uvh5(save_path)
    return

def _parse_filename_for_cmp(filename):
    """Infer the sky component from the provided filename."""
    cmp_re = re.compile("[a-z]+.uvh5")
    sky_cmp = cmp_re.findall(str(filename))
    if sky_cmp == []:
        raise ValueError(
            f"Simulation component could not be inferred from {filename}."
        )
    return sky_cmp[0][:-5]
