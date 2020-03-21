import pathlib
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

def adjust_sim_to_data(sim_file, data_files, save_dir, sky_cmp=None, clobber=True):
    """
    Modify simulated data to be consistent with an observation's metadata.
    
    Parameters
    ----------
    sim_file : str or :class:`pathlib.Path`
        Path to simulation file to modify.
        
    data_files : array-like of str or :class:`pathlib.Path`
        Collection of paths to observed data.
        
    save_dir : str or :class:`pathlib.Path`
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
    # make sure the data files are a list and load in their metadata
    data_files = _listify(data_files)
    ref_uvd = UVData()
    ref_uvd.read_uvh5(data_files, read_data=False)
    
    # get the sky component if not specified
    sky_cmp = sky_cmp or _parse_filename_for_cmp(sim_file)
    
    # downselect in polarization
    use_pols = [polstr2num(pol) for pol in ('xx', 'yy')]
    sim_uvd = UVData()
    sim_uvd.read_uvh5(sim_file, polarizations=use_pols)
    
    # modify the antenna labels and data array for the simulation
    sim_uvd, ref_ants, sim_ants, baseline_map = downselect_antennas(sim_uvd, ref_uvd)
    
    # downselect in time and rephase lsts
    sim_uvd = rephase_to_reference(sim_uvd, ref_uvd)
    
    # chunk the simulation data and save to disk
    chunk_sim_and_save(sim_uvd, data_files[0], save_dir, sky_cmp, clobber)
    
    return ref_ants, sim_ants, baseline_map

def downselect_antennas(sim_uvd, ref_uvd):
    """
    Downselect simulated array antennas to maximally match unique baselines between data sets.
    
    Parameters
    ----------
    sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the simulation data.
        
    ref_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the reference data. Only needs 
        to have the metadata loaded.
                
    Returns
    -------
    new_sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object whose antennas have been relabeled, and whose 
        data array has been rearranged to reflect the change in antenna labels. If 
        the number of antennas to keep is less than those in the simulated data, then 
        a downselection is performed on the simulated data.
        
    ref_ants : list of int
        List of antenna numbers that remain in the reference data.
        
    sim_ants : list of int
        List of antenna numbers that remain in the simulated data.
        
    baseline_map : dict
        Dictionary mapping unique baselines between data sets. The keys of this 
        dictionary are the antenna pairs from the reference data, and the values 
        are the antenna pairs from the simulated data.
    """
    # TODO: actually write the code
    return sim_uvd, None, None, None

def rephase_to_reference(sim_uvd, ref_uvd):
    """
    Rephase simulation LSTs to reference LSTs after downselection in time.
    
    After downselection and rephasing, this function overwrites the simulation LST 
    and time arrays to match the reference LSTs and times.
    
    Parameters
    ----------
    sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the simulation data.
        
    ref_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the reference data. Only the 
        metadata for this object needs to be read.
        
    Returns
    -------
    new_sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the modified simulation data. The 
        modified data has been rephased to match the reference LSTs and has had its 
        time and LST arrays adjusted to contain only the times and LSTs that are 
        present in the reference data.
    """
    # TODO: actually write the code
    return sim_uvd

def chunk_sim_and_save(sim_uvd, ref_file, save_dir, sky_cmp, clobber=True):
    """
    Chunk the simulation data to match the reference file and write to disk.
    
    Parameters
    ----------
    sim_uvd : :class:`pyuvdata.UVData`
        :class:`pyuvdata.UVData` object containing the simulation data to chunk 
        and write to disk.
        
    ref_file : str
        Path to a file to use for reference when chunking. It is assumed that all 
        files for the same JD have the same number of integrations as the reference 
        file.
        
    save_dir : str or :class:`pathlib.Path`
        Path to the directory where the chunked files will be saved.
        
    sky_cmp : str
        String denoting which sky component has been simulated. Should be one of 
        the following: ('foregrounds', 'eor', 'sum').
        
    clobber : bool, optional
        Whether to overwrite any existing files that share the new filenames. Default 
        is to overwrite files.
    """
    # TODO: actually write the code
    pass

def _parse_filename_for_cmp(filename):
    """
    Infer the sky component from the provided filename.
    
    Parameters
    ----------
    filename : str or :class:`pathlib.Path`
        Name of the file to parse; may be a path to the file.
        
    Returns
    -------
    sky_cmp : str
        Inferred sky component. Should be one of the following: ('foregrounds', 'eor', 'sum').
    """
    return None