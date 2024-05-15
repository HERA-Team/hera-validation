import numpy as np
import os
import copy
from pyuvdata import UVData
import hera_sim
from hera_cal.io import HERAData


# All functions here are copied from sim_prep.py from step 4.0
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
    ref_uvd_meta = ref_uvd.copy(metadata_only=True)

    # Find the optimal intersection, get the map between antennas, and downselect the data.
    sim_antpos = _get_antpos(sim_uvd, ENU=False)
    ref_antpos = _get_antpos(ref_uvd_meta, ENU=False)
    array_intersection = _get_array_intersection(sim_antpos, ref_antpos, tol)
    sim_to_ref_ant_map = _get_antenna_map(array_intersection, ref_antpos, tol)
    ref_to_sim_ant_map = {
        ref_ant : sim_ant for sim_ant, ref_ant in sim_to_ref_ant_map.items()
    }
    sim_uvd.select(antenna_nums=list(sim_to_ref_ant_map.keys()), keep_all_metadata=False)
    ref_uvd_meta.select(
        antenna_nums=list(sim_to_ref_ant_map.values()), keep_all_metadata=False
    )
    
    # Update some of the antenna metadata in a copy of the simulation.
    # This is to ensure that slicing is performed correctly.
    sim_uvd_copy = copy.deepcopy(sim_uvd)
    sim_uvd_copy.ant_1_array = np.asarray(
        [sim_to_ref_ant_map[sim_ant] for sim_ant in sim_uvd.ant_1_array]
    )
    sim_uvd_copy.ant_2_array = np.asarray(
        [sim_to_ref_ant_map[sim_ant] for sim_ant in sim_uvd.ant_2_array]
    )
    attrs_to_update = (
        "antenna_numbers",
        "antenna_names",
        "telescope_location",
        "telescope_location_lat_lon_alt",
        "telescope_location_lat_lon_alt_degrees",
    )
    for attr in attrs_to_update:
        setattr(sim_uvd_copy, attr, getattr(ref_uvd_meta, attr))

    # Update the antenna positions array to use shifted simulation positions.
    # array_intersection has simulation antenna numbers as keys and shifted 
    # simulation positions as values.
    sim_uvd_copy.antenna_positions = np.array([
        array_intersection[ref_to_sim_ant_map[ref_ant]]
        for ref_ant in ref_uvd_meta.antenna_numbers
    ])
    sim_uvd_copy.history += "\nAntennas adjusted to optimally match H1C antennas."

    # Prepare the new data array.
    for antpairpol, vis in sim_uvd.antpairpol_iter():
        ai, aj, pol = antpairpol
        ref_antpairpol = (sim_to_ref_ant_map[ai], sim_to_ref_ant_map[aj], pol)
        ref_bl = ref_antpairpol[:2]
        blts, conj_blts, pol_inds = sim_uvd_copy._key2inds(ref_antpairpol)
        
        # Correctly choose which slice to use depending on whether the 
        # reference baseline corresponding to (ai, aj) is conjugated.
        # (If it's conjugated in the reference, then blts is empty.)
        if len(blts) > 0:
            this_slice = (blts, 0, slice(None), pol_inds[0])
        else:
            this_slice = (conj_blts, 0, slice(None), pol_inds[1])
            vis = vis.conj()
            ref_bl = ref_bl[::-1]

        # Update the data-like parameters.
        sim_uvd_copy.data_array[this_slice] = vis
        sim_uvd_copy.flag_array[this_slice] = sim_uvd.get_flags(antpairpol)
        sim_uvd_copy.nsample_array[this_slice] = sim_uvd.get_nsamples(antpairpol)
        # Update the baseline array.
        old_bl_int = sim_uvd.antnums_to_baseline(ai, aj)
        new_bl_int = sim_uvd.antnums_to_baseline(*ref_bl)
        sim_uvd_copy.baseline_array[sim_uvd.baseline_array==old_bl_int] = new_bl_int
        
    # Update the last of the metadata.
    sim_uvd_copy.set_uvws_from_antenna_positions()

    return sim_uvd_copy


# ------- Helper Functions ------- #

def _sim_to_uvd(sim):
    """Update simulation object type."""
    if isinstance(sim, hera_sim.Simulator):
        sim = sim.data
    elif isinstance(sim, HERAData):
        uvd = UVData()
        for attr in uvd:
            setattr(uvd, attr, getattr(sim, attr))
        return uvd
    elif isinstance(sim, str):
        sim_ = UVData()
        sim_.read(sim)
        sim = sim_
    return sim


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


# ------- Main ------- #

if __name__ == '__main__':
    # Paths (uvfits needs some additional information from the original
    # uvfits ref sim data)
    data_path = '/lustre/aoc/projects/hera/Validation/test-neg1.1.0'
    ref_path = os.path.join(data_path, 'ref_sims/ref_1.1_uniform.uvh5')
    rimez_path = os.path.join(data_path, 'rimez_sims/rimez_1.1_uniform_lmax1768.uvh5')
    orig_path = ('/lustre/aoc/projects/hera/ref_sim/simulation_results/' +
                 'ref_1.1_uniform.uvfits')

    ref_uvd, rimez_uvd, uvd_orig = UVData(), UVData(), UVData()
    ref_uvd.read(ref_path)
    rimez_uvd.read(rimez_path)
    uvd_orig.read(orig_path)

    # Relabel the RIMEz antennas to match the ref sim antennas
    # This function was written for validation test 4.0
    relabeled_rimez_uvd = downselect_antennas(copy.deepcopy(rimez_uvd),
                                              copy.deepcopy(ref_uvd))

    # Conjugate the relabeled RIMEz data to the pyuvsim convention (ant2 < ant1)
    relabeled_rimez_uvd.conjugate_bls(convention='ant2<ant1')

    # Downselect the ref sim antennas to the antennas in the RIMEz simulation
    downsampled_ref_uvd = copy.deepcopy(ref_uvd)
    downsampled_ref_uvd.select(antenna_nums=relabeled_rimez_uvd.antenna_numbers,
                               keep_all_metadata=False)

    # Metadata for imaging
    for uvd in [downsampled_ref_uvd, relabeled_rimez_uvd]:
        uvd.object_name = ref_uvd.object_name
        uvd.dut1 = uvd_orig.dut1
        uvd.earth_omega = uvd_orig.earth_omega
        uvd.gst0 = uvd_orig.gst0
        uvd.rdate = uvd_orig.rdate
        uvd.timesys = uvd_orig.timesys
        uvd.phase_center_ra = uvd_orig.phase_center_ra
        uvd.phase_center_dec = uvd_orig.phase_center_dec
        uvd.phase_center_epoch = uvd_orig.phase_center_epoch
        uvd.phase_center_frame = uvd_orig.phase_center_frame
        uvd.channel_width = 1.0  # if this is 0, CASA breaks

    # Save paths
    ref_out = os.path.join(data_path, 'ref_sims/ref_1.1_uniform_downsampled')
    rimez_out = os.path.join(data_path, 'rimez_sims/' +
                             'rimez_1.1_uniform_lmax1768_relabeled')

    # uvh5
    downsampled_ref_uvd.write_uvh5(ref_out + '.uvh5', clobber=True)
    relabeled_rimez_uvd.write_uvh5(rimez_out + '.uvh5', clobber=True)

    # miriad
    ref_miriad_uvd = copy.deepcopy(downsampled_ref_uvd)
    rimez_miriad_uvd = copy.deepcopy(relabeled_rimez_uvd)
    ref_miriad_uvd.write_miriad(ref_out.replace('ref_sims', 'imaging') +
                                '.uv', clobber=True)
    rimez_miriad_uvd.write_miriad(rimez_out.replace('rimez_sims', 'imaging') +
                                  '.uv', clobber=True)

    # uvfits
    ref_uvfits_uvd = copy.deepcopy(downsampled_ref_uvd)
    rimez_uvfits_uvd = copy.deepcopy(relabeled_rimez_uvd)
    ref_uvfits_uvd.write_uvfits(ref_out.replace('ref_sims', 'imaging') +
                                '.uvfits', force_phase=True)
    rimez_uvfits_uvd.write_uvfits(rimez_out.replace('rimez_sims', 'imaging') +
                                  '.uvfits', force_phase=True)
