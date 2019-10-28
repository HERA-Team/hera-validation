#!/usr/bin/env python2
"""
skycal_pipe.py
-----------------------------------------
This script is used as an automatic calibration
and imaging pipeline in CASA for HERA data.

Run as python skycal_pipe.py -c skycal_params.yml <options>

See skycal_params.yml for relevant parameter selections.

Nicholas Kern
nkern@berkeley.edu
November, 2018
"""
from __future__ import print_function, division, absolute_import

import numpy as np
from pyuvdata import UVData
import pyuvdata.utils as uvutils
import casa_imaging
from casa_imaging import casa_utils as utils
import os
import sys
import glob
import yaml
from datetime import datetime
import json
import itertools
import shutil
from collections import OrderedDict as odict
from astropy.time import Time
from functools import reduce
import copy
import operator
import subprocess
import argparse

# get casa_imaging path
casa_scripts = casa_imaging.__path__[0]
casa_scripts = os.path.join('/'.join(casa_scripts.split('/')[:-1]), 'scripts')

#-------------------------------------------------------------------------------
# Define Optional Command-Line Parameters
#-------------------------------------------------------------------------------
args = argparse.ArgumentParser(description="skycal_pipe.py: run as python 2.2skycal_pipe.py -c 2.2skycal_params.yml <opts>")
args.add_argument("-c", "--param_file", type=str, help='Path to a YAML parameter file. See skycal_params.yml for details.')
# Optional Parameters that, if defined, overwrite their counterpart in param_file.yml
args.add_argument("--data_root", default=None, type=str, help="Root path to data files: overwrites skycal_params.yml")
args.add_argument("--data_file", default=None, type=str, help="Data file basename: overwrites skycal_params.yml")
args.add_argument("--source", default=None, type=str, help="Source name: overwrites skycal_params.yml")
args.add_argument("--source_ra", default=None, type=float, help="Source right ascension in J2000 degrees: overwrites skycal_params.yml ")
args.add_argument("--source_dec", default=None, type=float, help="Source declination in J2000 degrees: overwrites skycal_params.yml")
a = args.parse_args()

#-------------------------------------------------------------------------------
# Parse YAML Configuration File
#-------------------------------------------------------------------------------
# Get config and load dictionary
cf = utils.load_config(a.param_file)#yaml.load(a.param_file, Loader=yaml.FullLoader)#utils.load_config(a.param_file)

# Consolidate IO, data and analysis parameter dictionaries
params = odict(list(cf['io'].items()) + list(cf['data'].items()) + list(cf['analysis'].items())+ list(cf['obs'].items()))
#assert (len(params) == len(cf['io']) + len(cf['data']) + len(cf['analysis']) +list(cf['obs'].items())), "Repeated parameters found within the scope of io, data and analysis dicts"

# Get algorithm dictionary
algs = cf['algorithm']
datafile = os.path.join(params['data_root'], params['data_file'])

# Get parameters used globally in the pipeline
verbose = params['verbose']
overwrite = params['overwrite']
casa = params['casa'].split() + params['casa_flags'].split()
longitude = params['longitude']
latitude = params['latitude']
out_dir = params['out_dir']
source_ra = params['source_ra']
source_dec = params['source_dec']
source = params['source']

# Change to working dir
os.chdir(params['work_dir'])

# Open a logfile
logfile = os.path.join(out_dir, params['logfile'])
if os.path.exists(logfile) and params['overwrite'] == False:
    raise IOError("logfile {} exists and overwrite == False, quitting pipeline...".format(logfile))
lf = open(logfile, "w")
if params['joinlog']:
    ef = lf
else:
    ef = open(os.path.join(params['out_dir'], params['errfile']), "w")
casa += ['--logfile', logfile]
sys.stdout = lf
sys.stderr = ef

# Setup (Small) Global Variable Dictionary
varlist = ['datafile', 'verbose', 'overwrite', 'out_dir', 'casa', 'point_ra', 'longitude',
           'latitude', 'lf', 'gaintables']
def global_vars(varlist=[]):
    d = []
    for v in varlist:
        try:
            d.append((v, globals()[v]))
        except KeyError:
            continue
    return dict(d)

# Print out parameter header
time = datetime.utcnow()
utils.log("Starting sky calibration pipeline on {}\n{}\n".format(time, '-'*60),
          f=lf, verbose=verbose)
_cf = copy.copy(cf)
_cf.pop('algorithm')
utils.log(json.dumps(_cf, indent=1) + '\n', f=lf, verbose=verbose)

# Setup a dict->object converter
class Dict2Obj:
    def __init__(self, **entries):
        self.__dict__.update(entries)

#-------------------------------------------------------------------------------
# Search for a Source and Prepare Data for MS Conversion
#-------------------------------------------------------------------------------
if params['prep_data']:
    # start block
    time = datetime.utcnow()
    utils.log("\n{}\n...Starting PREP_DATA: {}\n".format("-"*60, time),
    f=lf, verbose=verbose)
    utils.log(json.dumps(algs['prep_data'], indent=1) + '\n', f=lf, verbose=verbose)
    p = Dict2Obj(**algs['prep_data'])

    # Check if datafile is already MS
    sys.path.insert(0,'/Users/tashaleebillings/casa_imaging/scripts/')
    import source2file

    if os.path.splitext(datafile)[1] in ['.ms', '.MS']:
        # get transit time of source
        if source_ra is not None:
            (lst, transit_jd, utc_range, utc_center, source_files,
            source_utc_range) = source2file.source2file(source_ra, lon=longitude, lat=latitude,
                                                   duration=p.duration, start_jd=p.start_jd, get_filetimes=False,
                                                   verbose=verbose)
            utils.log("...the file {} is already a CASA MS, skipping rest of PREP_DATA".format(datafile), f=lf, verbose=verbose)
            timerange = utc_range
        else:
            transit_jd = None

    else:
        # Iterate over polarizations
        if p.pols is None: p.pols = [None]
        uvds = []
        for pol in p.pols:
            if pol is None:
                pol = ''
            else:
                utils.log("...working on {} polarization".format(pol), f=lf, verbose=verbose)
                pol = '.{}.'.format(pol)
        
            # glob-parse the data file / template
            datafiles = [df for df in glob.glob(datafile) if pol in df]
            assert len(datafiles) > 0, "Searching for {} with pol {} but found no files...".format(datafile, pol)
            
            # get transit times
            if source_ra is not None:
                (lst, transit_jd, utc_range, utc_center, source_files,
                 source_utc_range) = source2file.source2file(source_ra, lon=longitude, lat=latitude,
                                                             duration=p.duration, start_jd=p.start_jd, get_filetimes=p.get_filetimes,
                                                             verbose=verbose, jd_files=copy.copy(datafiles))
                timerange = utc_range
                                                             
                # ensure source_utc_range and utc_range are similar
                if source_utc_range is not None:
                    utc_range_start = utc_range.split('~')[0].strip('"').split('/')
                    utc_range_start = map(int, utc_range_start[:-1] + utc_range_start[-1].split(':'))
                    utc_range_start = Time(datetime(*utc_range_start), format='datetime').jd
                    source_utc_range_start = source_utc_range.split('~')[0].strip('"').split('/')
                    source_utc_range_start = map(int, source_utc_range_start[:-1] + source_utc_range_start[-1].split(':'))
                    source_utc_range_start = Time(datetime(*source_utc_range_start), format='datetime').jd
                    # if difference is larger than 1 minute,
                    # then probably the correct files were not found
                    if np.abs(utc_range_start - source_utc_range_start) * 24 * 60 > 1:
                        utils.log("Warning: Difference between theoretical transit time and transit time " \
                               "deduced from files found is larger than 1-minute: probably the correct " \
                               "files were not found because the correct files did not exist under the " \
                               "data template {}".format(datafile), f=lf, verbose=verbose)
                    timerange = source_utc_range
            else:
                source_files = datafiles
        
            # load data into UVData
            utils.log("...loading data", f=lf, verbose=verbose)
            _uvds = []
            for sf in list(source_files):
                # read data
                _uvd = UVData()
                _uvd.read(sf, antenna_nums=p.antenna_nums)
                
                # read flagfile if fed
                if p.flag_ext != "":
                    flagfile = glob.glob("{}{}".format(sf, p.flag_ext))
                    if len(flagfile) == 1:
                        utils.log("...loading and applying flags {}".format(flagfile[0]), f=lf, verbose=verbose)
                        ff = np.load(flagfile[0])
                        _uvd.flag_array += ff['flag_array']
        
                # append to list
                _uvds.append(_uvd)
            
            # concatenate source files
            uvd = reduce(operator.add, _uvds)
            
            # isolate only relevant times
            times = np.unique(uvd.time_array)
            if source_ra is not None:
                transit_jd = np.mean(times)
            times = times[np.abs(times - transit_jd) < (p.duration / (24. * 60. * 2))]
            assert len(times) > 0, "No times found in source_files {} given transit JD {} and duration {}".format(source_files, transit_jd, p.duration)
            uvd.select(times=times)
            
            # append
            uvds.append(uvd)

        # concatenate uvds
        uvd = reduce(operator.add, uvds)

        # get output filepath w/o uvfits extension if provided
        outfile = os.path.join(params['out_dir'], p.outfile.format(uvd.time_array.min()))
        if os.path.splitext(outfile)[1] == '.uvfits':
            outfile = os.path.splitext(outfile)[0]
        
        # renumber antennas (and antenna names!) if above 254
        if uvd.antenna_numbers.max() > 254:
            large_ant_nums = sorted(list(uvd.antenna_numbers[np.where(uvd.antenna_numbers > 254)[0]]))
            new_nums = sorted(list(set(range(255)) - set(uvd.antenna_numbers)))
            if len(new_nums) < len(large_ant_nums):
                raise ValueError('too many antennas in dataset, cannot renumber all below 255')
            new_nums = new_nums[-1 * len(large_ant_nums):]
            renumber_dict = dict(list(zip(large_ant_nums, new_nums)))
            
            history = ''
            name_prefix = os.path.commonprefix(uvd.antenna_names)
            for ant_in, ant_out in renumber_dict.items():
                if verbose:
                    msg = "renumbering {a1} to {a2}".format(a1=ant_in, a2=ant_out)
                    print(msg)
                history += '{}\n'.format(msg)
                
                wh_ant_num = np.where(uvd.antenna_numbers == ant_in)[0]
                wh_ant1_arr = np.where(uvd.ant_1_array == ant_in)[0]
                wh_ant2_arr = np.where(uvd.ant_2_array == ant_in)[0]
                
                uvd.antenna_numbers[wh_ant_num] = ant_out
                uvd.antenna_names[wh_ant_num[0]] = "RN{:d}".format(ant_out)
                uvd.ant_1_array[wh_ant1_arr] = ant_out
                uvd.ant_2_array[wh_ant2_arr] = ant_out
            
            uvd.baseline_array = uvd.antnums_to_baseline(uvd.ant_1_array, uvd.ant_2_array)
            uvd.history = '{}\n{}'.format(history, uvd.history)
            uvd.check()
            
            # write renumbering dictionary to .npz
            np.savez("{}.renumber.npz".format(outfile),
                     renumber=dict(zip(renumber_dict.values(), renumber_dict.keys())),
                     history="Access dictionary via f['renumber'].item()")
        
        # write to file
        if True: #uvd.phase_type == 'phased':
            # write uvfits
            uvfits_outfile = outfile + '.uvfits'
            if not os.path.exists(uvfits_outfile) or overwrite:
                utils.log("...writing {}".format(uvfits_outfile), f=lf, verbose=verbose)
                uvd.write_uvfits(uvfits_outfile, spoof_nonessential=True, force_phase=True)
        # convert the uvfits file to ms
        ms_outfile = outfile + '.ms'
        utils.log("...converting to Measurement Set")
        if not os.path.exists(ms_outfile) or overwrite:
            if os.path.exists(ms_outfile):
                shutil.rmtree(ms_outfile)
            utils.log("...writing {}".format(ms_outfile), f=lf, verbose=verbose)
            ecode = subprocess.call(casa + ["-c", "importuvfits('{}', '{}')".format(uvfits_outfile, ms_outfile)])
        
        # overwrite relevant parameters for downstream analysis
        datafile = ms_outfile

        del uvds, uvd

    # overwrite downstream parameters
    algs['gen_model']['time'] = transit_jd

    # end block
    time2 = datetime.utcnow()
    utils.log("...finished PREP_DATA: {:d} sec elapsed".format(utils.get_elapsed_time(time, time2)), f=lf, verbose=verbose)


#-------------------------------------------------------------------------------
# Generate Flux Model
#-------------------------------------------------------------------------------

# pbcorr defined function
def pbcorr(modelname):
    
    import numpy as np
    import astropy.io.fits as fits
    from astropy import wcs
    from pyuvdata import UVBeam, utils as uvutils
    import os
    import sys
    import glob
    import argparse
    import shutil
    import copy
    import healpy
    import scipy.stats as stats
    from casa_imaging import casa_utils
    from scipy import interpolate
    from astropy.time import Time
    from astropy import coordinates as crd
    from astropy import units as u
    
    
    _fitsfiles = ["{}.fits".format(modelname)]
    
    # PB args
    _multiply = True
    _lon = p.longitude
    _lat = p.latitude
    _time = p.time
    
    # beam args
    _beamfile = p.beamfile
    _pols = -5, -6
    _freq_interp_kind = 'cubic'
    
    # IO args
    _ext = ''
    _outdir = p.out_dir
    _overwrite = True
    _silence = False
    _spec_cube = False
    
    def echo(message, type=0):
        if verbose:
            if type == 0:
                print(message)
            elif type == 1:
                print('\n{}\n{}'.format(message, '-'*40))

    verbose = _silence == False
    
    # load pb
    echo("...loading beamfile {}".format(_beamfile))
    # load beam
    uvb = UVBeam()
    uvb.read_beamfits(_beamfile)
    if uvb.pixel_coordinate_system == 'healpix':
        uvb.interpolation_function = 'healpix_simple'
    else:
        uvb.interpolation_function = 'az_za_simple'
    uvb.freq_interp_kind = _freq_interp_kind
    
    # get beam models and beam parameters
    beam_freqs = uvb.freq_array.squeeze() / 1e6
    Nbeam_freqs = len(beam_freqs)
    
    # iterate over FITS files
    for i, ffile in enumerate(_fitsfiles):
        
        # create output filename
        if _outdir is None:
            output_dir = os.path.dirname(ffile)
        else:
            output_dir = _outdir
        
        output_fname = os.path.basename(ffile)
        output_fname = os.path.splitext(output_fname)
        if _ext is not None:
            output_fname = output_fname[0] + '.pbcorr{}'.format(_ext) + output_fname[1]
        else:
            output_fname = output_fname[0] + '.pbcorr' + output_fname[1]
        output_fname = os.path.join(output_dir, output_fname)
        
        # check for overwrite
        if os.path.exists(output_fname) and _overwrite is False:
            raise IOError("{} exists, not overwriting".format(output_fname))

        # load hdu
        echo("...loading {}".format(ffile))
        hdu = fits.open(ffile)
    
        # get header and data
        head = hdu[0].header
        data = hdu[0].data
        
        # get polarization info
        ra, dec, pol_arr, data_freqs, stok_ax, freq_ax = casa_utils.get_hdu_info(hdu)
        Ndata_freqs = len(data_freqs)
        
        # get axes info
        npix1 = head["NAXIS1"]
        npix2 = head["NAXIS2"]
        nstok = head["NAXIS{}".format(stok_ax)]
        nfreq = head["NAXIS{}".format(freq_ax)]
        
        # replace with forced polarization if provided
        if _pols is not None:
            pol_arr = np.asarray(_pols, dtype=np.int)
    
        pols = [uvutils.polnum2str(pol) for pol in pol_arr]
        
        # make sure required pols exist in maps
        if not np.all([p in uvb.polarization_array for p in pol_arr]):
            raise ValueError("Required polarizationns {} not found in Beam polarization array".format(pol_arr))

        # convert from equatorial to spherical coordinates
        loc = crd.EarthLocation(lat=_lat*u.degree, lon=_lon*u.degree)
        time = Time(_time, format='jd', scale='utc')
        equatorial = crd.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5', location=loc, obstime=time)
        altaz = equatorial.transform_to('altaz')
        theta = np.abs(altaz.alt.value - 90.0)
        phi = altaz.az.value

        # convert to radians
        theta *= np.pi / 180
        phi *= np.pi / 180
        
        if i == 0 or _spec_cube is False:
            # evaluate primary beam
            echo("...evaluating PB")
            pb, _ = uvb.interp(phi.ravel(), theta.ravel(), polarizations=pols, reuse_spline=True)
            pb = np.abs(pb.reshape((len(pols), Nbeam_freqs) + phi.shape))

        # interpolate primary beam onto data frequencies
        echo("...interpolating PB")
        pb_shape = (pb.shape[1], pb.shape[2])
        pb_interp = interpolate.interp1d(beam_freqs, pb, axis=1, kind=_freq_interp_kind, fill_value='extrapolate')(data_freqs / 1e6)
        
        # data shape is [naxis4, naxis3, naxis2, naxis1]
        if freq_ax == 4:
            pb_interp = np.moveaxis(pb_interp, 0, 1)

        # divide or multiply by primary beam
        if _multiply is True:
            echo("...multiplying PB into image")
            data_pbcorr = data * pb_interp
        else:
            echo("...dividing PB into image")
            data_pbcorr = data / pb_interp

        # change polarization to interpolated beam pols
        head["CRVAL{}".format(stok_ax)] = pol_arr[0]
        if len(pol_arr) == 1:
            step = 1
        else:
            step = np.diff(pol_arr)[0]
        head["CDELT{}".format(stok_ax)] = step
        head["NAXIS{}".format(stok_ax)] = len(pol_arr)

        echo("...saving {}".format(output_fname))
        fits.writeto(output_fname, data_pbcorr, head, overwrite=True)

        output_pb = output_fname.replace(".pbcorr.", ".pb.")
        echo("...saving {}".format(output_pb))
        fits.writeto(output_pb, pb_interp, head, overwrite=True)
        
        return

# Make a Model Generation Function
p = Dict2Obj(**dict(list(algs['gen_model'].items()) + list(global_vars(varlist).items())+list(cf['obs'].items())))
def gen_model(**kwargs):
    p = Dict2Obj(**kwargs)
    utils.log("\n{}\n...Generating a Flux Model", f=p.lf, verbose=p.verbose)
    
    # compile complist_gleam.py command
    cmd = casa + ["-c", "{}/complist_gleam.py".format(casa_scripts)]
    cmd += ['--point_ra', p.source_ra, '--point_dec', p.latitude, '--outdir', p.out_dir,
            '--gleamfile', p.gleamfile, '--radius', p.radius, '--min_flux', p.min_flux,
            '--freqs', p.freqs, '--cell', p.cell, '--imsize', p.imsize]
    if p.image:
        cmd += ['--image']
    if p.use_peak:
        cmd += ['--use_peak']
    if p.overwrite:
        cmd += ['--overwrite']
    if hasattr(p, 'regions'):
        cmd += ['--regions', p.regions, '--exclude', '--region_radius', p.region_radius]
    if hasattr(p, 'file_ext'):
        cmd += ['--ext', p.file_ext]
    else:
        p.file_ext = ''
    cmd = map(str, cmd)
    ecode = subprocess.call(cmd)
    
    modelstem = os.path.join(p.out_dir, "gleam{}.cl".format(p.file_ext))
    model = modelstem
    if p.image:
        model += ".image"

    print(p.longitude)
    print(p.latitude)
    print(p.time)

    # beam args
    print(p.beamfile)
    print(p.out_dir)

    # pbcorrect
    if p.pbcorr:
        # Do primary beam correction
        pbcorr(modelname=modelstem)
        # importfits
        cmd = p.casa + ["-c", "importfits('{}', '{}', overwrite={})".format(modelstem + '.pbcorr.fits', modelstem + '.pbcorr.image', p.overwrite)]
        ecode = subprocess.call(cmd)
        model = modelstem + ".pbcorr.image"

    return model

if params['gen_model']:
    # start block
    time = datetime.utcnow()
    utils.log("\n{}\n...Starting GEN_MODEL: {}\n".format("-"*60, time),
    f=lf, verbose=verbose)
    utils.log(json.dumps(algs['gen_model'], indent=1) + '\n', f=lf, verbose=verbose)

    # Generate Model
    kwvalues = dict(list(algs['gen_model'].items()) + list(global_vars(varlist).items())+list(cf['obs'].items()))
    model = gen_model(**kwvalues)

    # update di_cal model path
    algs['di_cal']['model'] = model

    # end block
    time2 = datetime.utcnow()
    utils.log("...finished GEN_MODEL: {:d} sec elapsed".format(utils.get_elapsed_time(time, time2)), f=lf, verbose=verbose)


#-------------------------------------------------------------------------------
# Calibration and Imaging Functions
#-------------------------------------------------------------------------------
# Define a calibration Function
def calibrate(**cal_kwargs):
    kwargs = dict(list(algs['gen_cal'].items()) + list(global_vars(varlist).items())+list(cf['obs'].items()) + list(algs['di_cal'].items()))
    p = Dict2Obj(**kwargs)

    # compile command
    cmd = p.casa + ["-c", "{}/sky_cal.py".format(casa_scripts)]
    cmd += ["--msin", p.datafile, "--out_dir", p.out_dir, "--model", p.model,
            "--refant", p.refant, "--gain_spw", p.gain_spw, "--uvrange", p.uvrange, "--timerange",
            p.timerange, "--ex_ants", p.ex_ants, "--gain_ext", p.gain_ext, '--bp_spw', p.bp_spw]
    if p.source_ra is not None:
        cmd += ["--source_ra", p.source_ra]
    if p.source_dec is not None:
        cmd += ["--source_dec", p.source_dec]

    if isinstance(p.gaintables, list):
        gtables = p.gaintables
    else:
        if p.gaintables in ['', None, 'None', 'none']:
            gtables = []
        else:
            gtables = [p.gaintables]
    if len(gtables) > 0:
        cmd += ['--gaintables'] + gtables

    if p.rflag:
        cmd += ["--rflag"]
    if p.KGcal:
        cmd += ["--KGcal", "--KGsnr", p.KGsnr]
    if p.Acal:
        cmd += ["--Acal", "--Asnr", p.Asnr]
    if p.BPcal:
        cmd += ["--BPcal", "--BPsnr", p.BPsnr, "--bp_spw", p.bp_spw]
    if p.BPsolnorm:
        cmd += ['--BPsolnorm']
    if p.split_cal:
        cmd += ["--split_cal", "--cal_ext", p.cal_ext]
    if p.split_model:
        cmd += ["--split_model"]
    cmd = [' '.join(_cmd) if type(_cmd) == list else str(_cmd) for _cmd in cmd]

    utils.log("...starting calibration", f=p.lf, verbose=p.verbose)
    ecode = subprocess.call(cmd)

    # Gather gaintables
    gext = ''
    if p.gain_ext not in ['', None]:
        gext = '.{}'.format(p.gain_ext)
    gts = sorted(glob.glob("{}{}.?.cal".format(p.datafile, gext)) + glob.glob("{}{}.????.cal".format(p.datafile, gext)))

    # export to calfits if desired
    if p.export_gains:
        utils.log("...exporting\n{}\n to calfits and combining into a single cal table".format('\n'.join(gts)), f=p.lf, verbose=p.verbose)
        # do file checks
        mirvis = os.path.splitext(p.datafile)[0]
        gtsnpz = ["{}.npz".format(gt) for gt in gts]
        if not os.path.exists(mirvis):
            utils.log("...{} doesn't exist: cannot export gains to calfits".format(mirvis), f=p.lf, verbose=p.verbose)

        elif len(gts) == 0:
            utils.log("...no gaintables found, cannot export gains to calfits", f=p.lf, verbose=p.verbose)

        elif not np.all([os.path.exists(gtnpz) for gtnpz in gtsnpz]):
            utils.log("...couldn't find a .npz file for all input .cal tables, can't export to calfits", f=p.lf, verbose=p.verbose)

        else:
            calfits_fname = "{}.{}{}.calfits".format(os.path.basename(mirvis), p.source, p.gain_ext)
            cmd = ['skynpz2calfits.py', "--fname", calfits_fname, "--uv_file", mirvis, '--out_dir', p.out_dir]
            if p.overwrite:
                cmd += ["--overwrite"]
            # add a delay and phase solution
            matchK = ["K.cal.npz" in gt for gt in gtsnpz]
            matchGphs = ["Gphs.cal.npz" in gt for gt in gtsnpz]
            if np.any(matchK):
                cmd += ["--plot_dlys"]
                cmd += ["--dly_files"] + [gtsnpz[i] for i, b in enumerate(matchK) if b == True]
                if not np.any(matchGphs):
                    utils.log("...WARNING: A delay file {} was found, but no mean phase file, which is needed if a delay file is present.", f=lf, verbose=verbose)
            if np.any(matchGphs):
                cmd += ["--plot_phs"]
                cmd += ["--phs_files"] + [gtsnpz[i] for i, b in enumerate(matchGphs) if b == True]

            # add a mean amp solution
            matchGamp = ["Gamp.cal.npz" in gt for gt in gtsnpz]
            if np.any(matchGamp):
                cmd += ["--plot_amp"]
                cmd += ["--amp_files"] + [gtsnpz[i] for i, b in enumerate(matchGamp) if b == True]

            # add a bandpass solution
            matchB = ["B.cal.npz" in gt for gt in gtsnpz]
            if np.any(matchB):
                cmd += ["--plot_bp"]
                cmd += ["--bp_files"] + [gtsnpz[i] for i, b in enumerate(matchB) if b == True]

            # additional smoothing options
            if p.smooth:
                cmd += ["--bp_gp_smooth", "--bp_gp_max_dly", p.gp_max_dly]
            if p.medfilt:
                cmd += ["--bp_medfilt", "--medfilt_kernel", p.kernel]
            if p.bp_broad_flags:
                cmd += ["--bp_broad_flags", "--bp_flag_frac", p.bp_flag_frac]
            if not p.verbose:
                cmd += ['--silence']

            cmd = map(str, cmd)
            ecode = subprocess.call(cmd)

            # convert calfits back to a single Btotal.cal table
            if np.any(matchB):
                # convert to cal
                bfile = gts[matchB.index(True)]
                btot_file = os.path.join(out_dir, "{}{}.Btotal.cal".format(os.path.basename(p.datafile), gext))
                cmd = p.casa + ["-c", "{}/calfits_to_Bcal.py".format(casa_scripts), "--cfits", os.path.join(p.out_dir, calfits_fname), "--inp_cfile", bfile,"--out_cfile", btot_file]
                if overwrite:
                    cmd += ["--overwrite"]
                ecode = subprocess.call(cmd)
                # replace gaintables with Btotal.cal
                gts = [btot_file]

    # append to gaintables
    gtables += gts

    return gtables

# Define imaging functions
def img_cmd(**kwargs):
    p = Dict2Obj(**kwargs)
    cmd = p.casa + ["-c", "{}/sky_image.py".format(casa_scripts)]
    cmd += ["--out_dir", p.out_dir,
            "--pxsize", p.pxsize, "--imsize", p.imsize,
            "--uvrange", p.uvrange, "--timerange", p.timerange,
            "--stokes", p.stokes, "--weighting", p.weighting, "--robust", p.robust,
            "--pblimit", p.pblimit, "--deconvolver", p.deconvolver, "--niter",
            p.niter, '--cycleniter', p.cycleniter, '--threshold', p.threshold,
            '--mask', p.mask, '--gridder', p.gridder, '--wprojplanes', p.wpplanes]
    if p.source is not None:
        cmd += ["--source", p.source]
    if p.source_ra is not None:
        cmd += ["--source_ra", p.source_ra]
    if p.source_dec is not None:
        cmd += ["--source_dec", p.source_dec]
    cmd = [map(str, _cmd) if type(_cmd) == list else str(_cmd) for _cmd in cmd]
    cmd = reduce(operator.add, [i if type(i) == list else [i] for i in cmd])
    return cmd, p

def mfs_image(**kwargs):
    cmd, p = img_cmd(**kwargs)

    # Perform MFS imaging
    utils.log("...starting MFS image of {} data".format(p.mfstype), f=p.lf, verbose=p.verbose)
    icmd = cmd + ['--image_mfs', '--msin', p.datafile, '--spw', p.spw]
    if p.mfstype == 'resid':
        icmd += ['--uvsub']
    if p.source_ext == '':
        source_ext = ''
    else:
        source_ext = '{}_'.format(p.source_ext)
    if p.mfstype == '':
        mfstype = ''
    else:
        mfstype = '_{}'.format(p.mfstype)
    icmd += ['--source_ext', "{}{}".format(source_ext, mfstype)]

    ecode = subprocess.call(icmd)

    if p.mfstype == 'resid':
        # Apply gaintables to make CORRECTED column as it was
        utils.log("...reapplying gaintables to CORRECTED data", f=p.lf, verbose=p.verbose)
        cmd2 = p.casa + ["-c", "{}/sky_cal.py".format(casa_scripts), "--msin", p.datafile, "--gaintables"] + p.gaintables
        ecode = subprocess.call(cmd2)

def spec_image(**kwargs):
    cmd, p = img_cmd(**kwargs)

    # Perform Spectral Cube imaging
    utils.log("...starting {} spectral cube imaging".format(p.datafile), f=p.lf, verbose=p.verbose)
    icmd = cmd + ['--spec_cube', '--msin', p.datafile,
                  '--spec_start', str(p.spec_start), '--spec_end', str(p.spec_end),
                  '--spec_dchan', str(p.spec_dchan)]
    if p.source_ext == '':
        source_ext = ''
    else:
        source_ext = '{}_'.format(p.source_ext)
    icmd += ['--source_ext', "{}spec".format(source_ext)]

    ecode = subprocess.call(icmd)

    # Collate output images and run a source extraction
    img_cube_template = "{}.{}spec.chan????.image.fits".format(p.datafile, source_ext)
    img_cube = sorted(glob.glob(img_cube_template))
    if p.source_extract:
        utils.log("...extracting {} from {} cube".format(p.source, img_cube_template), f=p.lf, verbose=p.verbose)
        if len(img_cube) == 0:
            utils.log("...no image cube files found, cannot extract spectrum", f=p.lf, verbose=p.verbose)
        else:
            cmd = ["source_extract.py", "--source", p.source, "--source_ra", p.source_ra, "--source_dec", p.source_dec,
                   "--radius", p.radius, '--pols'] + p.pols \
                   + ["--outdir", p.out_dir, "--gaussfit_mult", p.gauss_mult, "--source_ext", p.source_ext]
            if p.overwrite:
                cmd += ["--overwrite"]
            if p.plot_fit:
                cmd += ["--plot_fit"]
            cmd += img_cube
            cmd = map(str, cmd)
            ecode = subprocess.call(cmd)

# generalized MFS + spectral imaging function
def image(**img_kwargs):
    # Perform MFS of corrected data
    if img_kwargs['image_mfs']:
        kwargs = dict(list(img_kwargs.items()) + list(global_vars(varlist).items()))
        kwargs['mfstype'] = 'corr'
        mfs_image(**kwargs)

    # Perform MFS of model data
    if img_kwargs['image_mdl']:
        kwargs = dict(list(img_kwargs.items()) + list(global_vars(varlist).items()))
        mfile = "{}.model".format(datafile)
        if not os.path.exists(mfile):
            utils.log("Didn't split model from datafile, which is required to image the model", f=lf, verbose=verbose)
        else:
            kwargs['datafile'] = mfile
            kwargs['mfstype'] = 'model'
            mfs_image(**kwargs)

    # Perform MFS of residual data
    if img_kwargs['image_res']:
        kwargs = dict(list(img_kwargs.items()) + list(global_vars(varlist).items()))
        kwargs['mfstype'] = 'resid'
        mfs_image(**kwargs)

    # Get spectral cube of corrected data
    if img_kwargs['image_spec']:
        kwargs = dict(list(img_kwargs.items()) + list(global_vars(varlist).items()))
        spec_image(**kwargs)

    # Get spectral cube of model data
    if img_kwargs['image_mdl_spec']:
        kwargs = dict(list(img_kwargs.items()) + list(global_vars(varlist).items()))
        mfile = "{}.model".format(datafile)
        if not os.path.exists(mfile):
            utils.log("Didn't split model from datafile, which is required to image the model", f=lf, verbose=verbose)
        else:
            kwargs['datafile'] = mfile
            spec_image(**kwargs)

#-------------------------------------------------------------------------------
# Direction Independent Calibration
#-------------------------------------------------------------------------------
if params['di_cal']:
    # start block
    time = datetime.utcnow()
    utils.log("\n{}\n...Starting DI_CAL: {}\n".format("-"*60, time), f=lf, verbose=verbose)
    cal_kwargs = dict(list(algs['gen_cal'].items()) + list(global_vars(varlist).items())+list(cf['obs'].items()) + list(algs['di_cal'].items()))
    #utils.log(json.dumps(cal_kwargs, indent=1) + '\n', f=lf, verbose=verbose)

    # Perform Calibration
    kwargs = global_vars(varlist)
    kwargs.update(cal_kwargs)
    gaintables = calibrate(**kwargs)

    # end block
    time2 = datetime.utcnow()
    utils.log("...finished DI_CAL: {:d} sec elapsed".format(utils.get_elapsed_time(time, time2)), f=lf, verbose=verbose)

#-------------------------------------------------------------------------------
# Imaging
#-------------------------------------------------------------------------------
if params['di_img']:
    # start block
    time = datetime.utcnow()
    utils.log("\n{}\n...Starting DI_IMG: {}\n".format("-"*60, time), f=lf, verbose=verbose)
    img_kwargs = dict(list(algs['gen_cal'].items()) + list(global_vars(varlist).items()) + list(cf['obs'].items()) + list(algs['imaging'].items()) + list(algs['di_cal'].items()))
    utils.log(json.dumps(img_kwargs, indent=1) + '\n', f=lf, verbose=verbose)

    # Peform Imaging
    image(**img_kwargs)

    # end block
    time2 = datetime.utcnow()
    utils.log("...finished DI_IMG: {:d} sec elapsed".format(utils.get_elapsed_time(time, time2)), f=lf, verbose=verbose)

#-------------------------------------------------------------------------------
# Direction Dependent Calibration
#-------------------------------------------------------------------------------
if params['dd_cal']:
    # start block
    time = datetime.utcnow()
    utils.log("\n{}\n...Starting DD_CAL: {}\n".format("-"*60, time), f=lf, verbose=verbose)
    cal_kwargs = copy.deepcopy(dict(list(algs['gen_cal'].items()) + list(algs['dd_cal'].items())))
    utils.log(json.dumps(cal_kwargs, indent=1) + '\n', f=lf, verbose=verbose)
    p = Dict2Obj(**cal_kwargs)

    # make a proper CASA spectral cube
    imname = '{}.{}'.format(datafile, p.model_ext)
    utils.log("...making a dummmy CASA spectral cube {}".format(imname), f=lf, verbose=verbose)
    cmd = casa + ['-c', "tclean(vis='{}', imagename='{}', niter=0, cell='{}arcsec', " \
                            "imsize={}, spw='', specmode='cube', start=0, width=300, stokes='{}')" \
                            "".format(datafile, imname, p.pxsize, p.imsize, p.stokes)]
    ecode = subprocess.call(cmd)

    # export to fits
    utils.log("...exporting to fits", f=lf, verbose=verbose)
    cmd = casa + ['-c', "exportfits('{}', '{}', overwrite={}, stokeslast=False)".format(imname+'.image', imname+".image.fits", overwrite)]
    ecode = subprocess.call(cmd)

    # erase all the original CASA files
    files = [f for f in glob.glob("{}*".format(imname)) if '.fits' not in f]
    for f in files:
        try:
            shutil.rmtree(f)
        except:
            os.remove(f)

    # make a spectral model of the imaged sources
    utils.log("...making spectral model of {}".format(p.inp_images), f=lf, verbose=verbose)
    inp_images = sorted(glob.glob(p.inp_images))
    assert len(inp_images) > 0, "nothing found under glob({})".format(p.inp_images)
    cmd = ["make_model_cube.py"] + inp_images \
        + ["--cubefile", imname + '.image.fits', '--sourcefile', p.sourcefile,
           '--outfname', imname + '.image.fits', "--makeplots", "--rb_Npix", p.rb_Npix,
           "--gp_ls", p.gp_ls, "--gp_nl", p.gp_nl,
           "--taper_alpha", p.taper_alpha, '--search_frac', p.search_frac]
    if overwrite:
        cmd += ['--overwrite']
    if p.fit_pl:
        cmd += ['--fit_pl']
    if p.fit_gp:
        cmd += ['--fit_gp']
    cmd += ['--exclude_sources'] + p.exclude_sources
    cmd = map(str, cmd)

    ecode = subprocess.call(cmd)

    # importfits
    utils.log("...importing from FITS", f=lf, verbose=verbose)
    cmd = casa +  ['-c', "importfits('{}', '{}', overwrite={})".format(imname+'.image.fits', imname+'.image', overwrite)]
    ecode = subprocess.call(cmd)
 
    # make a new flux model
    utils.log("...making new flux model for peeled visibilities, drawing parameters from gen_model", f=lf, verbose=verbose)
    utils.log(json.dumps(algs['gen_model'], indent=1) + '\n', f=lf, verbose=verbose)

    # First generate a clean_sources that excludes certain sources
    secondary_sourcefile = "{}_secondary.tab".format(os.path.splitext(p.sourcefile)[0])
    with open(secondary_sourcefile, "w") as f:
        f1 = open(p.sourcefile).readlines()
        f.write(''.join([l for i, l in enumerate(f1) if i-1 not in p.exclude_sources]))

    # Generate a New Model
    cal_kwargs['sourcefile'] = secondary_sourcefile
    model = gen_model(**dict(cal_kwargs + list(global_vars(varlist).items())))

    # uvsub model from corrected data
    utils.log("...uvsub CORRECTED - MODEL --> CORRECTED", f=lf, verbose=verbose)
    cmd = casa + ["-c", "uvsub('{}')".format(datafile)]
    ecode = subprocess.call(cmd)

    # split corrected
    split_datafile = "{}{}{}".format(os.path.splitext(datafile)[0], p.file_ext, os.path.splitext(datafile)[1])
    utils.log("...split CORRECTED to {}".format(split_datafile))
    cmd = casa + ["-c", "split('{}', '{}', datacolumn='corrected')".format(datafile, split_datafile)]
    ecode = subprocess.call(cmd)

    # Recalibrate
    utils.log("...recalibrating with peeled visibilities", f=lf, verbose=verbose)
    kwargs = dict(global_vars(list(varlist).items()) + list(cal_kwargs.items()))  # add order here is important
    dd_gtables = calibrate(**kwargs)

    # append new gaintables
    try:
        gaintables += dd_gtables
    except NameError:
        gaintables = dd_gtables

    # apply gaintables to datafile
    utils.log("...applying all gaintables \n\t{}\nto {}".format('\n\t'.join(gaintables), datafile), f=lf, verbose=verbose)
    cmd = casa + ['-c', '{}/sky_cal.py'.format(casa_scripts), '--msin', datafile, '--gaintables'] + gaintables
    ecode = subprocess.call(cmd)

    # end block
    time2 = datetime.utcnow()
    utils.log("...finished DD_CAL: {:d} sec elapsed".format(utils.get_elapsed_time(time, time2)), f=lf, verbose=verbose)

#-------------------------------------------------------------------------------
# Imaging
#-------------------------------------------------------------------------------
if params['dd_img']:
    # start block
    time = datetime.utcnow()
    utils.log("\n{}\n...Starting DD_IMG: {}\n".format("-"*60, time), f=lf, verbose=verbose)
    img_kwargs = copy.deepcopy(dict(list(algs['imaging'].items()) + list(algs['dd_img'].items())))
    utils.log(json.dumps(img_kwargs, indent=1) + '\n', f=lf, verbose=verbose)

    # Peform Imaging
    image(**img_kwargs)

    # end block
    time2 = datetime.utcnow()
    utils.log("...finished DD_IMG: {:d} sec elapsed".format(utils.get_elapsed_time(time, time2)), f=lf, verbose=verbose)

