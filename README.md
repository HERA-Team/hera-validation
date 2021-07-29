# HERA Validation Repository

**Archive of formal software pipeline validation tests**

This repository contains code and documentation for formal
tests of the HERA software pipeline. Tests are typically
performed and documented as Jupyter notebooks and are
archived to provide a long-standing account of the accuracy
of the pipeline as it evolves. Directory structures define the
broad kinds of tests performed.

## Mission Statement

The validation group seeks to validate the HERA data pipeline
software and algorithms by testing the specific software against
simulations where the expected output is well understood
theoretically.The group also helps to develop and define
increasingly sophisticated simulations on which to build an
end-to-end test and validation of the HERA pipeline.

## Structure of this repository

The validation effort seeks to verify the HERA software pipeline
through a number of well-defined steps of increasing complexity.
Each of these steps (called **major step**s or just **step**s in this
repository) reflects a broad validation concern or a specific 
element of the pipeline. For example, **step** 0 seeks to validate
just the ``hera_pspec`` software when given a well-known white-noise
P(k)-generated sky. 

Within each **step** exists the possibility of a set of variations 
(called **minor variation**s or just **variation**s in this repo). For 
example, variations for **step** 0 may be to generate flat-spectrum P(k)
and non-flat P(k). 

Finally, each combination of **step**-**variation** has the potential to incur
several staged tests or trials (we call them **trial**s in the repo). 

Importantly, failing **trial**s _will not be removed/overwritten_ in this
repo. Each formally-run trial is archived here for posterity. 

Thus the structure for this repo is as follows: Under the ``test-series``
directory, a number of directories labelled simply with their corresponding
**step** number are housed. Within each of these directories, each actual 
**trial** is presented as a notebook labelled ``test-<step>.<variation>.<trial>.ipynb``.

All **step**s, **variation**s and **trial**s are assigned increasing numerical
values. Generally, these values are increasing (from 0) in order of time/complexity.

In addition to the trial notebooks in these directories, each directory will
contain a ``README.md`` which lists the formal goals and conditions of each of
its **variation**s. 

Finally, each **variation** will be represented as a specific Github _project_,
in which the progress can be tracked and defined. Each project should receive 
a title which contains the **step**.**variation** identifier as well as a brief
description.

### Writing a validation test notebook

We have provided a template notebook which should serve as a starting
place for creating a validation notebook. The template is self-describing,
and has no intrinsic dependencies. All text in the notebook surrounded
by curly braces are meant to be replaced.

The template can be slightly improved/cleaned if you use jupyter notebook
extensions -- in particular the ToC and python-markdown extensions. The
first allows a cleaner way to build a table of contents (though one is
already included), and the latter allows using python variables in
markdown cells. This makes the writing of git hashes and versions simpler,
and means for example that the execution time/date can be written directly
into a markdown cell. 

## Project Plan
To create a simple tabulated version of the Project Plan, download the repo, save a
[personal access token](https://github.com/settings/tokens) to a file called `.pesonal-github-token`,
(ensure there is no trailing "\n" in the file)
and run `make_project_table.py` at the root directory. 
Note that you will need python 3.4+ and the `pygithub` code to run this script (`pip install pygithub`).
A semi-up-to-date version of this table is found at [project_table.md](./project_table.md).

## H1C Sims

The data for the H1C sims reported in the [H1C IDR2 Validation (Aguirre et al., 2021)](https://ui.adsabs.harvard.edu/abs/2021arXiv210409547A/abstract) are available upon 
reasonable request.  For collaboration members, the paths on the NRAO machines are listed below.

There are three main versions of the simulated data, for different levels of processing:

1. Daily raw visibilities (equivalent to raw HERA observations in terms of format/content)
2. Daily "processed" visibilities (processed by some portion of the analysis pipeline)
3. LST-binned data suitable for power spectrum analysis.

The simulated data has the following properties (see above linked paper for details):

* EoR is a Gaussian Random Field with power-law power spectrum
* Foregrounds are GLEAM + eGSM
* 10 full days were simulated (nominally 245 + [8098, 8099, 8101, 8102, 8103, 8106, 8107, 8108, 8110, 8111])
* Visibility Simulator was `RIMEz`
* Primary Beam was Fagnoni+2019 HERA beam model
* Simulations consist of a strict subset of H1C antennas / baselines
* Frequency range/resolution is the same as H1C data (i.e. 100-200 MHz in 1024 channels)
* Systematics applied to the data (where appropriate) are:
  * Thermal Noise: based on constant receiver noise plus simulated autos
  * Bandpass Gains: randomly perturbed around H1C bandpass measurements
  * Cable Reflections
  * Cross-coupling

Below, we specify for each of the three data versions where to find the data.

### Daily Raw Visibilities

These reside in `/lustre/aoc/projects/hera/Validation/test-4.0.0/data/visibilities/245*/` (one folder per day of mock-observation).
The files in each of these directories follow this naming convention: `zen.{jd_major}.{jd_minor}.{sky_component}.{state}.uvh5`.  
The `{sky_component}` indicates which sky models are present in the data, and can be one of

* `eor`: The EoR.
* `foregrounds`: The Foregrounds (both GLEAM and eGSM)
* `sum`: Both EoR and Foregrounds

The `{state}` indicates whether the state of the data in terms of systematics, and can be one of

* `true`: No systematics included (including noise)
* `corrupt`: All systematics included (see above list)
* `uncal`: only bandpass gains included.

For example, the command `ls zen.2458098.32685.*` yields:

```
zen.2458098.32685.eor.true.uvh5
zen.2458098.32685.foregrounds.corrupt.uvh5
zen.2458098.32685.foregrounds.true.uvh5
zen.2458098.32685.sum.corrupt.uvh5
zen.2458098.32685.sum.true.uvh5
zen.2458098.32685.sum.uncal.ref_uncal.uvh5
zen.2458098.32685.sum.uncal.uvh5
```

### Daily Processed Visibilities
We assume most users are not interested in this step, but it is included for completeness.

These partially-processed visibilities are in 
`/lustre/aoc/projects/hera/Validation/test-4.0.0/pipeline/245{jd}_{kind}/`,
where `{kind}` represents which data was input to the processing, and can be one of

* `foregrounds`: Just foregrounds (no EoR) with all systematics
* `sum`: FG+EoR with all systematics
* `sum_uncal`: FG+EoR with just bandpass gains but not other systematics

Other combinations of mock observations (eg. `eor.corrupt`) were not processed.

Each of these folders contains a number of files. Each file corresponds to a stage/product
of processing, such as firstcal, abscal, smoothcal etc. For example, the files in 
`/lustre/aoc/projects/hera/Validation/test-4.0.0/pipeline/2458098_foregrounds/` are:

```
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.abs.calfits
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.autos.uvh5
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.first.calfits
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.firstcal_metrics.hdf5
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.flagged_abs.calfits
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.flagged_abs_vis.uvh5
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.noise_std.uvh5
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.omni.calfits
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.omni_vis.uvh5
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.smooth_abs.calfits
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.smooth_abs_vis.uvh5
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.uvh5
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.calibrated.ms
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.calibrated.uvh5_image
```

The fully-processed files end with `.smooth_abs_vis.uvh5`.

### LST-Binned Data

The final LST-binned data suitable for power spectrum analysis is in 
`/lustre/aoc/projects/hera/Validation/test-4.0.0/pipeline/LSTBIN/`.  
Within this directory, the LST-binned data for different combinations
of input models and processing is kept in different directories, namely:

* `sum/`: FG+EoR, systematics applied, then calibrated out. 
* `foregrounds/`: FG-only (no EoR), systematics appplied, then calibrated out.
* `true_eor/`: EoR-only, no systematics applied
* `true_foregrounds/`: FG-only, no systematics applied
* `true_sum/`: FG+EoR, no systematics applied

In particular, no EoR-only with systematics was produced.

In each directory there are a great number of files. The most important files (i.e the LST-binned
and pre-processed visibilities) have the filename convention `zen.grp1.of1.LST.{lst}.HH.{processing_tags}.uvh5`.
Here the `lst` is a floating point number giving the LST in radians. The `{processing_tags}` are 
a group of single upper-case characters that indicate which processing steps have been applied. They are:

* O: omnical
* C: abscal
* R: RFI
* S: smoothcal
* L: LST binned
* P: in-painted
* X: cross-talk mitigated
* T: time averaged
* K: pseudo-Stokes

For example, doing `ls LSTBIN/foregrounds/zen.grp1.of1.LST.*.HH.O*.uvh5` gives

```
LSTBIN/foregrounds/zen.grp1.of1.LST.0.28190.HH.OCRSLP.uvh5     LSTBIN/foregrounds/zen.grp1.of1.LST.1.03362.HH.OCRSL.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.28190.HH.OCRSL.uvh5      LSTBIN/foregrounds/zen.grp1.of1.LST.1.03441.HH.OCRSLPX.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.28268.HH.OCRSLPX.uvh5    LSTBIN/foregrounds/zen.grp1.of1.LST.1.12759.HH.OCRSLP.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.28973.HH.OCRSLPXTK.uvh5  LSTBIN/foregrounds/zen.grp1.of1.LST.1.12759.HH.OCRSL.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.28973.HH.OCRSLPXT.uvh5   LSTBIN/foregrounds/zen.grp1.of1.LST.1.12837.HH.OCRSLPX.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.37586.HH.OCRSLP.uvh5     LSTBIN/foregrounds/zen.grp1.of1.LST.1.22156.HH.OCRSLP.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.37586.HH.OCRSL.uvh5      LSTBIN/foregrounds/zen.grp1.of1.LST.1.22156.HH.OCRSL.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.37665.HH.OCRSLPX.uvh5    LSTBIN/foregrounds/zen.grp1.of1.LST.1.22234.HH.OCRSLPX.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.46983.HH.OCRSLP.uvh5     LSTBIN/foregrounds/zen.grp1.of1.LST.1.22939.HH.OCRSLPXTK.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.46983.HH.OCRSL.uvh5      LSTBIN/foregrounds/zen.grp1.of1.LST.1.22939.HH.OCRSLPXT.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.47061.HH.OCRSLPX.uvh5    LSTBIN/foregrounds/zen.grp1.of1.LST.1.31552.HH.OCRSLP.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.56380.HH.OCRSLP.uvh5     LSTBIN/foregrounds/zen.grp1.of1.LST.1.31552.HH.OCRSL.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.56380.HH.OCRSL.uvh5      LSTBIN/foregrounds/zen.grp1.of1.LST.1.31631.HH.OCRSLPX.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.56458.HH.OCRSLPX.uvh5    LSTBIN/foregrounds/zen.grp1.of1.LST.1.40949.HH.OCRSLP.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.60295.HH.OCRSLPXTK.uvh5  LSTBIN/foregrounds/zen.grp1.of1.LST.1.40949.HH.OCRSL.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.60295.HH.OCRSLPXT.uvh5   LSTBIN/foregrounds/zen.grp1.of1.LST.1.41027.HH.OCRSLPX.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.65776.HH.OCRSLP.uvh5     LSTBIN/foregrounds/zen.grp1.of1.LST.1.50345.HH.OCRSLP.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.65776.HH.OCRSL.uvh5      LSTBIN/foregrounds/zen.grp1.of1.LST.1.50345.HH.OCRSL.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.65854.HH.OCRSLPX.uvh5    LSTBIN/foregrounds/zen.grp1.of1.LST.1.50424.HH.OCRSLPX.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.75173.HH.OCRSLP.uvh5     LSTBIN/foregrounds/zen.grp1.of1.LST.1.54261.HH.OCRSLPXTK.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.75173.HH.OCRSL.uvh5      LSTBIN/foregrounds/zen.grp1.of1.LST.1.54261.HH.OCRSLPXT.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.75251.HH.OCRSLPX.uvh5    LSTBIN/foregrounds/zen.grp1.of1.LST.1.59742.HH.OCRSLP.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.84569.HH.OCRSLP.uvh5     LSTBIN/foregrounds/zen.grp1.of1.LST.1.59742.HH.OCRSL.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.84569.HH.OCRSL.uvh5      LSTBIN/foregrounds/zen.grp1.of1.LST.1.59820.HH.OCRSLPX.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.84648.HH.OCRSLPX.uvh5    LSTBIN/foregrounds/zen.grp1.of1.LST.1.69139.HH.OCRSLP.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.91617.HH.OCRSLPXTK.uvh5  LSTBIN/foregrounds/zen.grp1.of1.LST.1.69139.HH.OCRSL.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.91617.HH.OCRSLPXT.uvh5   LSTBIN/foregrounds/zen.grp1.of1.LST.1.69217.HH.OCRSLPX.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.93966.HH.OCRSLP.uvh5     LSTBIN/foregrounds/zen.grp1.of1.LST.1.78535.HH.OCRSLP.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.93966.HH.OCRSL.uvh5      LSTBIN/foregrounds/zen.grp1.of1.LST.1.78535.HH.OCRSL.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.0.94044.HH.OCRSLPX.uvh5    LSTBIN/foregrounds/zen.grp1.of1.LST.1.78613.HH.OCRSLPX.uvh5
LSTBIN/foregrounds/zen.grp1.of1.LST.1.03362.HH.OCRSLP.uvh5
```
