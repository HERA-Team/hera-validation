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

The data for the H1C sims reported in https://ui.adsabs.harvard.edu/abs/2021arXiv210409547A/abstract are available upon reasonable request.  For collaboration members, the paths on the NRAO machines are listed below.

There are two versions on the data, daily visibilities for 10 days (nominally 245 + [8098, 8099, 8101, 8102, 8103, 8106, 8107, 8108, 8110, 8111])  in `/lustre/aoc/projects/hera/Validation/test-4.0.0/data/visibilities/245*/`.  These were the starting points of the simulations.

Daily data which has been processed by some portion of the analysis pipeline is in `/lustre/aoc/projects/hera/Validation/test-4.0.0/pipeline/245????_*/`.  We assume most users are not interested in this step, but it is included below for completeness.

The final LST-binned data suitable for power spectrum analysis is in `/lustre/aoc/projects/hera/Validation/test-4.0.0/pipeline/LSTBIN/`.  

All three groups of data have different filename extensions.

For daily data we have `/lustre/aoc/projects/hera/Validation/test-4.0.0/data/visibilities/245*/`:

The files in each of these directories follow this naming convention: `zen.{jd_major}.{jd_minor}.{sky_component}.{state}.uvh5`.  For example:

`zen.2458098.32685.eor.true.uvh5
zen.2458098.32685.foregrounds.corrupt.uvh5
zen.2458098.32685.foregrounds.true.uvh5
zen.2458098.32685.sum.corrupt.uvh5
zen.2458098.32685.sum.true.uvh5
zen.2458098.32685.sum.uncal.ref_uncal.uvh5
zen.2458098.32685.sum.uncal.uvh5`

As an unpacked example: `foregrounds.true.uvh5` corresponds to foreground-only data that has not had any systematics (including noise) applied. Alternatively, `sum.corrupt.uvh5` has both foreground and EoR emission, and has been "corrupted" by the following systematics: thermal noise, bandpass gains, cable reflections, cross-coupling.

For the processed days in `/lustre/aoc/projects/hera/Validation/test-4.0.0/pipeline/245????_*/`:

`2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.abs.calfits
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
2458098_foregrounds/zen.2458098.49089.foregrounds.corrupt.calibrated.uvh5_image`

For the LST-binned data in `/lustre/aoc/projects/hera/Validation/test-4.0.0/pipeline/LSTBIN/`:

`sum`: foregrounds + eor, systematics applied, then calibrated out. many different `*.OCRSL*.uvh5` files that correspond to different operations applied to the data

`foregrounds`: foregrounds only, systematics appplied, then calibrated out. same notes as above.

`true_eor`: eor only, no systematics applied

`true_foregrounds`: foregrounds only, no systematics applied

`true_sum`: foregrounds + eor, no systematics applied

All of these files have been LST-binned. There are various stages of post-processing, pre-pspec analysis applied to the data, with that information stored in the section of the file name immediately preceding the .uvh5 extension.  In particular, 

O: omnical
C: abscal
R: RFI
S: smoothcal
L: LST binned
P: in-painted
X: cross-talk mitigated
T: time averaged
K: pseudo-Stokes

