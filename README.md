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
Here follows a formal plan for envisaged tests. This list is *not final*: it may be extended or modified at any time. However, it provides a reasonable snapshot of the current status and plans of the validation team.

What the symbols mean:

Sym. | Meaning | Sym. | Meaning | Sym.  |Meaning
-------| ----- | ----- | ---- | ----- | ------
:egg:      | Idea in gestation  | :hammer:   | Currently being worked on | :thinking: | PR submitted... 
:heavy_check_mark: | Passed     | :x:        | Failed | :watch:    | Outside current scope

### [Step -1](test-series/-1/): Compare visibility simulators with `pyuvsim`.  
Use formal `pyuvsim` [reference simulations](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/master/reference_simulations)
  as strict comparison points with every visibility simulator used in any other step, using [this template]( https://github.com/RadioAstronomySoftwareGroup/pyuvsim/pull/211).

Status     | #   | Description | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| ----|-------------|--------------|-----------------|---------------------|----------|
:hammer:   | -1.0| `healvis`   | `healvis`    | `GSM` (?)       | None                | [@alanman][al] |
:hammer:   | -1.1| `RIMEz`     | `RIMEz`      | `GSM` (?)       | None                | [@zacharymartinot][zm] | 
:hammer:   | -1.2| `PRISim`    | `PRISim`     | GLEAM           | None                | [@nithyanandan][nt]   |
:egg:      | -1.3| `vis_cpu`   | `vis_cpu`    | GLEAM           | None                | [@steven-murray][sgm] |
  
### [Step 0](test-series/0/): Test `hera_pspec` directly, without foregrounds.
Test `hera_pspec`'s ability to reproduce known power spectra from EoR-only simulations of visibilities (which are sky-locked) and noise visibilities. Noise models are both white with frequency and time, and following a fiducial sky model.  Noise is taken from hera_sim and added in visibilities to the various simulators.  This should (eventually) be able to deal with different amounts of coherent and incoherent averaging.

Status     | #   | Description | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| ----|-------------|--------------|-----------------|---------------------|----------|
:thinking: | [0.0](https://github.com/HERA-Team/hera-validation/pull/5) | Uncorrelated _k_- and time-dep noise | `hera_sim` | `noise`     | `hera_pspec`        | [@nkern][nk]   |
:heavy_check_mark: | 0.1 | Flat P(k) | `healvis` | EoR | `hera_pspec` | [@r-pascua][rp] |
:hammer:   | 0.2 | Power-law P(k) | `RIMEz` | EoR | `hera_pspec` | [@zacharymartinot][zm] |
:egg:      | 0.3 | Sharp-feature P(k) | `RIMEz` | EoR | `hera_pspec` | [@zacharymartinot][zm] |
:egg:      | 0.4 | P(k) from 21cmFAST | `PRISim` | EoR | `hera_pspec` | [@nithyanandan][nt] |

### [Step 1](test-series/1/): Test `hera_pspec` directly, with foregrounds.
Test `hera_pspec`'s ability to recover EoR P(k) from visibility simulations including unpolarized foregrounds and noise. This includes tests with different amounts of coherent and incoherent averaging.  Error bars should correctly be predicted from noise and signal levels.  

Status     | #   | Description | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| ----|-------------|--------------|-----------------|---------------------|----------|
:thinking: | 1.0 | Power-law P(k) + Diffuse (GSM) | `RIMEz` | EoR, GSM | `hera_pspec` | [@zacharymartinot][zm] |
:egg:      | 1.1 | Power-law P(k) + point sources | `RIMEz` | EoR, GLEAM or Random Pt.Srcs. | `hera_pspec` | [@zacharymartinot][zm] |
:egg:      | 1.2 | Flat P(k) + GSM + point sources + noise | `healvis` (?), `hera_sim` | EoR, GSM, GLEAM, `noise` | `hera_pspec` | [@alanman][al], [@r-pascua][rp] |  

### [Step 2](test-series/2/): Test `hera_cal`'s effect on recovered P(k)
Test effect of `hera_cal` on recovered P(k) (i.e. `hera_pspec`), accummulated piecewise from bits of the analysis flowchart. The underlying assumptions of the calibration (ideal antenna positions, identical beams, smooth antenna-based gains) are respected.  

Status     | #   | Description | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| ----|-------------|--------------|-----------------|---------------------|----------|
:thinking: | [2.0](https://github.com/HERA-Team/hera-validation/pull/4) | Recover known gains | `healvis` | Pt.Srcs. |  `redcal`, `abscal` | [@jsdillon][jsd], [@jaguirre][ja] |
:egg:      | 2.1 | Recover flat P(k) with known gains | `healvis` |  EoR, Pt.Srcs. | `redcal`, `abscal`, `smoothcal`, `hera_pspec` | [@alanman][al], [@jaguirre][ja] |
:egg:      | 2.2 | Validation of reference model construction. | * | * | * | * |

### [Step 3](test-series/3/): Test effects of RFI (and `xRFI`)
Test the effects of both realistic RFI, and data-driven *flags* on various parts of the pipeline.

Status     | #   | Description | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| ----|-------------|--------------|-----------------|---------------------|----------|
:hammer:   | 3.0 | Freq-dep noise (according to flagged channels) direct to `pspec` | `healvis`, `hera_sim` | EoR, `noise` | `hera_pspec` | [@zacharymartinot][zm] |
:egg:      | 3.1 | Apply *data* RFI flags  | `healvis`(?) | EoR, GSM | `smoothcal`, `pspec` | ? | 
:egg:      | 3.2 | Apply *data* RFI flags w/systematics | `healvis`(?), `hera_sim` | EoR, GSM, `sigchain.gen_gains`, `sigchain.xtalk` | `smoothcal`, `pspec` | ? | 
:egg:      | 3.3 | Simulated RFI and `xRFI` | `healvis`, `hera_sim` | GSM, `rfi` | `xRFI` | [@steven-murray][sgm] |

### [Step 4](test-series/4): Test full end-to-end pipeline at modest realism
Test most or all components of the full pipeline (including `redcal`, `abscal`, `xRFI`, `smoothcal`, `hera_pspec`).

Status     | #   | Description | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| ----|-------------|--------------|-----------------|---------------------|----------|
:hammer:   | 4.0 |   | `healvis`, `hera_sim` | EoR, GSM, GLEAM, `rfi` | `abscal`, `redcal`, `xRFI`, `smoothcal`, `pspec` | [@steven-murray][sgm], [@r-pascua][rp] |
:egg:      | 4.1 | Test LST-binning and "fringe rate filtering" (time averaging).
:watch:    | 4.2 | Non-ideal antenna positions
:watch:    | 4.3 | Antenna-to-antenna beam variation
:watch:    | 4.4 | Beam real != beam model
:watch:    | 4.5 | Polarized sky

[rp]: https://github.com/r-pascua
[nk]: https://github.com/nkern
[sgm]: https://github.com/steven-murray
[ja]: https://github.com/jaguirre
[zm]: https://github.com/zacharymartinot
[jsd]: https://github.com/jsdillon
[al]: https://github.com/aelanman
[nt]: https://github.com/nithyanandan
