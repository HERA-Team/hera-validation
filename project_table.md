## Project Plan
Here follows a formal plan for envisaged tests, automatically generated from our Projects. 
This list is *not final*: it may be extended or 
modified at any time. However, it provides a reasonable snapshot of the current status and plans of 
the validation team.

What the symbols mean:

Sym. | Meaning | Sym. | Meaning | Sym.  |Meaning
-------| ----- | ----- | ---- | ----- | ------
:egg:      | Idea in gestation  | :hammer:   | Currently being worked on | :thinking: | PR submitted... 
:heavy_check_mark: | Passed     | :x:        | Failed | :watch:    | Not required for H1C IDR2

### [Step -1](https://api.github.com/projects/3274950)
Use formal pyuvsim [reference simulations](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/master/reference_simulations) as strict comparison points with every visibility simulator used in any other step, using [this template](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/pull/211).

Status     | #    | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| -----|--------------|-----------------|---------------------|----------|
:hammer:  |  [-1](https://api.github.com/repos/HERA-Team/hera-validation/issues/31):  notebook  |   |  |  | [@steven-murray](https://api.github.com/users/steven-murray), [@piyanatk](https://api.github.com/users/piyanatk)  |
:hammer:  |  [-1.0](https://api.github.com/repos/HERA-Team/hera-validation/issues/25)  | healvis, pyuvsim  | gsm |  |   |
:hammer:  |  [-1.1](https://api.github.com/repos/HERA-Team/hera-validation/issues/35):  Validation of RIMEz  | pyuvsim, rimez  | gsm |  | [@piyanatk](https://api.github.com/users/piyanatk), [@zacharymartinot](https://api.github.com/users/zacharymartinot)  |
:hammer:  |  [-1.2](https://api.github.com/repos/HERA-Team/hera-validation/issues/36):  PRISim Validation  | prisim, pyuvsim  | gleam, gsm |  | [@piyanatk](https://api.github.com/users/piyanatk), [@nithyanandan](https://api.github.com/users/nithyanandan)  |
:hammer:  |  [-1.3](https://api.github.com/repos/HERA-Team/hera-validation/issues/37):  vis_cpu validation  | pyuvsim, viscpu  | gleam |  | [@steven-murray](https://api.github.com/users/steven-murray), [@piyanatk](https://api.github.com/users/piyanatk), [@Jackmastr](https://api.github.com/users/Jackmastr)  |


### [Step 0](https://api.github.com/projects/3274969)
Test `hera_pspec`'s ability to reproduce known power spectra from EoR-only simulations of visibilities (which are sky-locked) and noise visibilities. Noise models are both white with frequency and time, and following a fiducial sky model. Noise is taken from `hera_sim` and added in visibilities to the various simulators. This should (eventually) be able to deal with different amounts of coherent and incoherent averaging.

Status     | #    | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| -----|--------------|-----------------|---------------------|----------|
:heavy_check_mark:  |  [0.0](https://api.github.com/repos/HERA-Team/hera-validation/issues/5):  Noise Pspec  | hera_sim  | noise | pspec | [@nkern](https://api.github.com/users/nkern)  |
:heavy_check_mark:  |  [0.1](https://api.github.com/repos/HERA-Team/hera-validation/issues/7)  | healvis  | flat | pspec | [@r-pascua](https://api.github.com/users/r-pascua)  |
:egg:  |  [0.2](https://api.github.com/repos/HERA-Team/hera-validation/issues/23)  | rimez  | powerlaw | pspec |   |
:hammer:  |  [0.3](https://api.github.com/repos/HERA-Team/hera-validation/issues/38):  P(k) from 21cmFAST  | prisim  | 21cmfast | pspec | [@nithyanandan](https://api.github.com/users/nithyanandan)  |
:hammer:  |  [0.4](https://api.github.com/repos/HERA-Team/hera-validation/issues/27)  | hera_sim, prisim  | flat, noise | pspec | [@acliu](https://api.github.com/users/acliu), [@taylordb](https://api.github.com/users/taylordb)  |
:hammer:  |  [0.5](https://api.github.com/repos/HERA-Team/hera-validation/issues/39):  Sharp-Feature P(k)  | rimez  | feature | pspec | [@zacharymartinot](https://api.github.com/users/zacharymartinot), [@JianrongTan](https://api.github.com/users/JianrongTan)  |
:hammer:  |  [0.6](https://api.github.com/repos/HERA-Team/hera-validation/issues/26):  Crazy Re-weighting  |   | flat | pspec | [@acliu](https://api.github.com/users/acliu), [@taylordb](https://api.github.com/users/taylordb)  |


### [Step 1](https://api.github.com/projects/3274994)
Test `hera_pspec`'s ability to recover EoR P(k) from visibility simulations including unpolarized foregrounds and noise. This includes tests with different amounts of coherent and incoherent averaging. Error bars should correctly be predicted from noise and signal levels.

Status     | #    | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| -----|--------------|-----------------|---------------------|----------|
:hammer:  |  [1.2](https://api.github.com/repos/HERA-Team/hera-validation/issues/40):  GSM and GLEAM  | healvis, rimez  | powerlaw, gleam, gsm | pspec | [@r-pascua](https://api.github.com/users/r-pascua)  |


### [Step 2](https://api.github.com/projects/3275007)
Test effect of `hera_cal` on recovered P(k) (i.e. `hera_pspec`), accummulated piecewise from bits of the analysis flowchart. The underlying assumptions of the calibration (ideal antenna positions, identical beams, smooth antenna-based gains) are respected.

Status     | #    | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| -----|--------------|-----------------|---------------------|----------|
:heavy_check_mark:  |  [2.0](https://api.github.com/repos/HERA-Team/hera-validation/issues/4):  `redcal` and `abscal`  | healvis  | gleam | abscal, redcal | [@jaguirre](https://api.github.com/users/jaguirre), [@jsdillon](https://api.github.com/users/jsdillon)  |
:hammer:  |  [2.1](https://api.github.com/repos/HERA-Team/hera-validation/issues/16)  | rimez  | powerlaw, gleam, gsm, gains | abscal, pspec, redcal, smoothcal |   |
:hammer:  |  [2.2](https://api.github.com/repos/HERA-Team/hera-validation/issues/28):  Validation of reference model construction  | rimez  | powerlaw, gains, noise | abscal, casa, redcal | [@TashaleeB](https://api.github.com/users/TashaleeB)  |


### [Step 3](https://api.github.com/projects/3275013)
Test the effects of realistic RFI, and data-driven flags, as well as other post-analysis systematics on various parts of the pipeline.

Status     | #    | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| -----|--------------|-----------------|---------------------|----------|
:hammer:  |  [3.0](https://api.github.com/repos/HERA-Team/hera-validation/issues/20)  | hera_sim, rimez  | powerlaw, noise | pspec | [@zacharymartinot](https://api.github.com/users/zacharymartinot)  |
:hammer:  |  [3.1](https://api.github.com/repos/HERA-Team/hera-validation/issues/21)  | rimez  | powerlaw, gsm | pspec, smoothcal | [@zacharymartinot](https://api.github.com/users/zacharymartinot), [@jburba](https://api.github.com/users/jburba)  |
:hammer:  |  [3.2](https://api.github.com/repos/HERA-Team/hera-validation/issues/22)  | hera_sim, rimez  | powerlaw, gsm, gains, xtalk | pspec, smoothcal | [@nkern](https://api.github.com/users/nkern), [@r-pascua](https://api.github.com/users/r-pascua)  |
:hammer:  |  [3.3](https://api.github.com/repos/HERA-Team/hera-validation/issues/41):  Test of xRFI  | hera_sim, rimez  | powerlaw, gleam, gsm, gains, rfi | pspec, redcal, smoothcal, xrfi | [@steven-murray](https://api.github.com/users/steven-murray), [@lwhitler](https://api.github.com/users/lwhitler)  |


### [Step 4](https://api.github.com/projects/3275024)
Test most or all components of the full pipeline (including redcal, abscal, xRFI, smoothcal, hera_pspec).

Status     | #    | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| -----|--------------|-----------------|---------------------|----------|
:hammer:  |  [4.0](https://api.github.com/repos/HERA-Team/hera-validation/issues/42):  Basic End-to-End  | hera_sim, rimez  | powerlaw, gleam, gsm, gains, reflections, rfi, xtalk | abscal, casa, pspec, redcal, smoothcal, xrfi | [@steven-murray](https://api.github.com/users/steven-murray), [@r-pascua](https://api.github.com/users/r-pascua)  |


