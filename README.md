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


## Large-scale Validation Steps

0. P(k) white, on sky -> vis - P_meas(k)
1. GSM + pt src + EoR -> vis -> P_meas(k) (tests hera_pspec w / FRF?)
2. (FG matching abscal model + EoR) x per antenna gains (tests hera_cal + hera_pspec, first end-to-end)
3. FG + EoR + gains + cross-talk + RFI + … (end-to-end, systematic oriented)


**Additional Variations**
- P(k) non-white
- FG not matching abscal
- antenna-to-antenna beam variation
- beam real != beam model
- polarized sky
