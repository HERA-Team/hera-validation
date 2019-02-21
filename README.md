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

## Plan for Validation

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
