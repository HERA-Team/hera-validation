This document specifies where various data products used in the validation of
the H1C IDR3 analysis can be found.

There are three base sets of raw simulation data. These can be thought of as
noiseless, perfectly-calibrated data sets.

diffuse: 2008 GSM + confused point sources
ptsrc: mock point source catalog + brights (e.g. Fornax, Pictor, etc.)
eor: realization of a reionization signal with 1/k^2 power spectrum

These data sets can be found in the following directories:
/lustre/aoc/projects/hera/Validation/H1C_IDR3/chunked_data/{sky_cmp}

Here, {sky_cmp} is any choice of ["diffuse", "ptsrc", "eor"].

The data products have the following file prototype:
zen.LST.{lst_major}.{lst_minor}.{sky_cmp}.uvh5

{lst_major}.{lst_minor} specifies the LST of the first integration in the file.
The filenames are chosen so that the noted LSTs cover 2\pi radians, starting from
roughly 4.7 radians. This coincides with the LST branch cut chosen for this
analysis (see HERA Memo 107, section 2.3), and arranging the data by increasing
LST ensures that the logic for time-interpolation remains simple.

The raw simulation data covers 24 hours of LST with an integration time of 5
seconds (for a total of 17280 integrations), spans 100 MHz to 200 MHz over 1024
frequency channels, and uses an array with 48 elements.

The array layout is stored in a pyuvsim-style csv, located here:
/lustre/aoc/projects/hera/Validation/H1C_IDR3/antenna_select/array_layout.csv

Daily data products are stored in the following directories:
/lustre/aoc/projects/hera/Validation/test-4.1.0/makeflow/{jd}

Configuration files for simulating systematic effects can be found here:
/lustre/aoc/projects/hera/Validation/H1C_IDR3/configs/{jd}.yaml

The following systematics are applied to the simulated data:
bandpass gains: randomized each day, constant in time for any given day
reflections: phase, amplitude, jitter randomized each day
mutual coupling: same cable properties as reflections, based on HERA Memo 104
thermal noise: fixed receiver temperature, all draws independent

Note that for reflection-like systematics, the cable delays are chosen once,
and they are slightly varied between each epoch. The jitter is randomized each
day to mimic daily-changing temperature distributions on site.
