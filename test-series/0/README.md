# Validation Step 0

* **Major Step Description**:
  * Simulate thermal noise, estimate PS with `hera_pspec`, calculate analytic and empirical uncertainty
  * Simulate sky-locked, white-spectrum Gaussian field and estimate PS with ``PSpec``.
* **Pipelines Tested**: ``hera_pspec.pspecdata``, ``hera_pspec.noise``, ``hera_pspec.grouping``, ``hera_sim.noise``.

## Minor Variations

### 0.0.x -- Thermal noise

* **Minor Step Description**: Generate thermal noise visibilities and calculate power spectra
* **Conditions for sucesss**: RMS of P(k) agrees with calculated analytic and empirical uncertainty up to sampling error of finite ensemble average
* **Expected First Trial**: 05/01/2019
* **First Trial**: 05/29/2019
* **Initial Pass Date**:

### 0.1.x -- Flat-k Gaussian P(k)

* **Minor Step Description**: Generate flat-k, sky-locked Gaussian field as per major step description.
* **Conditions for success**: Ensemble-averaged P(k) agrees with analytic expectation to within the predicted mean standard error. Potential inflation of variance from ``hera_pspec`` inaccuracies up to 1% allowable.
* **Expected First Trial**: 02/25/2019
* **First Trial**:
* **Initial Pass Date**:


### 0.2.x -- Non-Flat-k Gaussian P(k)

* **Minor Step Description**: Generate non-flat-k, sky-locked Gaussian field as per major step description.
* **Conditions for success**: Ensemble-averaged P(k) agrees with analytic expectation to within the predicted mean standard error. Potential inflation of variance from ``hera_pspec`` inaccuracies up to 1% allowable.
* **Expected First Trial**: 03/10/2019
* **First Trial**:
* **Initial Pass Date**:
