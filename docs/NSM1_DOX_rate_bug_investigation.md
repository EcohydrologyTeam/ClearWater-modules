# NSM1 DOX Rate Bug Investigation

**Date:** 2026-02-20
**Reported symptom:** DOX rate in NSM1 suddenly massively increases when T > 20 deg C.

## Root Cause

The default value `SOD_theta = 999` in `src/clearwater_modules/nsm1/constants.py` (line 323, `DEFAULT_GLOBALVARS`) is a placeholder sentinel value. When used in the Arrhenius temperature correction formula, it causes catastrophic numerical blow-up at temperatures above 20 deg C.

### The Arrhenius correction formula

Located in `src/clearwater_modules/nsm1/processes.py` (line 35):

```python
rc20 * theta ** (TwaterC - 20.0)
```

At 20 deg C, the exponent is zero, so `theta^0 = 1` regardless of the theta value. Above 20 deg C, the exponent becomes positive, and `999^(T-20)` explodes:

| Temperature (deg C) | Exponent | 999^exponent | SOD_tc (SOD_20 * 999^exp) |
|---|---|---|---|
| 19 | -1 | 0.001 | ~1 |
| 20 | 0 | 1 | 999 |
| 20.5 | 0.5 | 31.6 | ~31,600 |
| 21 | 1 | 999 | ~998,000 |
| 22 | 2 | 998,001 | ~10^9 |
| 25 | 5 | ~10^15 | ~10^18 |

The SOD temperature-corrected rate feeds directly into the DOX budget via `DOX_SOD = SOD_tc / depth` (`processes.py` line 3092), which is subtracted in the `dDOXdt` equation (`processes.py` line 3119). This astronomical SOD overwhelms all other terms in the DOX budget.

## All sentinel values in DEFAULT_GLOBALVARS

The following parameters in `DEFAULT_GLOBALVARS` (`constants.py` lines 317-351) use the value `999` as a placeholder:

| Parameter | Default | Description |
|---|---|---|
| `vsop` | 999 | Organic phosphorus settling velocity (m/d) |
| `vs` | 999 | Settling velocity (m/d) |
| `SOD_20` | 999 | Sediment oxygen demand at 20 deg C (mg-O2/m2/d) |
| `SOD_theta` | 999 | Arrhenius coefficient for SOD |
| `kaw_20_user` | 999 | User-defined wind reaeration velocity at 20 deg C (m/d) |
| `kah_20_user` | 999 | User-defined hydraulic reaeration rate at 20 deg C (1/d) |

Of these, `SOD_theta` is by far the most damaging because it appears as the **base of the exponential** in the Arrhenius correction. The others create unrealistically high base rates but do not compound exponentially with temperature (their theta values are reasonable: 1.024 or 1.047).

### Additional unrealistic defaults

| Parameter | Default | Physically realistic range |
|---|---|---|
| `pressure_mb` | 2026.5 | ~870-1013 mb (1013.25 = sea level) |
| `slope` | 2 | ~0.0001-0.01 for most channels |
| `shear_velocity` | 4 m/s | ~0.01-0.1 m/s typically |

## Why tests don't catch this

The test suite (`tests/test_11_nsm_DOX_calculations.py`) overrides all sentinel values with scenario-specific values. For example:

```python
SOD_20 = 0.5,       # instead of 999
SOD_theta = 1.047,   # instead of 999
kaw_20_user = 0,     # instead of 999
kah_20_user = 1,     # instead of 999
pressure_mb = 1013.25  # instead of 2026.5
```

The test only runs at T=25 and T=15 deg C, and never at T=20 to compare behavior across the 20 deg C threshold. Since tests always provide their own parameters, the broken defaults are never exercised.

## Recommended fix

The sentinel values are intentional -- they signal that users must provide site-specific values. However, nothing in the code validates that these have been replaced. The recommended fix is to add **input validation at model initialization** that raises an error or warning when sentinel values (999) are detected in parameters that feed into Arrhenius corrections or other exponential computations. Simply replacing the sentinels with "typical" values would mask configuration errors in real applications.
