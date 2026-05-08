# TSM heat-balance audit — Sep 2008 cool-T bias

**Date:** 2026-05-08
**Scope:** With BCs now properly transporting tracers (see
`design/bc_inflow_continuity_findings.md`), the persistent -2.68 °C
cool bias at the in-domain Salem cell is isolated to TSM. This memo
audits the runner's met-data preparation pipeline against the v3 TSM
module's expected conventions and identifies one strong candidate
explaining the bias.

## Conventions check

`08_run_coupled_v3_smoke.py` reads the KSLE hourly meteorology CSV and
populates the v3 registry forcings each master step (lines 977–981,
1126–1131):

| Registry variable | CSV column | Conversion | TSM expects |
|---|---|---|---|
| solar_radiation | `solar_W_m2` | passthrough | W/m² (used as net flux into water; no internal albedo) |
| air_temperature | `air_temp_C` | passthrough | °C |
| wind_speed | `wind_m_s` | **passthrough (no height correction)** | m/s **at 2 m above water surface** (Edinger 1974) |
| atmospheric_pressure | `pressure_mb` | passthrough | mb (= hPa) |
| atmospheric_vapor_pressure | `vapor_pressure_kPa × 10` | kPa→mb | mb |
| cloudiness | `cloud_frac` | passthrough | 0–1 |

All units are correct *except* wind. KSLE wind comes from an ASOS
anemometer at the WMO standard 10 m above ground; TSM's wind function
expects 2 m wind over the water surface.

## Wind-function sensitivity

TSM's wind function (`temperature.py:1197`):

```
f(W) = Ri × ( a/1e6  +  b/1e6 × W^c )
```

with **defaults `a = 0.3, b = 1.5, c = 3.0`**, calibrated against
Edinger, Brady & Geyer (1974). The cubic exponent makes the result
extremely sensitive to wind input height:

* Log-law profile over water (`z₀ ≈ 0.001 m`):
  `U₂ / U₁₀ = ln(2/z₀) / ln(10/z₀) ≈ 0.83`
* Combined ASOS-over-land → 2m-over-water (more conservative):
  `U₂_water / U₁₀_land ≈ 0.74–0.78`
* Effect on `f(W)` when 10 m wind is fed in unchanged:
  `(U₁₀/U₂)^c = (1/0.78)^3 ≈ 2.1`

So latent + sensible heat fluxes are ~**2× over-estimated** with the
current passthrough.

For Sep 2008 KSLE conditions (mean wind 1.4–2.1 m/s, dewpoint 8.8 °C
giving an air-water vapor-pressure deficit of ~7 mb):

* Latent heat at 2 m wind ~1.5 m/s:    ~50 W/m² out
* Latent heat at 10 m wind ~1.5 m/s used as 2 m: ~110 W/m² out
* Net forcing at 2 m wind:    +180 (solar) − 50 (latent) − 50 (LW) − 30 (sens) ≈ **+50 W/m² (warming)**
* Net forcing at 10 m wind:   +180 (solar) − 110 (latent) − 50 (LW) − 30 (sens) ≈ **−10 W/m² (slight cooling)**

That switches the equilibrium tendency from warming to cooling.
Given Sep 2008's air T mean of 18.1 °C and a healthy heat budget that
would equilibrate water near 17 °C (matching obs), the 60 W/m² extra
cooling drops the equilibrium by ~2 °C — close to the observed bias.

## Other forcings — secondary suspects

* **Solar:** `solar_W_m2` is plugged into `q_solar` as a passthrough
  (`temperature.py:643`). TSM does **not** apply an albedo
  internally; the input is treated as net shortwave at the surface.
  KSLE reports incoming shortwave (downward, top-of-canopy or
  near-surface). The albedo-not-corrected difference (~6% of
  ~180 W/m² mean = ~11 W/m²) goes the *wrong direction* (TSM gets
  ~11 W/m² *extra* solar gain), partially offsetting the wind error.
  Not a candidate for the cool bias; if anything, it *masks* it.

* **Cloud cover:** `cloud_frac` ∈ [0, 1] matches TSM convention. KSLE
  Sep 2008 mean = 0.06 (mostly clear), which produces strong net
  longwave loss (typical for clear nights). If the convention were
  flipped (0 = overcast, 1 = clear) we'd be feeding TSM "always
  overcast" → warmer equilibrium → no cool bias. The current
  convention is correct (verified: Sep 2008 ranges 0 to 1 with high
  mean values for cloudy days).

* **Vapor pressure:** the kPa→mb conversion is explicit and correct
  (line 983, `× 10.0`). Verified against the eair0 init at line 751.

* **Atmospheric pressure:** mb passthrough. KSLE reports station
  pressure already adjusted to sea level; sea-level vs. site pressure
  difference at Salem (elevation ~50 m) is ~6 mb, which barely
  affects the mixing-ratio formula. Not material.

* **Air temperature:** passthrough °C. No conversion issue.

## Recommended runner-side fix

Add a wind height correction at the case-study-script level (no
modules touch). Apply once at met-data ingest, then both the day-0
init (`wind_speed_init`) and the per-step update (`wind_arr`) carry
the corrected 2 m wind through to TSM. Suggested implementation:

```python
# At met-data load (after df_met_real is read):
WIND_HEIGHT_FACTOR = 0.78  # 10 m -> 2 m, log-law over water z0=0.001 m
if "wind_m_s" in df_met_real.columns:
    df_met_real["wind_m_s_2m"] = df_met_real["wind_m_s"] * WIND_HEIGHT_FACTOR
    # then read wind_m_s_2m in place of wind_m_s for both _met0 and _met_series
```

Or, more conservative, expose the factor as a CLI flag
`--wind-height-factor 0.78` (default 1.0 to preserve existing
behaviour, but document the recommended value for ASOS 10 m wind).

## Verification plan

1. Run the 15-day Sep 2008 simulation with `--wind-height-factor 0.78`
   and `--continuity-correction all_edges` (and `--no-reconstruct-newly-wet`,
   matching the alledges baseline), into a fresh
   `..._wind2m_alledges/` output dir.
2. Run Stage F validation against the existing 158k mesh + same
   obs-T BC drop-in.
3. Read residual at the Salem cell. Hypotheses:
   * Bias drops from -2.68 °C toward -1.0 °C: wind-height
     correction explains the lion's share of the bias.
   * Bias drops only modestly: there's an additional contributor
     (could be solar albedo, sediment temperature parameterisation,
     or a subtler issue in the heat-flux integration).
   * Bias becomes positive (warm bias): factor is too aggressive;
     try 0.83 (smoother water) or 0.74 (more conservative
     land-to-water transition).
4. Iterate factor or chase the next suspect.

## Notes for future TSM work

The `c = 3.0` cubic wind exponent is unusual. Edinger, Brady & Geyer
(1974) and many derived parameterisations use `c = 1.0` to `1.5`. The
v3 default of 3.0 amplifies any wind-height misuse and any wind
magnitude bias. Worth verifying against the original Edinger calibration
data before locking in. (This is a modules-side concern, deferred to
the v2→v3 merge.)

For HEC-RAS-style applications it's common to expose
`wind_a, wind_b, wind_c` as YAML-overridable parameters and document
the height assumption alongside. The runner currently uses defaults;
no override path exposed.

## Status

* Audit complete; wind-height passthrough identified as the strongest
  candidate.
* Runner-side fix is straightforward (~10 lines + 1 CLI flag).
* Verification run not yet executed — pending user decision on
  whether to implement the fix and re-run.
