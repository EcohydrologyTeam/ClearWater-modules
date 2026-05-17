"""v3 NSM1 benthic algae parameter defaults.

Consumed by: ``BenthicAlgae`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``BalgaeStaticVariables`` /
``DEFAULT_BALGAE``.

Corrections applied:

* ``BWa`` corrected from ``3500.0`` to ``1000.0`` ug-Chla per stoichiometric
  unit (Phase 9.E follow-up). Pre-9.E v3 inherited v1's ``3500`` (Fortran
  uses ``5000``); the derived ``rab = BWa / BWd`` Chla:DW ratio for v1/v3
  was 35 mg-Chla/g-DW (3.5x WASP7 canonical) and Fortran's was 50 mg/g
  (5x WASP7). Phase 9.E harmonized v3 to WASP7's documented canonical
  benthic Chla:DW = 10 mg-Chla/g-DW (computed from WASP7 Benthic Algae
  User's Guide Table 1: Chla:C = 0.025 mg-Chla/mg-C, DW:C = 2.5 mg-DW/mg-C
  -> Chla:DW = 0.01 mg/mg = 10 mg/g). With ``BWd = 100``, this requires
  ``BWa = 1000``, giving ``rab = 10`` mg-Chla/g-DW. The corrected value
  also matches NSM1's own floating-algae ratio (``AWa/AWd = 1000/100 =
  10``), bringing benthic and floating algae onto the same Chla:DW
  basis -- consistent with the WASP7 convention where benthic and
  floating algae share the same stoichiometry. See
  ``parameter_defaults_corrections.md`` Section 1.13.
"""

DEFAULTS: dict[str, float | int | bool] = {
    # BWd/BWc/BWn/BWp: raw stoichiometric weights "per stoichiometric
    # unit". Mass-fraction ratios are the *quotients* rdb = BWd/BWd = 1,
    # rcb = BWc/BWd, rnb = BWn/BWd, rpb = BWp/BWd (mg-X per mg-D);
    # see ``benthic_algae.py:336`` and ``carbon.py:403``. With BWd = 100
    # these give rcb = 0.4 mg-C/mg-D = 40% C, rnb = 0.072 mg-N/mg-D =
    # 7.2% N, rpb = 0.01 mg-P/mg-D = 1% P (typical algal mass fractions).
    'BWd': 100.0,                       # mg-D per stoichiometric unit (denominator for rcb/rnb/rpb; produces 100% dry-weight mass fraction by definition)
    'BWc': 40.0,                        # mg-C per stoichiometric unit (rcb = BWc/BWd = 0.4 mg-C/mg-D)
    'BWn': 7.2,                         # mg-N per stoichiometric unit (rnb = BWn/BWd = 0.072 mg-N/mg-D)
    'BWp': 1.0,                         # mg-P per stoichiometric unit (rpb = BWp/BWd = 0.01 mg-P/mg-D)
    'BWa': 1000.0,                      # ug-Chla per stoichiometric unit; gives rab = BWa/BWd = 10 mg-Chla/g-DW (matches WASP7 canonical and NSM1 floating-algae AWa/AWd; was 3500 in v1/v3 pre-Phase-9.E and 5000 in Fortran, both above canonical)
    'KLb': 10.0,                        # W/m^2 PAR; benthic light limitation half-saturation (PAR-scale)
    'Fr_PAR': 0.47,                     # unitless; PAR fraction of total broadband shortwave. NSM1-SCI-A3 (spec B1): restores v1 PAR=q_solar*Fr_PAR for benthic algae; matches floating-algae + pathogen Fr_PAR=0.47.
    'KsNb': 0.25,                       # mg-N/L; nitrogen half-saturation for benthic algae
    'KsPb': 0.125,                      # mg-P/L; phosphorus half-saturation for benthic algae
    'Ksb': 10.0,                        # g-D/m^2; benthic biomass half-saturation (FSb space limitation)
    'mub_max_20': 0.4,                  # 1/d; max benthic algae growth rate at 20 C
    'krb_20': 0.2,                      # 1/d; benthic algae respiration rate at 20 C
    'kdb_20': 0.3,                      # 1/d; benthic algae death rate at 20 C
    'mub_max_theta': 1.047,             # unitless; Arrhenius coefficient for growth
    'krb_theta': 1.06,                  # unitless; Arrhenius coefficient for respiration (note: 1.06 differs from typical 1.047)
    'kdb_theta': 1.047,                 # unitless; Arrhenius coefficient for death
    'b_growth_rate_option': 1,          # selector: 1=Multiplicative, 2=Min, 3=Harmonic
    'b_light_limitation_option': 1,     # selector: 1=Half-Saturation, 2=Smith, 3=Steele
    'Fw': 0.9,                          # unitless; fraction of benthic photosynthesis released as CO2 (1-Fw routed to DOC)
    'Fb': 0.9,                          # unitless; fraction of benthic biomass released in death (1-Fb routed to DOC)
}
