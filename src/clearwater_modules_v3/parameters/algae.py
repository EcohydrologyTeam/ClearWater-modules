"""v3 NSM1 floating algae parameter defaults.

Consumed by: ``FloatingAlgae`` Process.

Defaults are drawn from published literature and from established
multi-site-validated water-quality models. Site-specific calibration
(e.g., for a particular river system, lake, or estuary) belongs in
the per-application YAML config, not here.

Source: v1 ``clearwater_modules/nsm1/constants.py`` ``AlgaeStaticVariables`` /
``DEFAULT_ALGAE``, with ``mu_max_20``, ``kdp_20``, ``krp_20`` updated to
the published consensus midpoint (see comment block below).
"""

DEFAULTS: dict[str, float | int | bool] = {
    'AWd': 100.0,                   # mg-D/ug-Chla; algal dry weight ratio (stoichiometry)
    'AWc': 40.0,                    # mg-C/ug-Chla; algal carbon ratio
    'AWn': 7.2,                     # mg-N/ug-Chla; algal nitrogen ratio
    'AWp': 1.0,                     # mg-P/ug-Chla; algal phosphorus ratio
    'AWa': 1000.0,                  # ug-Chla/ug-Chla; chlorophyll per algal unit (rda = AWd/AWa)
    'KL': 10.0,                     # W/m^2; light limitation half-saturation
    'KsN': 0.04,                    # mg-N/L; nitrogen half-saturation for algal growth
    'KsP': 0.0012,                  # mg-P/L; phosphorus half-saturation for algal growth
    #
    # mu_max_20, kdp_20, krp_20 — moved from v1 NSM1 defaults toward the
    # published literature consensus for mesotrophic-river algal kinetics.
    # The original v1 values (mu=1.0/d, kdp=0.15/d, krp=0.2/d) sit below
    # the typical published range and predict net algal *death* at typical
    # mesotrophic-river conditions (T~17 C, depth~2 m, mid-range light and
    # nutrient limitation give a daily-averaged net rate of ~-0.18/d),
    # which is inconsistent with chlorophyll observations in real
    # mesotrophic systems. The values below are central within published
    # ranges from three independently validated sources:
    #
    #   * Bowie et al. (1985), "Rates, Constants, and Kinetics
    #     Formulations in Surface Water Quality Modeling," 2nd ed.,
    #     EPA/600/3-85/040, Tables 6-1 and 6-13 (compilation of multi-site
    #     literature ranges).
    #   * Chapra, Pelletier, & Tao (2008), QUAL2K User Manual v2.11
    #     (default values used in the U.S. EPA-supported QUAL2K model,
    #     extensively validated across U.S. river systems).
    #   * Cole & Wells, CE-QUAL-W2 v4.5 User Manual (default values used
    #     in the U.S. ACE-supported W2 reservoir/river model).
    #
    # Site-specific calibration values (e.g., from individual case studies
    # like the Willamette mainstem) belong in the per-application YAML
    # ``parameters:`` override block, not in these defaults.
    #
    'mu_max_20': 2.0,               # 1/d; max algal growth rate at 20 C (QUAL2E typical 1.5-3.0; midpoint ~2.0)
    'kdp_20': 0.05,                 # 1/d; algal death rate at 20 C (QUAL2K default; Bowie typical 0.01-0.1)
    'krp_20': 0.10,                 # 1/d; algal respiration rate at 20 C (Bowie typical 0.05-0.10)
    'mu_max_theta': 1.047,          # unitless; Arrhenius coefficient for growth
    'kdp_theta': 1.047,             # unitless; Arrhenius coefficient for death
    'krp_theta': 1.047,             # unitless; Arrhenius coefficient for respiration
    'vsap': 0.15,                   # m/d; algal settling velocity (v1/Bowie consensus; site-specific values via YAML override)
    'growth_rate_option': 1,        # selector: 1=Multiplicative, 2=Min, 3=Harmonic
    'light_limitation_option': 1,   # selector: 1=Half-Saturation, 2=Smith, 3=Steele
}
