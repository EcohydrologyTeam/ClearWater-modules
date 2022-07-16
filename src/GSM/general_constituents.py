'''
=======================================================================================
General Constituent Simulation Module (GSM): General Constituent Kinetics Algorithm
=======================================================================================

Developed by:
* Dr. Todd E. Steissberg (ERDC-EL)
* Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

This module computes the water quality of a single computational cell. The algorithms 
and structure of this program were adapted from the Fortran 95 version of this module, 
developed by:
* Dr. Billy E. Johnson (ERDC-EL)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

Version 1.0

Initial Version: April 10, 2021
Last Revision Date: April 11, 2021
'''

# %%
from EnvironmentalSystems.ClearWater.water_quality_modules.python.clearwater_modules.src import water_quality_functions


def GeneralConstituentKinetics(GC: float, TwaterC: float, order: float, k_rc20: float = 1.0, k_theta: float = 1.047,
    rgc_rc20: float = 1.0, rgc_theta: float = 1.047, release: bool = True, settling: bool = False, 
    settling_rate: float = 0.1, depth: float = 1.0):

    '''
    Compute a single general constituent

    Parameters
    ----------
    GC : float
        General constituent concentration
    TwaterC : float
        Water temperature in degrees Celsius
    order : int
        Order of reaction kinetics (0, 1, or 2)
    k_rc20 : float, optional
        Reaction (decay) rate at 20 degrees Celsius
    k_theta : float, optional
        Arrhenius temperature correction factor for decay rate
    rgc_rc20 : float, optional
        Sediment release rate at 20 degrees Celsius
    rgc_theta : float, optional
        Arrhenius temperature correction factor for sediment release
    release : bool, optional
        Compute resuspension (on/off)
    settling : bool, optional
        Compute setting (i.e., bed loss, on/off)
    settling_rate : float, optional
        Settling rate (m/s)
    depth : float, optional
        Depth of the bed (meters)

    Returns
    ----------
    float
        Rate of change of general constituent concentration
    '''

    # Temperature corrections
    # k0_rc20 = k1_rc20 = k2_rc20 = 1.0
    # theta_rc20 = theta_rc20 = theta_rc20 = 1.047

    # Compute concentration changes
    gc_zero_order_decay: float = 0.0    # Rate of zero-order decay
    gc_first_order_decay: float = 0.0   # Rate of first-order decay
    gc_second_order_decay: float = 0.0  # Rate of second-order decay
    gc_from_bed: float = 0.0            # Sediment release rate (gain)
    gc_settling: float = 0.0            # Settling rate (loss)

    # Correct reaction rate for current temperature
    k_corr: float  = ArrheniusCorrection(TwaterC, k_rc20, k_theta)

    if order == 0:
        # Zero-order decay rate (mg/L/d)
        gc_zero_order_decay = k_corr
    elif order == 1:
        # First-order decay rate (1/d)
        gc_first_order_decay = k_corr * GC
    elif order == 2:
        gc_second_order_decay = k_corr * GC**2
    if release:
        rgc_corr: float  = ArrheniusCorrection(TwaterC, rgc_rc20, rgc_theta)
        gc_from_bed = rgc_corr / depth
    if settling:
        gc_settling = settling_rate * GC / depth

    # Compute net rate of change of general constituent concentration
    dGCdt = - gc_zero_order_decay - gc_first_order_decay - gc_second_order_decay + gc_from_bed - gc_settling

    return dGCdt


# %%
def test():
    GC = 10.0               # Initial concentration
    TwaterC = 25.0          # Water temperature
    order = 1               # Compute 1st order kinetics
    k_rc20 = 0.5            # Reaction rate at 20 degrees Celsius, decay
    k_theta = 1.047         # Arrhenius temperature correction factor, decay
    rgc_rc20 = 0.5          # Reaction rate at 20 degrees Celsius, settling
    rgc_theta = 1.047       # Arrhenius temperature correction factor, settling
    release = True          # Turn suspension on
    settling = True         # Turn settling on
    depth = 1.0             # Bed depth
    settling_rate = 0.0002  # Settling rate

    # Compute change of concentration
    dGCdt = GeneralConstituentKinetics(GC, TwaterC, order, k_rc20 = k_rc20, k_theta = k_theta, 
        rgc_rc20 = rgc_rc20, rgc_theta = rgc_theta, release = release, settling = settling, 
        depth = depth, settling_rate = settling_rate)

    GC_new = GC + dGCdt

    print("Initial concentration: %.2f" % GC)
    print("Change rate:           %.2f" % dGCdt)
    print("Final concentration:   %.2f" % GC_new)

# %%
# test()

# if __name__ == '__main__':
#     test()