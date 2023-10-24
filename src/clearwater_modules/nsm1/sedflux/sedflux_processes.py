"""
File contains process to calculate sediment flux and associated dependent variables
"""

import math
from clearwater_modules.shared.processes import arrhenius_correction
import numba 
import numpy as np
import xarray as xr


# Note: OrgN, NH4, NO3, OrgP, TIP, POC and DOX must be on for SedFlux.

@numba.njit
def DOX_max(
    DOX: xr.DataArray,
) -> xr.DataArray:
    """Calculate DOX_max just to make sure DOX is 0.01 or greater (mg-O2/L).

    Args:
        DOX: Dissolve oxygen (mg-O2//L)
    """
    return xr.where(DOX>0.01, DOX, 0.01)


"""
    POC2=[0]*3        #concentration of sediment particulate organic carbon
    PON2=[0]*3        #concentration of sediment particulate organic nitrogen
    POP2=[0]*3        #concentration of sediment particulate organic phosphrous 

    POC2[1] = self.global_vars_sedflux['POC2_1']
    POC2[2] = self.global_vars_sedflux['POC2_2']
    POC2[3] = self.global_vars_sedflux['POC2_3']

    PON2[1] = self.global_vars_sedflux['PON2_1']
    PON2[2] = self.global_vars_sedflux['PON2_2']
    PON2[3] = self.global_vars_sedflux['PON2_3']

    POP2[1] = self.global_vars_sedflux['POP2_1']
    POP2[2] = self.global_vars_sedflux['POP2_2']
    POP2[3] = self.global_vars_sedflux['POP2_3']

    KPON_tc=[0]*3
    KPOP_tc=[0]*3
    KPOC_tc = [0]*3
  
"""

@numba.njit
def Dd_tc(
    Dd: xr.DataArray,
    TsedC: xr.DataArray,
) -> xr.DataArray:
    """Calculate Dd_tc: pore-water diffusion coefficient between layer 1 and 2 temperature corrected (m2/d).

    Args:
        Dd: pore-water diffusion coefficient between layer 1 and 2 (m2/d)
        TsedC: sediment temperature (C)
    """
    return arrhenius_correction(TsedC, Dd, 1.08)

@numba.njit
def Dp_tc(
    Dp: xr.DataArray,
    TsedC: xr.DataArray,
) -> xr.DataArray:
    """Calculate Dp_tc: particle mixing diffusion coefficient between layer 1 and 2 temperature corrected (m2/d).

    Args:
        Dd: particle mixing diffusion coefficient between layer 1 and 2 (m2/d)
        TsedC: sediment temperature (C)
    """
    return arrhenius_correction(TsedC, Dp, 1.117)

@numba.njit
def vnh41_tc(
    vnh41: xr.DataArray,
    TsedC: xr.DataArray,
) -> xr.DataArray:
    """Calculate vnh41_tc: nitrification reaction velocity in sediment layer 1 temperature corrected (m/d).

    Args:
        vnh41: nitrification reaction velocity in sediment layer 1 (m/d)
        TsedC: sediment temperature (C)
    """
    return vnh41 * 1.123 * ((TsedC-20)/2)

@numba.njit
def vno31_tc(
    vno31: xr.DataArray,
    TsedC: xr.DataArray,
) -> xr.DataArray:
    """Calculate vnh31_tc: denitrification reaction velocity in sediment layer 1 temperature corrected (m/d).

    Args:
        vnh41: denitrification reaction velocity in sediment layer 1 (m/d)
        TsedC: sediment temperature (C)
    """
    return vno31 * 1.08 * ((TsedC-20)/2)

@numba.njit
def vch41_tc(
    vch41: xr.DataArray,
    TsedC: xr.DataArray,
) -> xr.DataArray:
    """Calculate vch41_tc: methane oxidation reaction velocity in sediment layer 1 temperature corrected (m/d).

    Args:
        vch41: methane oxidation reaction velocity in sediment layer 1 (m/d)
        TsedC: sediment temperature (C)
    """
    return vch41 * 1.079 * ((TsedC-20)/2)

@numba.njit
def vh2sd_tc(
    vh2sd: xr.DataArray,
    TsedC: xr.DataArray,
) -> xr.DataArray:
    """Calculate vh2sd_tc: dissolve sulfide oxidation reaction velocity in sediment layer 1 temperature corrected (m/d).

    Args:
        vh2sd: dissolve sulfide oxidation reaction velocity in sediment layer 1 (m/d)
        TsedC: sediment temperature (C)
    """
    return vh2sd * 1.079 * ((TsedC-20)/2)

@numba.njit
def vh2sp_tc(
    vh2sp: xr.DataArray,
    TsedC: xr.DataArray,
) -> xr.DataArray:
    """Calculate vh2sp_tc: particulate sulfide oxidation reaction veolcoity in sediment layer 1 (m/d).

    Args:
        vh2sp: particulate sulfide oxidation reaction veolcoity in sediment layer 1 (m/d)
        TsedC: sediment temperature (C)
    """
    return vh2sp * 1.079 * ((TsedC-20)/2)

@numba.njit
def vno32_tc(
    vno32: xr.DataArray,
    TsedC: xr.DataArray,
) -> xr.DataArray:
    """Calculate vno32_tc: enitrification reaction velocity in sediment layer 3 temperature correction (m/d).

    Args:
        vh2sp: enitrification reaction velocity in sediment layer 2 (m/d)
        TsedC: sediment temperature (C)
    """

    return arrhenius_correction(TsedC, vno32, 1.08) #TODO why is this different from vno31 

@numba.njit
def KPON_1_tc(
  KPONG1: xr.DataArray,
  TsedC: xr.DataArray,

) -> xr.DataArray:
    """Calculate KPON_1_tc: diagenesis rate of PON G1 in sediment layer 2 temperature corrected (1/d)

    Args:
      KPONG1: diagenesis rate of PON G1 in sediment layer 2 (1/d)
      TsedC: sediment temperature (C)
    """

    return arrhenius_correction(TsedC, KPONG1, 1.1)

@numba.njit
def KPON_2_tc(
  KPONG2: xr.DataArray,
  TsedC: xr.DataArray,

) -> xr.DataArray:
    """Calculate KPON2__tc: diagenesis rate of PON G2 in sediment layer 2 temperature corrected (1/d)

    Args:
      KPONG2: diagenesis rate of PON G2 in sediment layer 2 (1/d)
      TsedC: sediment temperature (C)
    """

    return arrhenius_correction(TsedC, KPONG2, 1.1)

@numba.njit
def KPON_3_tc(

) -> xr.DataArray:
    """Calculate KPON_3_tc: diagenesis rate of PON G1 and G2 in sediment layer 2 temperature corrected (1/d)

    Args:

    """

    return 0 #TODO should this be a static variable

@numba.njit
def KPOP_1_tc(
  KPOPG1: xr.DataArray,
  TsedC: xr.DataArray,

) -> xr.DataArray:
    """Calculate KPOP_1_tc: diagenesis rate of POP G1 in sediment layer 2 temperature corrected (1/d)

    Args:
      KPOPG1: diagenesis rate of POP G1 in sediment layer 2 (1/d)
      TsedC: sediment temperature (C)
    """

    return arrhenius_correction(TsedC, KPOPG1, 1.1)

@numba.njit
def KPOP_2_tc(
  KPOPG2: xr.DataArray,
  TsedC: xr.DataArray,

) -> xr.DataArray:
    """Calculate KPOP_tc: diagenesis rate of POP G2 in sediment layer 2 temperature corrected (1/d)

    Args:
      KPOPG2: diagenesis rate of POP G2 in sediment layer 2 (1/d)
      TsedC: sediment temperature (C)
    """
    return arrhenius_correction(TsedC, KPOPG2, 1.15)

@numba.njit
def KPOP_3_tc(


) -> xr.DataArray:
    """Calculate KPOP_3_tc: diagenesis rate of POP in sediment layer 2 temperature corrected (1/d)

    Args:

    """
    return 0 #TODO should this be a static variable

@numba.njit
def KPOC_1_tc(
  KPOCG1: xr.DataArray,
  TsedC: xr.DataArray,

) -> xr.DataArray:
    """Calculate KPOC_1_c: diagenesis rate of POC G1 in sediment layer 2 temperature corrected (1/d)

    Args:
      KPOCG1: diagenesis rate of POC G1 in sediment layer 2 (1/d)
      TsedC: sediment temperature (C)
    """

    return arrhenius_correction(TsedC, KPOCG1, 1.1)

@numba.njit
def KPOC_2_tc(
  KPOCG2: xr.DataArray,
  TsedC: xr.DataArray,

) -> xr.DataArray:
    """Calculate KPOC_2_tc: diagenesis rate of POC G2 in sediment layer 2 temperature corrected (1/d)

    Args:
      KPOCG2: diagenesis rate of POC G2 in sediment layer 2 (1/d)
      TsedC: sediment temperature (C)
    """

    return arrhenius_correction(TsedC, KPOCG2, 1.15)

@numba.njit
def KPOC_3_tc(

) -> xr.DataArray:
    """Calculate KPOC_tc: diagenesis rate of POC G3 in sediment layer 2 temperature corrected (1/d)

    Args:

    """

    return 0 #TODO should this be a static variable

@numba.njit
def JPOC_1(
vsoc: xr.DataArray,
POC: xr.DataArray,
FPOC1: xr.DataArray,
use_POC: xr.DataArray,

ksbod_tc: xr.DataArray,
FCBOD1: xr.DataArray,
CBOD: xr.DataArray,
roc: xr.DataArray,

use_Algae:xr.DataArray,
FAP1:xr.DataArray,
rca:xr.DataArray,
Ap:xr.DataArray,
vsap:xr.DataArray,

use_Balgae: xr.DataArray,
AbDeath: xr.DataArray,
Fw: xr.DataArray,
rcb: xr.DataArray,
FAB1: xr.DataArray,

use_CBOD: xr.DataArray

) -> xr.DataArray:
    """Calculate JPOC_1: total depositional flux to sediment of particulate organic matter G1 (g-C/m^2/d)

    Args:
      vsoc: settling velocity of POC (m/d)
      POC: Particulate organic carbon concentration (mg-C/L)
      FPOC1: fraction of settled RPOC to sediment PCO G1 (unitless)
      FPOC2: fraction of settled RPOC to sediment PCO G2 (unitless)
      use_POC: true/false to use POC module (t/f)

      ksbod_tc: sedimentation rate (1/d)
      FCBOD1: fraction of CBOD sedimentation to G1 (unitless)
      FCBOD2: fraction of CBOD sedimentation to G2 (unitless)
      CBOD: carbonaceous biochemical oxygen demand (mg-O2/L)
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)

      use_Algae: true/false to use algae module (t/f)
      FAP1: fraction of settled algae to G1 (unitless)
      FAP2: fraction of settled algae to G2 (unitless)
      rca: Algal C " Chla ratio (mg-C/ugChla)
      Ap: Algal concentration (ugChla/L)
      vsap: Algal settling velocity (m/d)

      use_Balgae: true/false to use Balgae module (t/f)
      AbDeath: Balgal death rate (g/m^2/d)
      Fw: Fraction of benthic algae mortality into water column (unitless)
      rcb: Balgal C : Dry ratio (mg-C/mg-D)
      FAB1: fraction of benthic algae death to G1 (unitless)
      FAB2: fraction of benthic algae death to G1 (unitless)

      use_CBOD: true/false to use CBOD module (t/f)

    """

    #TODO add use_CBOD make sure that is appropriate.

    return xr.where(use_POC, vsoc * POC * FPOC1 , 0) + \
      xr.where(use_CBOD, ksbod_tc * FCBOD1 * CBOD / roc, 0) + \
      xr.where(use_Algae, FAP1 * rca * Ap * vsap, 0) + \
      xr.where(use_Balgae, AbDeath * (1.0 - Fw) * rcb * FAB1, 0)

@numba.njit
def JPOC_2(
vsoc: xr.DataArray,
POC: xr.DataArray,
FPOC2: xr.DataArray,
use_POC: xr.DataArray,

ksbod_tc: xr.DataArray,
FCBOD2: xr.DataArray,
CBOD: xr.DataArray,
roc: xr.DataArray,

use_Algae:xr.DataArray,
FAP2: xr.DataArray,
rca:xr.DataArray,
Ap:xr.DataArray,
vsap:xr.DataArray,

use_Balgae: xr.DataArray,
AbDeath: xr.DataArray,
Fw: xr.DataArray,
rcb: xr.DataArray,
FAB2: xr.DataArray,

use_CBOD: xr.DataArray

) -> xr.DataArray:
    """Calculate JPOC_2: total depositional flux to sediment of particulate organic matter G2 (g-C/m^2/d)

    Args:
      vsoc: settling velocity of POC (m/d)
      POC: Particulate organic carbon concentration (mg-C/L)
      FPOC1: fraction of settled RPOC to sediment PCO G1 (unitless)
      FPOC2: fraction of settled RPOC to sediment PCO G2 (unitless)
      use_POC: true/false to use POC module (t/f)

      ksbod_tc: sedimentation rate (1/d)
      FCBOD1: fraction of CBOD sedimentation to G1 (unitless)
      FCBOD2: fraction of CBOD sedimentation to G2 (unitless)
      CBOD: carbonaceous biochemical oxygen demand (mg-O2/L)
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)

      use_Algae: true/false to use algae module (t/f)
      FAP1: fraction of settled algae to G1 (unitless)
      FAP2: fraction of settled algae to G2 (unitless)
      rca: Algal C " Chla ratio (mg-C/ugChla)
      Ap: Algal concentration (ugChla/L)
      vsap: Algal settling velocity (m/d)

      use_Balgae: true/false to use Balgae module (t/f)
      AbDeath: Balgal death rate (g/m^2/d)
      Fw: Fraction of benthic algae mortality into water column (unitless)
      rcb: Balgal C : Dry ratio (mg-C/mg-D)
      FAB1: fraction of benthic algae death to G1 (unitless)
      FAB2: fraction of benthic algae death to G1 (unitless)

      use_CBOD: true/false to use CBOD module (t/f)

    """
    #TODO add use_CBOD make sure that is appropriate.
    return xr.where(use_POC, vsoc * POC * FPOC2 , 0) + \
      xr.where(use_CBOD, ksbod_tc * FCBOD2 * CBOD / roc, 0) + \
      xr.where(use_Algae, FAP2 * rca * Ap * vsap, 0) + \
      xr.where(use_Balgae, AbDeath * (1.0 - Fw) * rcb * FAB2, 0)

@numba.njit
def JPOC_3(
vsoc: xr.DataArray,
POC: xr.DataArray,
FPOC1: xr.DataArray,
FPOC2: xr.DataArray,
use_POC: xr.DataArray,

ksbod_tc: xr.DataArray,
FCBOD1: xr.DataArray,
FCBOD2: xr.DataArray,
CBOD: xr.DataArray,
roc: xr.DataArray,

use_Algae:xr.DataArray,
FAP1:xr.DataArray,
FAP2: xr.DataArray,
rca:xr.DataArray,
Ap:xr.DataArray,
vsap:xr.DataArray,

use_Balgae: xr.DataArray,
AbDeath: xr.DataArray,
Fw: xr.DataArray,
rcb: xr.DataArray,
FAB1: xr.DataArray,
FAB2: xr.DataArray,

use_CBOD: xr.DataArray

) -> xr.DataArray:
    """Calculate JPOC_3: total depositional flux to sediment of particulate organic matter g3 (g-C/m^2/d)

    Args:
      vsoc: settling velocity of POC (m/d)
      POC: Particulate organic carbon concentration (mg-C/L)
      FPOC1: fraction of settled RPOC to sediment PCO G1 (unitless)
      FPOC2: fraction of settled RPOC to sediment PCO G2 (unitless)
      use_POC: true/false to use POC module (t/f)

      ksbod_tc: sedimentation rate (1/d)
      FCBOD1: fraction of CBOD sedimentation to G1 (unitless)
      FCBOD2: fraction of CBOD sedimentation to G2 (unitless)
      CBOD: carbonaceous biochemical oxygen demand (mg-O2/L)
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)

      use_Algae: true/false to use algae module (t/f)
      FAP1: fraction of settled algae to G1 (unitless)
      FAP2: fraction of settled algae to G2 (unitless)
      rca: Algal C " Chla ratio (mg-C/ugChla)
      Ap: Algal concentration (ugChla/L)
      vsap: Algal settling velocity (m/d)

      use_Balgae: true/false to use Balgae module (t/f)
      AbDeath: Balgal death rate (g/m^2/d)
      Fw: Fraction of benthic algae mortality into water column (unitless)
      rcb: Balgal C : Dry ratio (mg-C/mg-D)
      FAB1: fraction of benthic algae death to G1 (unitless)
      FAB2: fraction of benthic algae death to G1 (unitless)

      use_CBOD: true/false to use CBOD module (t/f)

    """
    #TODO add use_CBOD make sure that is appropriate.
    return xr.where(use_POC, vsoc * POC * (1 - FPOC1 - FPOC2), 0) + \
      xr.where(use_CBOD, ksbod_tc * (1.0 - FCBOD1 - FCBOD2) * CBOD / roc, 0) + \
      xr.where(use_Algae, (1.0 - FAP1 - FAP2) * rca * Ap * vsap, 0) + \
      xr.where(use_Balgae, AbDeath * (1.0 - Fw) * rcb * (1.0 - FAB1 - FAB2), 0)

@numba.njit
def JPON_1(
vson: xr.DataArray,
OrgN: xr.DataArray,
FPON1: xr.DataArray,
use_OrgN: xr.DataArray,

use_Algae: xr.DataArray,
FAP1: xr.DataArray,
rna: xr.DataArray, 
Ap: xr.DataArray,
vsap: xr.DataArray,

use_Balgae: xr.DataArray,
AbDeath: xr.DataArray,
Fw: xr.DataArray, 
rnb: xr.DataArray,
FAB1: xr.DataArray,


) -> xr.DataArray:
    """Calculate JPON_2: total depositional flux of PON G1 from overlying water column to bed sediment (g-N/m^2/d)

    Args:
      vson: settling velocity of Organic Nitrogen (m/d)
      OrgN: Organic nitrogen concentration (mg-N/L)
      FPON1: fraction of settled RPON to sediment PON G1 (unitless)
      FPON2: fraction of settled RPON to sediment PON G2 (unitless)
      use_OrgN: true/false to use OrgN module (t/f)

      use_Algae: true/false to use algae module (t/f),
      FAP1: fraction of settled algae to G1 (unitless)
      FAP2: fraction of settled algae to G2 (unitless)
      rna: Algal N: Chla ratio (mg-N/ugChla)
      Ap: Algae concentration (ug-Chla/L)
      vsap: algal settling velocity (m/d)

      use_Balgae: true/false to use Balgae module (t/f),
      AbDeath: Balgal death rate (g/m2/d)
      Fw: fraction of benthic algae mortality into water column (unitless)
      rnb: Balgal N: Dry ratio (mg-N/mg-D)
      FAB1: fraction of benthic algae death to G1 (unitless)
      FAB2: fraction of benthic algae death to G2 (unitless)

    """

    return xr.where(use_OrgN, vson * OrgN * FPON1, 0) + \
      xr.where(use_Algae, FAP1 * rna * Ap * vsap, 0) + \
      xr.where(use_Balgae, AbDeath * (1.0 - Fw) * rnb * FAB1)

@numba.njit
def JPON_2(
vson: xr.DataArray,
OrgN: xr.DataArray,
FPON2: xr.DataArray,
use_OrgN: xr.DataArray,

use_Algae: xr.DataArray,
FAP2: xr.DataArray,
rna: xr.DataArray, 
Ap: xr.DataArray,
vsap: xr.DataArray,

use_Balgae: xr.DataArray,
AbDeath: xr.DataArray,
Fw: xr.DataArray, 
rnb: xr.DataArray,
FAB2: xr.DataArray,

) -> xr.DataArray:
    """Calculate JPON_2: total depositional flux of PON G2 from overlying water column to bed sediment (g-N/m^2/d)

    Args:
      vson: settling velocity of Organic Nitrogen (m/d)
      OrgN: Organic nitrogen concentration (mg-N/L)
      FPON1: fraction of settled RPON to sediment PON G1 (unitless)
      FPON2: fraction of settled RPON to sediment PON G2 (unitless)
      use_OrgN: true/false to use OrgN module (t/f)

      use_Algae: true/false to use algae module (t/f),
      FAP1: fraction of settled algae to G1 (unitless)
      FAP2: fraction of settled algae to G2 (unitless)
      rna: Algal N: Chla ratio (mg-N/ugChla)
      Ap: Algae concentration (ug-Chla/L)
      vsap: algal settling velocity (m/d)

      use_Balgae: true/false to use Balgae module (t/f),
      AbDeath: Balgal death rate (g/m2/d)
      Fw: fraction of benthic algae mortality into water column (unitless)
      rnb: Balgal N: Dry ratio (mg-N/mg-D)
      FAB1: fraction of benthic algae death to G1 (unitless)
      FAB2: fraction of benthic algae death to G2 (unitless)

    """

    return xr.where(use_OrgN, vson * OrgN * FPON2, 0) + \
      xr.where(use_Algae, FAP2 * rna * Ap * vsap, 0) + \
      xr.where(use_Balgae, AbDeath * (1.0 - Fw) * rnb * FAB2)

@numba.njit
def JPON_3(
vson: xr.DataArray,
OrgN: xr.DataArray,
FPON1: xr.DataArray,
FPON2: xr.DataArray,
use_OrgN: xr.DataArray,

use_Algae: xr.DataArray,
FAP1: xr.DataArray,
FAP2: xr.DataArray,
rna: xr.DataArray, 
Ap: xr.DataArray,
vsap: xr.DataArray,

use_Balgae: xr.DataArray,
AbDeath: xr.DataArray,
Fw: xr.DataArray, 
rnb: xr.DataArray,
FAB1: xr.DataArray,
FAB2: xr.DataArray,

) -> xr.DataArray:
    """Calculate JPON_3: total depositional flux of PON G3 from overlying water column to bed sediment (g-N/m^2/d)

    Args:
      vson: settling velocity of Organic Nitrogen (m/d)
      OrgN: Organic nitrogen concentration (mg-N/L)
      FPON1: fraction of settled RPON to sediment PON G1 (unitless)
      FPON2: fraction of settled RPON to sediment PON G2 (unitless)
      use_OrgN: true/false to use OrgN module (t/f)

      use_Algae: true/false to use algae module (t/f),
      FAP1: fraction of settled algae to G1 (unitless)
      FAP2: fraction of settled algae to G2 (unitless)
      rna: Algal N: Chla ratio (mg-N/ugChla)
      Ap: Algae concentration (ug-Chla/L)
      vsap: algal settling velocity (m/d)

      use_Balgae: true/false to use Balgae module (t/f),
      AbDeath: Balgal death rate (g/m2/d)
      Fw: fraction of benthic algae mortality into water column (unitless)
      rnb: Balgal N: Dry ratio (mg-N/mg-D)
      FAB1: fraction of benthic algae death to G1 (unitless)
      FAB2: fraction of benthic algae death to G2 (unitless)

    """
    return xr.where(use_OrgN, vson * OrgN * (1 - FPON1 - FPON2), 0) + \
      xr.where(use_Algae, (1.0 - FAP1 - FAP2) * rna * Ap * vsap, 0) + \
      xr.where(use_Balgae, AbDeath * (1.0 - Fw) * rnb * (1.0 - FAB1 - FAB2), 0)
        
@numba.njit
def JPOP_1(
vsop: xr.DataArray,
OrgP: xr.DataArray,
FPOP1: xr.DataArray,
use_OrgP: xr.DataArray,

use_Algae: xr.DataArray,
FAP1: xr.DataArray,
rpa: xr.DataArray, 
Ap: xr.DataArray, 
vsap: xr.DataArray,

use_Balgae: xr.DataArray,
AbDeath: xr.DataArray, 
Fw: xr.DataArray, 
rpb: xr.DataArray, 
FAB1: xr.DataArray,

) -> xr.DataArray:
    """Calculate JPOP_1: total depositional flux of POP G1 from overlying water column to bed sediment (g-P/m^2/d)

    Args:
      vsop: settling velocity of Organic Phosphorus (m/d)
      OrgP: Organic phosphorus concentration (mg-P/L)
      FPOP1: fraction of settled RPOP to sediment POP G1 (unitless)
      FPOP2: fraction of settled RPOP to sediment POP G2 (unitless)
      use_OrgP: true/false to use OrgP module (t/f)

      use_Algae: true/false to use algae module (t/f),
      FAP1: fraction of settled algae to G1 (unitless)
      FAP2: fraction of settled algae to G2 (unitless)
      rPa: Algal p: Chla ratio (mg-p/ugChla)
      Ap: Algae concentration (ug-Chla/L)
      vsap: algal settling velocity (m/d)

      use_Balgae: true/false to use Balgae module (t/f),
      AbDeath: Balgal death rate (g/m2/d)
      Fw: fraction of benthic algae mortality into water column (unitless)
      rpb: Balgal P: Dry ratio (mg-P/mg-D)
      FAB1: fraction of benthic algae death to G1 (unitless)
      FAB2: fraction of benthic algae death to G2 (unitless)

    """

    return xr.where(use_OrgP, vsop * OrgP * FPOP1, 0) + \
      xr.where(use_Algae, FAP1 * rpa * Ap * vsap, 0) + \
      xr.where(use_Balgae, AbDeath * (1.0 - Fw) * rpb * FAB1, 0)

@numba.njit
def JPOP_2(
vsop: xr.DataArray,
OrgP: xr.DataArray,
FPOP2: xr.DataArray,
use_OrgP: xr.DataArray,

use_Algae: xr.DataArray,
FAP2: xr.DataArray,
rpa: xr.DataArray, 
Ap: xr.DataArray, 
vsap: xr.DataArray,

use_Balgae: xr.DataArray,
AbDeath: xr.DataArray, 
Fw: xr.DataArray, 
rpb: xr.DataArray, 
FAB2: xr.DataArray,

) -> xr.DataArray:
    """Calculate JPOP_2: total depositional flux of POP G1 from overlying water column to bed sediment (g-P/m^2/d)

    Args:
      vsop: settling velocity of Organic Phosphorus (m/d)
      OrgP: Organic phosphorus concentration (mg-P/L)
      FPOP1: fraction of settled RPOP to sediment POP G1 (unitless)
      FPOP2: fraction of settled RPOP to sediment POP G2 (unitless)
      use_OrgP: true/false to use OrgP module (t/f)

      use_Algae: true/false to use algae module (t/f),
      FAP1: fraction of settled algae to G1 (unitless)
      FAP2: fraction of settled algae to G2 (unitless)
      rPa: Algal p: Chla ratio (mg-p/ugChla)
      Ap: Algae concentration (ug-Chla/L)
      vsap: algal settling velocity (m/d)

      use_Balgae: true/false to use Balgae module (t/f),
      AbDeath: Balgal death rate (g/m2/d)
      Fw: fraction of benthic algae mortality into water column (unitless)
      rpb: Balgal P: Dry ratio (mg-P/mg-D)
      FAB1: fraction of benthic algae death to G1 (unitless)
      FAB2: fraction of benthic algae death to G2 (unitless)

    """
    return xr.where(use_OrgP, vsop * OrgP * FPOP2, 0)+ \
      xr.where(use_Algae, FAP2 * rpa * Ap * vsap, 0) + \
      xr.where(use_Balgae, AbDeath * (1.0 - Fw) * rpb * FAB2, 0)

@numba.njit
def JPOP_3(
vsop: xr.DataArray,
OrgP: xr.DataArray,
FPOP1: xr.DataArray,
FPOP2: xr.DataArray,
use_OrgP: xr.DataArray,

use_Algae: xr.DataArray,
FAP1: xr.DataArray,
FAP2: xr.DataArray,
rpa: xr.DataArray, 
Ap: xr.DataArray, 
vsap: xr.DataArray,

use_Balgae: xr.DataArray,
AbDeath: xr.DataArray, 
Fw: xr.DataArray, 
rpb: xr.DataArray, 
FAB1: xr.DataArray,
FAB2: xr.DataArray,

) -> xr.DataArray:
    """Calculate JPOP_3: total depositional flux of POP G3 from overlying water column to bed sediment (g-P/m^2/d)

    Args:
      vsop: settling velocity of Organic Phosphorus (m/d)
      OrgP: Organic phosphorus concentration (mg-P/L)
      FPOP1: fraction of settled RPOP to sediment POP G1 (unitless)
      FPOP2: fraction of settled RPOP to sediment POP G2 (unitless)
      use_OrgP: true/false to use OrgP module (t/f)

      use_Algae: true/false to use algae module (t/f),
      FAP1: fraction of settled algae to G1 (unitless)
      FAP2: fraction of settled algae to G2 (unitless)
      rPa: Algal p: Chla ratio (mg-p/ugChla)
      Ap: Algae concentration (ug-Chla/L)
      vsap: algal settling velocity (m/d)

      use_Balgae: true/false to use Balgae module (t/f),
      AbDeath: Balgal death rate (g/m2/d)
      Fw: fraction of benthic algae mortality into water column (unitless)
      rpb: Balgal P: Dry ratio (mg-P/mg-D)
      FAB1: fraction of benthic algae death to G1 (unitless)
      FAB2: fraction of benthic algae death to G2 (unitless)

    """

    return xr.where(use_OrgP, vsop * OrgP * (1 - FPOP1 - FPOP2), 0)+ \
      xr.where(use_Algae, (1.0 - FAP1 - FAP2) * rpa * Ap * vsap, 0) + \
      xr.where(use_Balgae, AbDeath * (1.0 - Fw) * rpb * (1.0 - FAB1 - FAB2), 0)

"""
  #compute POC2/PON2/POP2, POC2/PON2/POP2 pathways and JC/JN/JP
    POC2=[0]*3        #concentration of sediment particulate organic carbon
    PON2=[0]*3        #concentration of sediment particulate organic nitrogen
    POP2=[0]*3        #concentration of sediment particulate organic phosphrous 
    POC2_Diagenesis=[0]*3
    PON2_Diagenesis=[0]*3
    POP2_Diagenesis=[0]*3
    POC2_Burial=[0]*3
    PON2_Burial=[0]*3
    POP2_Burial=[0]*3
"""

@numba.njit
def POC2_1_new(
SedFlux_solution_option: xr.DataArray,
JPOC_1: xr.DataArray,
KPOC_1_tc: xr.DataArray, 
h2: xr.DataArray,
vb: xr.DataArray, 
POC2_1: xr.DataArray, 
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate POC2_1_new: POC G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JPOC_1: otal depositional flux to sediment of particulate organic matter G1 (g-C/m^2/d)
      KPOC_1_tc: diagenesis rate of POC G1 in sediment layer 2 temperature corrected (1/d) 
      h2: active sediment layer (m)
      vb: burial velocity of POM2 in bed sediment (m/d)
      POC2_1: sediment POC layer 1 (mg-C/L)
      dt: time step (d)

    """

    return xr.where(SedFlux_solution_option==1, max(JPOC_1 / (KPOC_1_tc * h2 + vb),0), xr.where(SedFlux_solution_option==2, max((JPOC_1 + POC2_1 * h2 / dt),0) / (h2 / dt + KPOC_1_tc * h2 + vb),0))

@numba.njit
def POC2_2_new(
SedFlux_solution_option: xr.DataArray,
JPOC_2: xr.DataArray,
KPOC_2_tc: xr.DataArray, 
h2: xr.DataArray,
vb: xr.DataArray, 
POC2_2: xr.DataArray, 
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate POC2_2_new: POC G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JPOC_2: total depositional flux to sediment of particulate organic matter G2 (g-C/m^2/d)
      KPOC_2_tc: diagenesis rate of POC G2 in sediment layer 2 temperature corrected (1/d) 
      h2: active sediment layer (m)
      vb: burial velocity of POM2 in bed sediment (m/d)
      POC2_2: sediment POC layer 2 (mg-C/L)
      dt: time step (d)

    """

    return xr.where(SedFlux_solution_option==1, max(JPOC_2 / (KPOC_2_tc * h2 + vb),0), xr.where(SedFlux_solution_option==2, max((JPOC_2 + POC2_2 * h2 / dt),0) / (h2 / dt + KPOC_2_tc * h2 + vb),0))

@numba.njit
def POC2_3_new(
SedFlux_solution_option: xr.DataArray,
JPOC_3: xr.DataArray,
KPOC_3_tc: xr.DataArray, 
h2: xr.DataArray,
vb: xr.DataArray, 
POC2_3: xr.DataArray, 
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate POC2_3_new: POC G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JPOC_3: total depositional flux to sediment of particulate organic matter G3 (g-C/m^2/d)
      KPOC_3_tc: diagenesis rate of POC G3 in sediment layer 2 temperature corrected (1/d) 
      h2: active sediment layer (m)
      vb: burial velocity of POM2 in bed sediment (m/d)
      POC2_3: sediment POC layer 2 (mg-C/L)
      dt: time step (d)

    """

    return xr.where(SedFlux_solution_option==1, max(JPOC_3 / (KPOC_3_tc * h2 + vb),0), xr.where(SedFlux_solution_option==2, max((JPOC_3 + POC2_3 * h2 / dt),0) / (h2 / dt + KPOC_3_tc * h2 + vb),0))

@numba.njit
def PON2_1_new(
SedFlux_solution_option: xr.DataArray,
JPON_1: xr.DataArray,
KPON_1_tc: xr.DataArray, 
h2: xr.DataArray,
vb: xr.DataArray, 
PON2_1: xr.DataArray, 
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate PON2_1_new: PON G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JPON_1: total depositional flux to sediment of particulate organic nitrogen G1 (g-N/m^2/d)
      KPON_1_tc: diagenesis rate of PON G1 in sediment layer 2 temperature corrected (1/d) 
      h2: active sediment layer (m)
      vb: burial velocity of POM2 in bed sediment (m/d)
      PON2_1: sediment PON layer 2 (mg-N/L)
      dt: time step (d)

    """

    return xr.where(SedFlux_solution_option==1, max(JPON_1 / (KPON_1_tc * h2 + vb),0), xr.where(SedFlux_solution_option==2, max((JPON_1 + PON2_1 * h2 / dt),0) / (h2 / dt + KPON_1_tc * h2 + vb),0))

@numba.njit
def PON2_2_new(
SedFlux_solution_option: xr.DataArray,
JPON_2: xr.DataArray,
KPON_2_tc: xr.DataArray, 
h2: xr.DataArray,
vb: xr.DataArray, 
PON2_2: xr.DataArray, 
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate PON2_2_new: PON G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JPON_2: total depositional flux to sediment of particulate organic nitrogen G2 (g-N/m^2/d)
      KPON_2_tc: diagenesis rate of PON G2 in sediment layer 2 temperature corrected (1/d) 
      h2: active sediment layer (m)
      vb: burial velocity of POM2 in bed sediment (m/d)
      PON2_2: sediment PON layer 2 (mg-N/L)
      dt: time step (d)

    """

    return xr.where(SedFlux_solution_option==1, max(JPON_2 / (KPON_2_tc * h2 + vb),0), xr.where(SedFlux_solution_option==2, max((JPON_2 + PON2_2 * h2 / dt),0) / (h2 / dt + KPON_2_tc * h2 + vb),0))

@numba.njit
def PON2_3_new(
SedFlux_solution_option: xr.DataArray,
JPON_3: xr.DataArray,
KPON_3_tc: xr.DataArray, 
h2: xr.DataArray,
vb: xr.DataArray, 
PON2_3: xr.DataArray, 
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate PON2_3_new: PON G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JPON_3: total depositional flux to sediment of particulate organic nitrogen G3 (g-N/m^2/d)
      KPON_3_tc: diagenesis rate of PON G3 in sediment layer 2 temperature corrected (1/d) 
      h2: active sediment layer (m)
      vb: burial velocity of POM2 in bed sediment (m/d)
      PON2_3: sediment PON layer 2 (mg-N/L)
      dt: time step (d)

    """

    return xr.where(SedFlux_solution_option==1, max(JPON_3 / (KPON_3_tc * h2 + vb),0), xr.where(SedFlux_solution_option==2, max((JPON_3 + PON2_3 * h2 / dt),0) / (h2 / dt + KPON_3_tc * h2 + vb),0))

@numba.njit
def POP2_1_new(
SedFlux_solution_option: xr.DataArray,
JPOP_1: xr.DataArray,
KPOP_1_tc: xr.DataArray, 
h2: xr.DataArray,
vb: xr.DataArray, 
POP2_1: xr.DataArray, 
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate POP2_1_new: POP G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JPOP_1: total depositional flux to sediment of particulate organic phosphorus G1 (g-P/m^2/d)
      KPOP_1_tc: diagenesis rate of POP G1 in sediment layer 2 temperature corrected (1/d) 
      h2: active sediment layer (m)
      vb: burial velocity of POM2 in bed sediment (m/d)
      PON2_1: sediment POP layer 2 (mg-P/L)
      dt: time step (d)

    """

    return xr.where(SedFlux_solution_option==1, max(JPOP_1 / (KPOP_1_tc * h2 + vb),0), xr.where(SedFlux_solution_option==2, max((JPOP_1 + POP2_1 * h2 / dt),0) / (h2 / dt + KPOP_1_tc * h2 + vb),0))

@numba.njit
def POP2_2_new(
SedFlux_solution_option: xr.DataArray,
JPOP_2: xr.DataArray,
KPOP_2_tc: xr.DataArray, 
h2: xr.DataArray,
vb: xr.DataArray, 
POP2_2: xr.DataArray, 
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate POP2_2_new: POP G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JPOP_2: total depositional flux to sediment of particulate organic phosphorus G2 (g-P/m^2/d)
      KPOP_2_tc: diagenesis rate of POP G2 in sediment layer 2 temperature corrected (1/d) 
      h2: active sediment layer (m)
      vb: burial velocity of POM2 in bed sediment (m/d)
      PON2_2: sediment POP layer 2 (mg-P/L)
      dt: time step (d)

    """

    return xr.where(SedFlux_solution_option==1, max(JPOP_2 / (KPOP_2_tc * h2 + vb),0), xr.where(SedFlux_solution_option==2, max((JPOP_2 + POP2_2 * h2 / dt),0) / (h2 / dt + KPOP_2_tc * h2 + vb),0))

@numba.njit
def POP2_3_new(
SedFlux_solution_option: xr.DataArray,
JPOP_3: xr.DataArray,
KPOP_3_tc: xr.DataArray, 
h2: xr.DataArray,
vb: xr.DataArray, 
POP2_3: xr.DataArray, 
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate POP2_3_new: POP G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JPOP_3: total depositional flux to sediment of particulate organic phosphorus G3 (g-P/m^2/d)
      KPOP_3_tc: diagenesis rate of POP G3 in sediment layer 2 temperature corrected (1/d) 
      h2: active sediment layer (m)
      vb: burial velocity of POM2 in bed sediment (m/d)
      PON2_3: sediment POP layer 2 (mg-P/L)
      dt: time step (d)

    """

    return xr.where(SedFlux_solution_option==1, max(JPOP_3 / (KPOP_3_tc * h2 + vb),0), xr.where(SedFlux_solution_option==2, max((JPOP_3 + POP2_3 * h2 / dt),0) / (h2 / dt + KPOP_3_tc * h2 + vb),0))

@numba.njit
def POC2_Diagenesis_1(
KPOC_1_tc: xr.DataArray,
h2: xr.DataArray,
POC2_1_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POC2_Diagenesis_1: sediment diagenesis of POC G1 in sediment layer (g-C/m2/d)

    Args:
      KPOC_1_tc: diagenesis rate of POC G1 in sediment layer 2 temperature corrected (1/d)
      h2: active sediment layer (m)
      POC2_1_new: POC G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)

    """

    return KPOC_1_tc * h2 * POC2_1_new

@numba.njit
def POC2_Diagenesis_2(
KPOC_2_tc: xr.DataArray,
h2: xr.DataArray,
POC2_2_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POC2_Diagenesis_2: sediment diagenesis of POC G2 in sediment layer (g-C/m2/d)

    Args:
      KPOC_2_tc: diagenesis rate of POC G2 in sediment layer 2 temperature corrected (1/d)
      h2: active sediment layer (m)
      POC2_2_new: POC G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)

    """

    return KPOC_2_tc * h2 * POC2_2_new

@numba.njit
def POC2_Diagenesis_3(
KPOC_3_tc: xr.DataArray,
h2: xr.DataArray,
POC2_3_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POC2_Diagenesis_3: sediment diagenesis of POC G3 in sediment layer (g-C/m2/d)

    Args:
      KPOC_3_tc: diagenesis rate of POC G3 in sediment layer 2 temperature corrected (1/d)
      h2: active sediment layer (m)
      POC2_3_new: POC G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)

    """

    return KPOC_3_tc * h2 * POC2_3_new

@numba.njit
def PON2_Diagenesis_1(
KPON_1_tc: xr.DataArray,
h2: xr.DataArray,
PON2_1_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate PON2_Diagenesis_1: sediment diagenesis of PON G1 in sediment layer (g-N/m2/d)

    Args:
      KPON_1_tc: diagenesis rate of PON G1 in sediment layer 2 temperature corrected (1/d)
      h2: active sediment layer (m)
      PON2_1_new: PON G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)

    """

    return KPON_1_tc * h2 * PON2_1_new

@numba.njit
def PON2_Diagenesis_2(
KPON_2_tc: xr.DataArray,
h2: xr.DataArray,
PON2_2_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate PON2_Diagenesis_2: sediment diagenesis of PON G2 in sediment layer (g-N/m2/d)

    Args:
      KPON_2_tc: diagenesis rate of PON G2 in sediment layer 2 temperature corrected (1/d)
      h2: active sediment layer (m)
      PON2_2_new: PON G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)

    """

    return KPON_2_tc * h2 * PON2_2_new

@numba.njit
def PON2_Diagenesis_3(
KPON_3_tc: xr.DataArray,
h2: xr.DataArray,
PON2_3_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate PON2_Diagenesis_3: sediment diagenesis of PON G3 in sediment layer (g-N/m2/d)

    Args:
      KPON_3_tc: diagenesis rate of PON G3 in sediment layer 2 temperature corrected (1/d)
      h2: active sediment layer (m)
      PON2_3_new: PON G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)

    """

    return KPON_3_tc * h2 * PON2_3_new

@numba.njit
def POP2_Diagenesis_1(
KPOP_1_tc: xr.DataArray,
h2: xr.DataArray,
POP2_1_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POP2_Diagenesis_1: sediment diagenesis of PON G1 in sediment layer (g-P/m2/d)

    Args:
      KPOP_1_tc: diagenesis rate of POP G1 in sediment layer 2 temperature corrected (1/d)
      h2: active sediment layer (m)
      POP2_1_new: POP G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)

    """

    return KPOP_1_tc * h2 * POP2_1_new

@numba.njit
def POP2_Diagenesis_2(
KPOP_2_tc: xr.DataArray,
h2: xr.DataArray,
POP2_2_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POP2_Diagenesis_2: sediment diagenesis of PON G2 in sediment layer (g-P/m2/d)

    Args:
      KPOP_2_tc: diagenesis rate of POP G2 in sediment layer 2 temperature corrected (1/d)
      h2: active sediment layer (m)
      POP2_2_new: POP G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)

    """

    return KPOP_2_tc * h2 * POP2_2_new

@numba.njit
def POP2_Diagenesis_3(
KPOP_3_tc: xr.DataArray,
h2: xr.DataArray,
POP2_3_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POP2_Diagenesis_3: sediment diagenesis of PON G3 in sediment layer (g-P/m2/d)

    Args:
      KPOP_3_tc: diagenesis rate of POP G3 in sediment layer 2 temperature corrected (1/d)
      h2: active sediment layer (m)
      POP2_3_new: POP G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)

    """

    return KPOP_3_tc * h2 * POP2_3_new

@numba.njit
def POC2_Burial_1(
vb: xr.DataArray,
POC2_1_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POC2_Burial_1: sediment diagenesis of POC G1 in sediment layer (g-C/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      POC2_1_new: POC G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)

    """

    return vb * POC2_1_new

@numba.njit
def POC2_Burial_2(
vb: xr.DataArray,
POC2_2_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POC2_Burial_1: sediment diagenesis of POC G2 in sediment layer (g-C/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      POC2_2_new: POC G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)

    """

    return vb * POC2_2_new

@numba.njit
def POC2_Burial_3(
vb: xr.DataArray,
POC2_3_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POC2_Burial_3: sediment diagenesis of POC G3 in sediment layer (g-C/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      POC2_3_new: POC G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)

    """

    return vb * POC2_3_new

@numba.njit
def PON2_Burial_1(
vb: xr.DataArray,
PON2_1_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate PON2_Burial_1: sediment diagenesis of PON G1 in sediment layer (g-N/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      PON2_1_new: PON G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)

    """

    return vb * PON2_1_new

@numba.njit
def PON2_Burial_2(
vb: xr.DataArray,
PON2_2_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate PON2_Burial_2: sediment diagenesis of PON G2 in sediment layer (g-N/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      PON2_2_new: PON G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)

    """

    return vb * PON2_2_new

@numba.njit
def PON2_Burial_3(
vb: xr.DataArray,
PON2_3_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate PON2_Burial_3: sediment diagenesis of PON G3 in sediment layer (g-N/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      PON2_3_new: PON G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)

    """

    return vb * PON2_3_new

@numba.njit
def POP2_Burial_1(
vb: xr.DataArray,
POP2_1_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POP2_Burial_1: sediment diagenesis of POP G1 in sediment layer (g-P/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      POP2_1_new: POP G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)

    """

    return vb * POP2_1_new

@numba.njit
def POP2_Burial_2(
vb: xr.DataArray,
POP2_2_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POP2_Burial_2: sediment diagenesis of POP G2 in sediment layer (g-P/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      POP2_2_new: POP G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)

    """

    return vb * POP2_2_new

@numba.njit
def POP2_Burial_3(
vb: xr.DataArray,
POP2_3_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate POP2_Burial_3: sediment diagenesis of POP G3 in sediment layer (g-P/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      POP2_3_new: POP G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)

    """

    return vb * POP2_3_new

@numba.njit
def JC(
POC2_Diagenesis_1: xr.DataArray,
POC2_Diagenesis_2: xr.DataArray,

) -> xr.DataArray:
    """Calculate JC:  total sediment diagenesis flux of POC  (g-C/m2/d)

    Args:
      POC2_Diagenesis_1: sediment diagenesis of POC G1 in sediment layer (g-C/m2/d)
      POC2_Diagenesis_2: sediment diagenesis of POC G2 in sediment layer (g-C/m2/d)

    """

    return POC2_Diagenesis_1 + POC2_Diagenesis_2

@numba.njit
def JN(
PON2_Diagenesis_1: xr.DataArray,
PON2_Diagenesis_2: xr.DataArray,

) -> xr.DataArray:
    """Calculate JN:  total sediment diagenesis flux of PON  (g-N/m2/d)

    Args:
      PON2_Diagenesis_1: sediment diagenesis of PON G1 in sediment layer (g-N/m2/d)
      PON2_Diagenesis_2: sediment diagenesis of PON G2 in sediment layer (g-N/m2/d)

    """

    return PON2_Diagenesis_1 + PON2_Diagenesis_2

@numba.njit
def JP(
POP2_Diagenesis_1: xr.DataArray,
POP2_Diagenesis_2: xr.DataArray,

) -> xr.DataArray:
    """Calculate JP:  total sediment diagenesis flux of POP  (g-P/m2/d)

    Args:
      POP2_Diagenesis_1: sediment diagenesis of POP G1 in sediment layer (g-P/m2/d)
      POP2_Diagenesis_2: sediment diagenesis of POP G2 in sediment layer (g-P/m2/d)

    """

    return POP2_Diagenesis_1 + POP2_Diagenesis_2

@numba.njit
def SOD_Bed(
JC: xr.DataArray,
roc: xr.DataArray,
JN: xr.DataArray,

) -> xr.DataArray:
    """Calculate SOD_Bed: SedFlux sediment oxygen demand (g-O2/m2/d)

    Args:
      JC: total sediment diagenesis flux of POC (g-C/m2/d)
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)
      JN: total sediment diagenesis flux of PON (g-N/m2/d)
    """

    return JC * roc + 1.714 * JN

"""
    # compute SOD
    SOD_old =0
    # coefficients of a quadratic equation
    ra2, ra1, ra0 =0
    #root of the quadratic equation      
    sn1, disc, r1, r2 =0     
"""  

@numba.njit
def BFORmax(
TsedC: xr.DataArray,
TempBen: xr.DataArray,
KsDp: xr.DataArray,
DOX: xr.DataArray,

) -> xr.DataArray:
    """Calculate BFORmax: maximum benthic stress oxygen correction coefficient (unitless)

    Args:
      TsedC: sediment temperature (C)
      TempBen: critical temperature for benthic stress (C)
      KsDp: half-saturation constant for oxygen in particle mixing (mg-O2/L)
      DOX: dissolved oxygen (mg-O/L)
    """
    #TODO should this be 'new'
    return xr.where (TsedC<TempBen, 0.0, max(0.0, KsDp / KsDp + DOX))

@numba.njit  
def ST_new(
SedFlux_solution_option: xr.DataArray,
TsedC: xr.DataArray,
TempBen: xr.DataArray,
ST: xr.DataArray,
dt: xr.DataArray,
KsDp: xr.DataArray,
DOX: xr.DataArray,
kst: xr.DataArray,
BFORmax: xr.DataArray,

) -> xr.DataArray:
    """Calculate ST_new: Sediment Benthic Stress: Low DO will eliminate bioturbation. Particle phase mixing coefficient is modified (d)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      TsedC: sediment temperature (C)
      TempBen: critical temperature for benthic stress (C)
      ST: Sediment Benthic Stress (d)
      dt: timestep (d)
      KsDp: half-saturation constant for oxygen in particle mixing (mg-O2/L)
      DOX: dissolved oxygen (mg-O/L)
      kst: decay rate of benthic stress (1/d)
      BFORmax: maximum benthic stress oxygen correction coefficient (unitless)
    """
    #TODO maybe a better way to do this. i do not even know if it is possible to get an nan value?
    to_return = xr.where(SedFlux_solution_option == 1, 0.0, xr.where (TsedC<TempBen, (ST + dt * KsDp / (KsDp + DOX)) / (1.0 + kst * dt), (ST + dt * BFORmax) / (1.0 + kst * dt),0))
    to_return = xr.where (math.isnan(to_return), 0, to_return)
    return to_return

@numba.njit  
def w12(
SedFlux_solution_option: xr.DataArray,
Dp_tc: xr.DataArray,
h2: xr.DataArray,
POC2_1_new: xr.DataArray,
POCr: xr.DataArray,
Css2: xr.DataArray,
TsedC: xr.DataArray,
TempBen: xr.DataArray,
kst: xr.DataArray,
ST_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      Dp_tc: particle mixing diffusion coefficient between layer 1 and 2 temperature corrected (m2/d)
      h2: active sediment layer (m)
      POC2_1_new: POC G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)
      POCr: reference POC G1 concentration for bioturbation (mg-C/L)
      Css2: solids concentration in sediment layer 2 (kg/L)
      TsedC: sediment temperature (C)
      TempBen: critical temperature for benthic stress (C)
      kst: decay rate of benthic stress (1/d)
      ST_new: Sediment Benthic Stress: Low DO will eliminate bioturbation. Particle phase mixing coefficient is modified (d)

    """

    to_return = xr.where(SedFlux_solution_option==1, Dp_tc/(0.5 *h2), xr.where(SedFlux_solution_option==2, Dp_tc / (0.5 * h2) * (POC2_1_new / (1000.0 * POCr * Css2),0)))
    to_return = xr.where(TsedC>TempBen, to_return * (1-kst * ST_new),to_return)
    to_return = xr.where(to_return <0, 0, to_return)

    #TODO the 1000 does not cancel with anything if that is what it is supposed to do
    #TODO probably a better way to write this
    return to_return

@numba.njit  
def KL12(
Dd_tc: xr.DataArray,
h2: xr.DataArray,

) -> xr.DataArray:
    """Calculate KL12: Dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)

    Args:
      Dd_tc: pore-water diffusion coefficient between layer 1 and 2 temperature corrected (m2/d)
      h2: active sediment layer (m)

    """

    return xr.where(math.isnan(Dd_tc / (0.5 * h2)), 0 , Dd_tc / (0.5 * h2))

@numba.njit  
def KL01(
SOD_Bed: xr.DataArray,
DOX: xr.DataArray,

) -> xr.DataArray:
    """Calculate KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)

    Args:
      SOD_Bed: SedFlux sediment oxygen demand (g-O2/m2/d)
      DOX: dissolved oxygen (mg-O/L) 

    """

    return xr.where(math.isnan(SOD_Bed / DOX) or (SOD_Bed / DOX) ==0, 1.0E-8, SOD_Bed / DOX)

@numba.njit  
def fds1(
Css1: xr.DataArray,
kdh2s2: xr.DataArray,

) -> xr.DataArray:
    """Calculate fds1: dissolved fraction for H2S1 and H2S2 in layer 1 (unitless)

    Args:
      Css1: solids concentration in sediment layer 1 (kg/L)
      kdh2s2: partition coefficient for sulfide in sediment layer 2 (L/kg)

    """

    return 1.0 / (1.0 + Css1 * kdh2s2)

@numba.njit  
def fds2(
Css2: xr.DataArray,
kdh2s2: xr.DataArray,

) -> xr.DataArray:
    """Calculate fds2: dissolved fraction for H2S1 and H2S2 in layer 2 (unitless)

    Args:
      Css2: solids concentration in sediment layer 1 (kg/L)
      kdh2s2: partition coefficient for sulfide in sediment layer 2 (L/kg)

    """

    return 1.0 / (1.0 + Css2 * kdh2s2)

@numba.njit  
def fdp1(
fds1: xr.DataArray,

) -> xr.DataArray:
    """Calculate fdp1: particulate fraction for H2S1 and H2S2 in layer 1 (unitless)

    Args:
      fds1: dissolved fraction for H2S1 and H2S2 in layer 1 (unitless)

    """

    return 1.0 - fds1

@numba.njit  
def fdp2(
fds2: xr.DataArray,

) -> xr.DataArray:
    """Calculate fdp2: particulate fraction for H2S1 and H2S2 in layer 2 (unitless)

    Args:
      fds2: dissolved fraction for H2S1 and H2S2 in layer 2 (unitless)

    """

    return 1.0 - fds2

@numba.njit  
def fdp2(
fds2: xr.DataArray,

) -> xr.DataArray:
    """Calculate fdp2: particulate fraction for H2S1 and H2S2 in layer 2 (unitless)

    Args:
      fds2: dissolved fraction for H2S1 and H2S2 in layer 2 (unitless)

    """

    return 1.0 - fds2

@numba.njit  
def fd1(
Css1: xr.DataArray,
kdnh42: xr.DataArray,

) -> xr.DataArray:
    """Calculate fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)

    Args:
      Css1: solids concentration in sediment layer 1 (kg/L)
      kdnh42: partition coefficient for ammonium in sediment layer 2 (L/kg)

    """

    return 1.0 / (1.0 + Css1 * kdnh42)

@numba.njit  
def fd2(
Css2: xr.DataArray,
kdnh42: xr.DataArray,

) -> xr.DataArray:
    """Calculate fd2: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 2 (unitless)

    Args:
      Css2: solids concentration in sediment layer 2 (kg/L)
      kdnh42: partition coefficient for ammonium in sediment layer 2 (L/kg)

    """

    return 1.0 / (1.0 + Css2 * kdnh42)

@numba.njit  
def fp1(
fd1: xr.DataArray,

) -> xr.DataArray:
    """Calculate fp1: fraction of inorganic matter (ammonia, phosphate) in particulate form in sediment layer 1 (unitless)

    Args:
      fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)

    """

    return 1.0 - fd1

@numba.njit  
def fp2(
fd2: xr.DataArray,

) -> xr.DataArray:
    """Calculate fp2: fraction of inorganic matter (ammonia, phosphate) in particulate form in sediment layer 2 (unitless)

    Args:
      fd2: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 2 (unitless)

    """
    return 1.0 - fd2

@numba.njit  
def TNH41(
NH41: xr.DataArray,
fd1: xr.DataArray,

) -> xr.DataArray:
    """Calculate THN41: total concentration NH4 dissolved layer 1 (mg-N/L)

    Args:
      NH41: NH4 seidment layer 1 (mg-N/L)
      fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)

    """
    return NH41 / fd1

@numba.njit  
def TNH42(
NH42: xr.DataArray,
fd2: xr.DataArray,

) -> xr.DataArray:
    """Calculate THN42: total concentration NH4 dissolved layer 2 (mg-N/L)

    Args:
      NH42: NH4 seidment layer 2 (mg-N/L)
      fd2: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 2 (unitless)

    """
    return NH42 / fd2

@numba.njit  
def FOxna(
DOX: xr.DataArray,
KsOxna1: xr.DataArray,

) -> xr.DataArray:
    """Calculate FOxna: nitrification attenuation due to low oxygen in layer 1 (unitless)

    Args:
      DOX: dissolved oxgyen concentration (mg-O2/L)
      KsOxna1: half-saturation oxygen constant for sediment nitrification (mg-O2/L)

    """
    return xr.where(math.isnan(DOX / (KsOxna1 * 2.0 + DOX)),0, DOX / (KsOxna1 * 2.0 + DOX))

@numba.njit  
def con_nit(
vnh41_tc: xr.DataArray,
FOxna: xr.DataArray,
fd1: xr.DataArray,

) -> xr.DataArray:
    """Calculate con_nit: something nitrogen (m2/d2)

    Args:
      vnh41_tc: nitrification reaction velocity in sediment layer 1 (m/d)
      FOxna: nitrification attenuation due to low oxygen in layer 1 (unitless)
      fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)

    """
    return vnh41_tc * vnh41_tc * FOxna * fd1 
 
@numba.njit  
def a12_TNH4(
w12: xr.DataArray,
fp2: xr.DataArray,
KL12: xr.DataArray,
fd2: xr.DataArray,

) -> xr.DataArray:
    """Calculate a12_TNH4: coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2) (m/d)

    Args:
      w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)
      fp2: fraction of inorganic matter (ammonia, phosphate) in particulate form in sediment layer 2 (unitless)
      KL12: KL12: Dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      fd2: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 2 (unitless)

    """
    return -w12 * fp2 - KL12 * fd2

@numba.njit  
def a21_TNH4(
w12: xr.DataArray,
fp1: xr.DataArray,
KL12: xr.DataArray,
fd1: xr.DataArray,
vb: xr.DataArray,

) -> xr.DataArray:
    """Calculate a12_TNH4: coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2) (m/d)

    Args:
      w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)
      fp1: fraction of inorganic matter (ammonia, phosphate) in particulate form in sediment layer 1 (unitless)
      KL12: KL12: Dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
    """
    return -w12 * fp1 - KL12 * fd1 - vb

@numba.njit  
def a22_TNH4(
SedFlux_solution_option: xr.DataArray,
a12_TNH4: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,
vb: xr.DataArray,

) -> xr.DataArray:
    """Calculate a22_TNH4: coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2) (m/d)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      a12_TNH4: coefficents for implicit finite difference form for TNH4 (m/d)
      h2: active sediment layer (m)
      dt: timestep (d)
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
    """
    return xr.where(SedFlux_solution_option == 1, -a12_TNH4 + vb, -a12_TNH4 + vb + h2 / dt)

@numba.njit  
def b2_TNH4(
SedFlux_solution_option: xr.DataArray,
JN: xr.DataArray,
TNH42: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,
vb: xr.DataArray,

) -> xr.DataArray:
    """Calculate b2_TNH4: coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2) (g-N/m2/d)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JN: total sediment diagenesis flux of PON (g-N/m2/d)
      THN42: total concentration NH4 dissolved layer 2 (mg-N/L)
      h2: active sediment layer (m)
      dt: timestep (d)
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
    """
    return xr.where(SedFlux_solution_option == 1, JN, JN + h2 * TNH42 / dt)

@numba.njit  
def a12_NO3(
KL12: xr.DataArray,

) -> xr.DataArray:
    """Calculate a12_NO3: coefficents for implicit finite difference form for NO3 (a11, a12, a21, a22, b1, b2) (m/d)

    Args:
      KL12: issolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)

    """
    return -KL12

@numba.njit  
def a21_NO3(
KL12: xr.DataArray,

) -> xr.DataArray:
    """Calculate a21_NO3: coefficents for implicit finite difference form for NO3 (a11, a12, a21, a22, b1, b2) (m/d)

    Args:
      KL12: issolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)

    """
    return -KL12

@numba.njit  
def a22_NO3(
SedFlux_solution_option: xr.DataArray,
KL12: xr.DataArray,
vno32_tc: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate a22_NO3: coefficents for implicit finite difference form for NO3 (a11, a12, a21, a22, b1, b2) (m/d))

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      h2: active sediment layer (m)
      vno32_tc: nitrification reaction velocity in sediment layer 3 temperature correction (m/d). 
      dt: timestep (d)
      KL12: issolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)

    """
    return xr.where(SedFlux_solution_option == 1, KL12 + vno32_tc, KL12 + vno32_tc + h2 / dt)

@numba.njit  
def b2_NO3(
SedFlux_solution_option: xr.DataArray,
NO32: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate b2_NO3: coefficents for implicit finite difference form for No3 (a11, a12, a21, a22, b1, b2) (g-N/m2/d)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      h2: active sediment layer (m)
      dt: timestep (d)
      NO32: nitrate concentration sediment layer 2 (mg-N/L)

    """
    return xr.where(SedFlux_solution_option == 1, 0.0, h2 * NO32 / dt)

@numba.njit  
def SO4(
Salinity: xr.DataArray,
roso4: xr.DataArray,
SO4_fresh: xr.DataArray,

) -> xr.DataArray:
    """Calculate SO4: estimate based on salinity for salt water (mg-O2/L)
    
    Args:
      Salinity: water salinity (ppt)
      roso4: oxygen stoichiometric coeff for sulfate SO4 (g-O2/g-SO4)
      SO4_fresh: SO4 concentration of overlaying water column in freshwater. User-defined if not simulated for fresh water (mg-O2/L)

    """
    #TODO make sure this does not have to be SO4_new and that the units are correct
    #TODO find the formula for if true
    return xr.where(Salinity > 0.01, (20.0 + 27.0 / 190.0 * 607.445 * Salinity) * roso4, SO4_fresh)

@numba.njit  
def TH2S2_prev(
TH2S2: xr.DataArray,

) -> xr.DataArray:
    """Calculate TH2S2_prev: previous timestep TH2S2 concentration in sediment layer 2 (mg-O/L)
    
    Args:
      TH2S2: TH2S sediment layer 2 (mg-O/L)

    """
    #TODO should this be a constant/will they need to used with the _new variables.
    return TH2S2

@numba.njit  
def SO42_prev(
SO42: xr.DataArray,

) -> xr.DataArray:
    """Calculate SO42_prev: previous timestep SO42 concentration in sediment layer 2 (mg-O/L)
    
    Args:
      SO42: SO42 sediment layer 2 (mg-O/L)

    """
    #TODO should this be a constant/will they need to used with the _new variables.
    return SO42

@numba.njit  
def HSO4_new (
POCdiagenesis_part_option: xr.DataArray,
t: xr.DataArray,
Dd_tc: xr.DataArray,
SO4: xr.DataArray,
h2: xr.DataArray,
roc: xr.DataArray,
JC: xr.DataArray,
JCc: xr.DataArray,

) -> xr.DataArray:
    """Calculate HSO4_new: previous timestep SO42 concentration in sediment layer 2 (mg-O/L)
    
    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless)
      t: #TODO what is this time something?
      Dd_tc: pore-water diffusion coefficient between layer 1 and 2 temperature corrected (m2/d)
      SO4: SO4: estimate based on salinity for salt water (mg-O2/L)
      h2: active sediment layer (m)
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)
      JC: total sediment diagenesis flux of POC (g-C/m2/d)
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)

    """
    #TODO possibly better way to code

    to_return = xr.where(POCdiagenesis_part_option == 2 and t< 1.0E-10, math.sqrt(Dd_tc * SO4 * h2 / (max(roc * JC, 1.0E-10))), math.sqrt(2.0 * Dd_tc * SO4 * h2 / JCc))
    to_return = xr.where(to_return > h2, h2, to_return)
    return to_return

    # Half-saturation method
    if self.sedFlux_constants['POCdiagenesis_part_option'] == 1 :  

            con_sox = (vh2sd_tc * vh2sd_tc * fds1 + vh2sp_tc *
                       vh2sp_tc * fps1) * DOX / 2.0 / KSh2s
            if math. isnan(con_sox):
                con_sox = (vh2sd_tc * vh2sd_tc * fds1 +
                           vh2sp_tc * vh2sp_tc * fps1)

        # Sulfate reduction depth method
        elif POCdiagenesis_part_option == 2:
            # Set initial value for HSO4
            if (t < 1.0E-10):
                # sulfate penetration into layer 2 from layer 1
                HSO4 = math.sqrt(Dd_tc * SO4 * h2 / (max(roc * JC, 1.0E-10)))
                if (HSO4 > h2):
                    HSO4 = h2

            # SO41 and SO42
            TH2S1_prev = TH2S1
            HSO4_prev = HSO4
            SO42_prev = SO42
            a12_SO4 = - KL12
            a21_SO4 = KL12

            if SedFlux_solution_option == 1:
                a22_SO4 = - KL12
            elif SedFlux_solution_option == 2:
                a22_SO4 = - KL12 - h2 / dt

            # TH2S1 and TH2S2
            fds1 = 1.0 / (1.0 + Css1 * kdh2s2)
            fds2 = 1.0 / (1.0 + Css2 * kdh2s2)
            fps1 = 1.0 - fds1
            fps2 = 1.0 - fds2
            TH2S2_prev = TH2S2
            con_sox = (vh2sd_tc * vh2sd_tc * fds1 + vh2sp_tc *
                       vh2sp_tc * fps1) * DOX / 2.0 / KSh2s
            if math.isnan(con_sox):
                con_sox = (vh2sd_tc * vh2sd_tc * fds1 +
                           vh2sp_tc * vh2sp_tc * fps1)

            a12_TH2S = -w12 * fps2 - KL12 * fds2
            a21_TH2S = -w12 * fps1 - KL12 * fds1 - vb
            if SedFlux_solution_option == 1:
                a22_TH2S = -a12_TH2S + vb
            elif SedFlux_solution_option == 2:
                a22_TH2S = -a12_TH2S + vb + h2 / dt

    # CH41 and CH42
    CH4sat = 100.0 * (1.0 + depth / 10.0) * 1.024**(20.0 - TsedC)
    if Methane_solution_option ==2 :
      CH42_prev = CH42
      FOxch = DOX / (KsOxch * 2.0 + DOX)
      if math.isnan(FOxch) :
        FOxch = 0.0
      con_cox = vch41_tc * vch41_tc * FOxch
      a12_CH4 = -KL12
      a21_CH4 = -KL12
      if SedFlux_solution_option == 1 :
        a22_CH4 = KL12
      elif SedFlux_solution_option == 2 :
        a22_CH4 = KL12 + h2 / dt

    #compute SOD
    for i in (1, int(maxit)) :
      #TNH41 and TNH42
      if (KsNh4 > 0.0) :
        FNH4 = KsNh4 / (KsNh4 + fd1 * TNH41)
      else :
        FNH4  = 1.0

      a11 = -a21_TNH4 + con_nit * FNH4 / KL01 + KL01 * fd1
      b1  = KL01 * NH4
      TNH41, TNH42 = MatrixSolution(TNH41, TNH42, a11, a12_TNH4, b1, a21_TNH4, a22_TNH4, b2_TNH4)
      TNH41 = max(TNH41, 0.0)
      TNH42 = max(TNH42, 0.0)
      
      # NO31 and NO32
      a11 = -a21_NO3 + vno31_tc * vno31_tc / KL01 + KL01
      b1  = con_nit * FNH4 / KL01 * TNH41 + KL01 * NO3
      NO31, NO32 = MatrixSolution(NO31, NO32, a11, a12_NO3, b1, a21_NO3, a22_NO3, b2_NO3)
      NO31 = max(NO31, 0.0)
      NO32 = max(NO32, 0.0)
      
      JC_dn = rondn * (vno31_tc * vno31_tc / KL01 * NO31 + vno32_tc * NO32)
      JCc = max(roc * JC - JC_dn, 1.0E-10)
      
      # SO41, SO42 and TH2S1, TH2S2
      #coefficients for SO41, SO42, H2S1 and H2S2 equations

      bx=[0]*4
      ad=[(0,0,0,0),(0,0,0,0), (0,0,0,0), (0,0,0,0)] # 4x4
      
      #coefficients for SO41 and SO42 equations    
      g=[0]*2
      h=[(0,0),(0,0)] #2x2          
  
      # Half-saturation method
      if POCdiagenesis_part_option == 1 :
        #compute HSO4
        HSO4 = math.sqrt(2.0 * Dd_tc * SO4 * h2 / JCc) 
        if (HSO4 > h2):
          HSO4 = h2
        if HSO4 == 0.0:
          KL12SO4 = 1.0       # set a large value (1) for KL12SO4. 
        else:
          KL12SO4 = Dd_tc / (0.5 * HSO4)

        # four equations
        # 0 = bx(0) + ad(0,0) * SO41 + ad(0,1) * SO42 + ad(0,2) * TH2S1
        # 0 = bx(1) + ad(1,0) * SO41 + ad(1,1) * SO42 - JCc * SO42 / (SO42 + KsSO4)
        # 0 = bx(2) + ad(2,2) * TH2S1 + ad(2,3) * TH2S2
        # 0 = bx(3) + ad(3,2) * TH2S1 + ad(3,3) * TH2S2 + JCc * SO42 / (SO42 + KsSO4)

        bx[0] = KL01 * SO4                     
        bx[2] = KL01 * H2S                     
        ad[0][0] = - KL01 - KL12SO4
        ad[0][1] = KL12SO4
        ad[0][2] = con_sox / KL01
        ad[2][2] = - vb - fps1 * w12 - fds1 * KL01 - fds1 * KL12SO4 - con_sox / KL01
        ad[2][3] = w12 * fps2 + KL12SO4 * fds2
        
        ad[1][0] = KL12SO4
        ad[3][2] = vb + fps1 * w12 + fds1 * KL12SO4
        if SedFlux_solution_option == 1 :
          ad[1][1] = - KL12SO4
          bx[1]   = 0.0
          ad[3][3] = - vb - KL12SO4 * fds2 - w12 * fps2
          bx[3]   = 0.0
        else:
          ad[1][1] = - KL12SO4 - h2 / dt
          bx[1]   = h2 * SO42_prev / dt
          ad[3][3] = - vb - KL12SO4 * fds2 - w12 * fps2 - h2 / dt 
          bx[3]   = h2 * TH2S2_prev / dt

        #eliminate H2ST1 and H2ST2 from above equation sets, get two equations
        # 0 = g(0) + h(0,0) * SO41 + h(0,1) * SO42
        # 0 = g(1) + h(1,0) * SO41 + h(1,1) * SO42 + JCc * SO42 / (SO42 + KsSO4)
        
        g[1] = ((bx[0] * ad[2][2] - ad[0][2] * bx[2]) * ad[3][3] - bx[0] * ad[2][3] * ad[3][2]) / (ad[0][2] * ad[2][3]) + bx[3]   
        g[0] = g[1] + bx[1]  
        h[1][0] = ad[0][0] * (ad[2][2] * ad[3][3] - ad[2][3] * ad[3][2]) / (ad[0][2] * ad[2][3])
        h[1][1] = ad[0][1] * (ad[2][2] * ad[3][3] - ad[2][3] * ad[3][2]) / (ad[0][2] * ad[2][3])
        h[0][0] = h[1][0] + ad[1][0] 
        h[0][1] = h[1][1] + ad[1][1]
        
        #eliminate SO41 and get a quadratic equation of SO42
        # ra2 * SO42 * SO42 + ra1 * SO42 + ra0 = 0.0 
        ra0 = (h[0][0] * g[1] - h[1][0] * g[0]) * KsSO4
        ra1 = h[0][0] * g[1] - h[1][0] * g[0] + (h[0][0] * h[1][1] - h[0][1] * h[1][0]) * KsSO4 + h[0][0] * JCc
        ra2 = h[0][0] * h[1][1] - h[0][1] * h[1][0]
      
        #solve SO42
        sn1 = 1.0
        if (ra1 <= 0.0):
          sn1 = - 1.0
        disc = (- ra1 - sn1 * math.sqrt(ra1 * ra1 - 4.0 * ra2 * ra0)) / 2.0
        if (disc != 0.0) :
          r1 = disc / ra2
          r2 = ra0 / disc
        else:
          if (ra2 == 0.0) :
            r1 = - ra0 / ra1
            r2 = r1
          else:   #TODO what does this calculate? This is when ra0 = 0 
            r1 = -ra1 / ra2
  
        SO42 = r1
        if (SO42 < 0.0) :
          SO42 = r2
        
        #solve SO41, H2ST1, H2ST2
        SO41  = - (g[0] + h[0][1] * SO42) / h[0][0]
        TH2S1 = - (bx[0] + ad[0][0] * SO41 + ad[0][1] * SO42) / ad[0][2]
        TH2S2 = - (bx[2] + ad[2][2] * TH2S1) / ad[2][3]
        CSOD_H2S = con_sox / KL01 * TH2S1
        JCc_CH4  = JCc * KsSO4 / (SO42 + KsSO4)       #TODO find this formula. See page 170 is SO42 supposed to be on top as well?
        
      #sulfate reduction depth method
      elif POCdiagenesis_part_option == 2 :
        
        if (SO4 <= 0.1) :
          HSO4    = 0.0
          KL12SO4 = 1.0
          SO41    = SO4
          SO42    = SO4
        else:
          
          # first compute HSO4 TODO I do not know where the equation comes from: see page 159 seems to be close
          if SedFlux_solution_option == 1:
            ra2 = 2.0 * KL01 * JCc / h2
            ra1 = 2.0 * Dd_tc * JCc / h2
            ra0 = - 2.0 * Dd_tc * (KL01 * SO4 + con_sox / KL01 * TH2S1_prev)      
          elif SedFlux_solution_option == 2 :
            ra2 = 2.0 * KL01 * JCc / h2 + (KL01 * SO4 + con_sox / KL01 * TH2S1_prev) / dt
            ra1 = 2.0 * Dd_tc * JCc / h2 - 2.0 * KL01 * HSO4_prev * SO42_prev / dt
            ra0 = - 2.0 * Dd_tc * (KL01 * SO4 + con_sox / KL01 * TH2S1_prev + HSO4_prev * SO42_prev / dt)

          HSO4 = (- ra1 + math.sqrt(ra1 * ra1 - 4.0 * ra2 * ra0)) / 2.0 / ra2
          if (HSO4 > h2) :
            HSO4 = h2
          
          # TH2S1, TH2S2
          a11 = -a21_TH2S + con_sox / KL01 + KL01 * fds1
          b1  = KL01 * H2S
          if SedFlux_solution_option == 1 :  
            b2  = JCc * HSO4 / h2
          elif SedFlux_solution_option == 2:
            b2  = JCc * HSO4 / h2 + TH2S2_prev * h2 / dt

          TH2S1, TH2S2 = MatrixSolution (TH2S1, TH2S2, a11, a12_TH2S, b1, a21_TH2S, a22_TH2S, b2)
          TH2S1 = max(TH2S1, 0.0)
          TH2S2 = max(TH2S2, 0.0)

          CSOD_H2S = con_sox / KL01 * TH2S1
          JCc_CH4  = JCc * (h2 - HSO4) / h2     #TODO should it be JCc_CH4 = JCc * (1-(HSO4/H2)) page 170 in PDF
          
          # SO41, SO42
          if (HSO4 == h2) :
            KL12SO4 = KL12
            a11  = KL01 + KL12
            b1   = KL01 * SO4 + CSOD_H2S
            if SedFlux_solution_option == 1 : #TODO check thes b2
              b2 = JCc
            elif SedFlux_solution_option == 2 :
              b2 = JCc - SO42_prev * HSO4_prev / dt

            #TODO TODO  call MatrixSolution(SO41, SO42, a11, a12_SO4, b1, a21_SO4, a22_SO4, b2) 
              SO41 = max(SO41, 0.0)
              SO42 = max(SO42, 0.0)  

          else :
            KL12SO4 = Dd_tc / (0.5 * HSO4)
            SO41 = (KL01 * SO4 + CSOD_H2S) / (KL01 + KL12SO4 * 0.5)
            SO42 = SO41 / 2.0

        # TH2S1, TH2S2
        a11 = -a21_TH2S + con_sox / KL01 + KL01 * fds1
        b1  = KL01 * H2S
        if SedFlux_solution_option == 1 :
          b2  = JCc * HSO4 / h2
        elif SedFlux_solution_option == 2 :
          b2  = JCc * HSO4 / h2 + TH2S2_prev * h2 / dt

        TH2S1, TH2S2 = MatrixSolution(TH2S1, TH2S2, a11, a12_TH2S, b1, a21_TH2S, a22_TH2S, b2)
        TH2S1 = max(TH2S1, 0.0)
        TH2S2 = max(TH2S2, 0.0)
        CSOD_H2S = con_sox / KL01 * TH2S1
        JCc_CH4  = JCc * (h2 - HSO4) / h2

      # CH41 and CH42

      # analytical solutions
      if Methane_solution_option == 1 :
        CSODmax  = min(math.sqrt(2.0 * KL12 * CH4sat * JCc_CH4), JCc_CH4)
        
        if vch41_tc / KL01 < 100.0 :
          CSOD_CH4 = CSODmax * (1.0 - 2.0 / (math.exp(vch41_tc / KL01) + math.exp(-vch41_tc / KL01))) 
        else:
          CSOD_CH4 = CSODmax

      #numerical solutions
      elif Methane_solution_option == 2 : 
        a11 = KL12 + con_cox / KL01 + KL01
        b1  = KL01 * CH4
        if SedFlux_solution_option == 1 : 
          b2 = JCc_CH4
        elif SedFlux_solution_option == 2 :
          b2 = JCc_CH4 + CH42_prev * h2 / dt
    
        CH41, CH42 = MatrixSolution(CH41, CH42, a11, a12_CH4, b1, a21_CH4, a22_CH4, b2)
        
        if CH42 > CH4sat :
          CH42 = CH4sat
          CH41 = (b1 - a12_CH4 * CH42) / a11

        CH41     = max(CH41, 0.0)
        CH42     = max(CH42, 0.0)
        CSOD_CH4 = con_cox / KL01 * CH41

      #update SOD
      NSOD    = ron * con_nit * FNH4 / KL01 * TNH41     #TODO potentially missing a dn1 equation 5.37 page 162
      SOD_old = SOD_Bed
      SOD_Bed  = (CSOD_CH4 + CSOD_H2S + NSOD + SOD_Bed) / 2.0
      if (abs(SOD_Bed - SOD_old) / SOD_Bed * 100.0 < res) :
        break    
      KL01 = SOD_Bed / DOX
      if math.isnan(KL01) or KL01 == 0.0:
        KL01 = 1.0E-8

    #TODO not sure if i is defined outside of the loop
    #determine whether SOD is converged
    if i > int(maxit) :
      print('SOD iterations exceeded.')

    # inorganic species reactions and mass transfers
    kdpo41=0
    hsat=0    #depth where methane reaches saturation
    
    KL01 = SOD_Bed / DOX
      
    # pathways of TNH41/2, NO31/2, CH41/2, SO41/2, TH2S1/2, DIC1/2, TIP1/2
    #TNH41 and TNH42
    NH41= fd1 * TNH41
    NH42= fd2 * TNH42
    
    JNH4= KL01 * (NH41 - NH4)
    TNH41_Burial= vb * TNH41
    NH41_Nitrification = con_nit * FNH4 / KL01 * TNH41
    #PNH41_PNH42 = vb * (fp2 * TNH42 - fp1 * TNH41)
    NH41_NH42 = KL12 * (NH42 - NH41)
    TNH42_Burial = vb * TNH42
    
    # NO31 and NO32
    JNO3 = KL01 * (NO31 - NO3)
    NO31_Denit = vno31_tc * vno31_tc / KL01 * NO31
    NO31_NO32 = KL12 * (NO32 - NO31)
    NO32_Denit = vno32_tc * NO32
    
    # SO41 and SO42
    if POCdiagenesis_part_option == 1 :
      JCc_SO4 = JCc * SO42 / (SO42 + KsSO4)
    elif POCdiagenesis_part_option == 2 :
      JCc_SO4 = JCc * HSO4 / h2

    JSO4 = KL01 * (SO41 - SO4)
    SO41_SO42 = KL12SO4 * (SO42 - SO41)
    
    # TH2S1 and TH2S2
    H2S1 = fds1 * TH2S1
    H2S2 = fds2 * TH2S2
    
    JH2S = KL01 * (H2S1 - H2S)
    H2S1_Oxidation = con_sox / KL01 * TH2S1
    TH2S1_Burial = vb * TH2S1
    H2S1_H2S2 = KL12 * (H2S2 - H2S1)
    
    #PH2S1_PH2S2 = vb  * (TH2S2 * fps2 - TH2S1 * fps1)
    TH2S2_Burial = vb * TH2S2
    
    # CH41 and CH42
    if Methane_solution_option == 1 :
      CH41_Oxidation = CSOD_CH4
      JCH4 = CSODmax - CSOD_CH4
      JCH4g = JCc_CH4 - CSODmax
      if vch41_tc <= 0 :
        CH41 = 0.0
      else :
        CH41 = CSOD_CH4 / (vch41_tc * vch41_tc / KL01)

      #CH42 is not computed???
      CH42 = 0.0                                        
      
    elif Methane_solution_option == 2 :
      CH41_Oxidation = con_cox / KL01 * CH41
      JCH4 = KL01 * (CH41 - CH4)
      if CH42 == CH4sat :
        JCH4g = JCc_CH4 - JCH4 - CH41_Oxidation - (CH42 - CH42_prev) / dt * h2
      else:
        JCH4g = 0.0

    # DIC1 and DIC2
    a11 = KL01 + KL12
    a12 = -KL12
    b1 = KL01 * DIC * 12000.0 + CH41_Oxidation / 2.0 / roc + rcdn * NO31_Denit
    a21   = -KL12
    if SedFlux_solution_option == 1 :
      a22 = KL12
      b2  = (JCc_CH4 / 2.0 + JCc_SO4) / roc  + rcdn * NO32_Denit
    elif SedFlux_solution_option == 2 :
      a22 = KL12 + h2 / dt
      b2  = (JCc_CH4 / 2.0 + JCc_SO4) / roc + rcdn * NO32_Denit + DIC2 * h2 / dt

    DIC1, DIC2 = MatrixSolution(DIC1, DIC2, a11, a12, b1, a21, a22, b2)
    DIC1  = max(DIC1, 0.0)
    DIC2  = max(DIC2, 0.0)
    
    JDIC = KL01 * (DIC1 - DIC * 12000.0)
    DIC1_CH41_Oxidation = CH41_Oxidation / 2.0 / roc
    DIC1_NO31_Denit = rcdn * NO31_Denit
    DIC1_DIC2 = KL12 * (DIC2 - DIC1)
    DIC2_POC2_SO42 = JCc_SO4 / roc
    DIC2_CH42 = JCc_CH4 / 2.0 / roc
    DIC2_NO32_Denit = rcdn * NO32_Denit
    
    # TIP1 and TIP2
    if DOX >= DOcr : 
      kdpo41 = kdpo42 * d_kpo41
    else :
      kdpo41 = kdpo42 * d_kpo41**(DOX / DOcr)

    fd1 = 1.0 / (1.0 + Css1 * kdpo41)
    fd2 = 1.0 / (1.0 + Css2 * kdpo42)
    fp1 = 1.0 - fd1
    fp2 = 1.0 - fd2
    
    a21 = -w12 * fp1 - KL12 * fd1 - vb
    a11 = -a21 + KL01 * fd1
    a12 = -w12 * fp2 - KL12 * fd2
    b1  = KL01 * fdp * TIP

    if SedFlux_solution_option == 1 :
      a22 = -a12 + vb
      b2  = JP + TIP * (1.0 - fdp) * vs
    elif SedFlux_solution_option == 2 :
      a22 = -a12 + vb + h2 / dt
      b2  =  JP + TIP * (1.0 - fdp) * vs + h2 * TIP2 / dt

    TIP1, TIP2 = MatrixSolution(TIP1, TIP2, a11, a12, b1, a21, a22, b2)
    TIP1 = max(TIP1, 0.0)
    TIP2 = max(TIP2, 0.0)
    
    DIP1 = TIP1 * fd1
    DIP2 = TIP2 * fd2
    
    JDIP = KL01 * (DIP1 - fdp * TIP)
    TIP_TIP2 = TIP * (1.0 - fdp) * vs
    TIP1_Burial = vb * TIP1
    DIP1_DIP2 = KL12 * (DIP2 - DIP1)
    #PIP1_PIP2 = w12  * (TIP2 * fp2 - TIP1 * fp1)
    TIP2_Burial = vb * TIP2

  # solve mass balance equations by a matrix solution

  '''
# output sediment diagenesis pathways
  subroutine SedFluxPathwayOutput(na, a)
    integer  :: na
    real(R8) :: a(na)
    !
    do i = 1, 3
      if (JPOC_index(i) > 0)            a(JPOC_index(i))            = JPOC(i)
      if (JPON_index(i) > 0)            a(JPON_index(i))            = JPON(i)
      if (JPOP_index(i) > 0)            a(JPOP_index(i))            = JPOP(i)
      if (POC2_Diagenesis_index(i) > 0) a(POC2_Diagenesis_index(i)) = POC2_Diagenesis(i)
      if (PON2_Diagenesis_index(i) > 0) a(PON2_Diagenesis_index(i)) = PON2_Diagenesis(i)
      if (POP2_Diagenesis_index(i) > 0) a(POP2_Diagenesis_index(i)) = POP2_Diagenesis(i)
      if (POC2_Burial_index(i) > 0)     a(POC2_Burial_index(i))     = POC2_Burial(i)
      if (PON2_Burial_index(i) > 0)     a(PON2_Burial_index(i))     = PON2_Burial(i)
      if (POP2_Burial_index(i) > 0)     a(POP2_Burial_index(i))     = POP2_Burial(i)
    end do
    !
    if (JC_index > 0)    a(JC_index)    = JC
    if (JC_dn_index > 0) a(JC_dn_index) = JC_dn
    if (JN_index > 0)    a(JN_index)    = JN
    if (JP_index > 0)    a(JP_index)    = JP
    !
    if (w12_index > 0)      a(w12_index)      = w12
    if (KL12_index > 0)     a(KL12_index)     = KL12
    if (KL01_index > 0)     a(KL01_index)     = KL01
    if (SOD_Bed_index > 0)  a(SOD_Bed_index)  = SOD_Bed
    !
    if (JNH4_index > 0)                 a(JNH4_index)               = JNH4
    if (TNH41_Burial_index > 0)         a(TNH41_Burial_index)       = TNH41_Burial             
    if (NH41_Nitrification_index > 0)   a(NH41_Nitrification_index) = NH41_Nitrification
    !if (PNH41_PNH42_index > 0)          a(PNH41_PNH42_index)        = PNH41_PNH42
    if (NH41_NH42_index > 0)            a(NH41_NH42_index)          = NH41_NH42
    if (TNH42_Burial_index > 0)         a(TNH42_Burial_index)       = TNH42_Burial
    !
    if (JNO3_index > 0)                 a(JNO3_index)               = JNO3
    if (NO31_Denit_index > 0)           a(NO31_Denit_index)         = NO31_Denit
    if (NO31_NO32_index > 0)            a(NO31_NO32_index)          = NO31_NO32
    if (NO32_Denit_index > 0)           a(NO32_Denit_index)         = NO32_Denit
    !
    if (CH4sat_index > 0)               a(CH4sat_index)             = CH4sat
    if (JCH4_index > 0)                 a(JCH4_index)               = JCH4
    if (CH41_Oxidation_index > 0)       a(CH41_Oxidation_index)     = CH41_Oxidation
    if (JCc_CH4_index > 0)              a(JCc_CH4_index)            = JCc_CH4
    if (JCH4g_index > 0)                a(JCH4g_index)              = JCH4g
    !
    if (JSO4_index > 0)                 a(JSO4_index)               = JSO4
    if (JCc_SO4_index > 0)              a(JCc_SO4_index)            = JCc_SO4
    if (SO41_SO42_index > 0)            a(SO41_SO42_index)          = SO41_SO42
    !
    if (JH2S_index > 0)                 a(JH2S_index)               = JH2S
    if (H2S1_Oxidation_index > 0)       a(H2S1_Oxidation_index)     = H2S1_Oxidation
    if (TH2S1_Burial_index > 0)         a(TH2S1_Burial_index)       = TH2S1_Burial
    if (H2S1_H2S2_index > 0)            a(H2S1_H2S2_index)          = H2S1_H2S2
    !if (PH2S1_PH2S2_index > 0)          a(PH2S1_PH2S2_index)        = PH2S1_PH2S2
    if (TH2S2_Burial_index > 0)         a(TH2S2_Burial_index)       = TH2S2_Burial
    !
    if (JDIC_index > 0)                 a(JDIC_index)               = JDIC
    if (DIC1_CH41_Oxidation_index > 0)  a(DIC1_CH41_Oxidation_index)= DIC1_CH41_Oxidation
    if (DIC1_NO31_Denit_index > 0)      a(DIC1_NO31_Denit_index)    = DIC1_NO31_Denit
    if (DIC1_DIC2_index > 0)            a(DIC1_DIC2_index)          = DIC1_DIC2
    if (DIC2_POC2_SO42_index > 0)       a(DIC2_POC2_SO42_index)     = DIC2_POC2_SO42
    if (DIC2_CH42_index > 0)            a(DIC2_CH42_index)          = DIC2_CH42
    if (DIC2_NO32_Denit_index > 0)      a(DIC2_NO32_Denit_index)    = DIC2_NO32_Denit
    !
    if (JDIP_index > 0)                 a(JDIP_index)               = JDIP
    if (TIP_TIP2_index > 0)             a(TIP_TIP2_index)           = TIP_TIP2
    if (TIP1_Burial_index > 0)          a(TIP1_Burial_index)        = TIP1_Burial
    !if (PIP1_PIP2_index > 0)            a(PIP1_PIP2_index)          = PIP1_PIP2
    if (DIP1_DIP2_index > 0)            a(DIP1_DIP2_index)          = DIP1_DIP2
    if (TIP2_Burial_index > 0)          a(TIP2_Burial_index)        = TIP2_Burial
  ''' 

# compute sediment diagenesis derived variables

  #real(R8) :: kdpo41 TODO why is this here?
  TPOC2 = 0.0
  TPON2 = 0
  TPOP2 = 0.0
  for i in (1, 3) :
    TPOC2 = TPOC2 + POC2(i)
    TPON2 = TPON2 + PON2(i)
    TPOP2 = TPOP2 + POP2(i)

  POM2 = TPOC2 / focm2
  
  # TH2S1 and TH2S2
  H2S1 = TH2S1 / (1.0 + Css1 * kdh2s2)
  H2S2 = TH2S2 / (1.0 + Css2 * kdh2s2)
  
  # TIP1 and TIP2
  if DOX >= DOcr :
    kdpo41 = kdpo42 * d_kpo41
  else :
    kdpo41 = kdpo42 * d_kpo41**(DOX / DOcr)

  DIP1  = TIP1 / (1.0 + Css1 * kdpo41)
  DIP2  = TIP2 / (1.0 + Css2 * kdpo42)


def MatrixSolution(x1, x2, a11, a12, b1, a21, a22, b2):
    x1 = x1
    x2 = x2
    a11 = a11
    a12 = a12
    b1 = b1
    a21 = a21
    a22 = a22
    b2 = b2

    if (a11 * a22 - a12 * a21 == 0.0):
        print('Twod is singular: A11,A12,A21,A22')
        print('a11, a12, a21, a22')

    x1 = (a22 * b1 - a12 * b2) / (a11 * a22 - a12 * a21)
    x2 = (a11 * b2 - a21 * b1) / (a11 * a22 - a12 * a21)
    return x1, x2
