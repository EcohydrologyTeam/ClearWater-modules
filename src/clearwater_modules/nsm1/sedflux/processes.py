"""
File contains process to calculate sediment flux and associated dependent variables
"""

import numba
import xarray as xr
from clearwater_modules.shared.processes import arrhenius_correction
import math

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
    """Calculate vno31_tc: denitrification reaction velocity in sediment layer 1 temperature corrected (m/d).

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
def fps1(
fds1: xr.DataArray,

) -> xr.DataArray:
    """Calculate fps1: particulate fraction for H2S1 and H2S2 in layer 1 (unitless)

    Args:
      fds1: dissolved fraction for H2S1 and H2S2 in layer 1 (unitless)

    """

    return 1.0 - fds1

@numba.njit  
def fps2(
fds2: xr.DataArray,

) -> xr.DataArray:
    """Calculate fps2: particulate fraction for H2S1 and H2S2 in layer 2 (unitless)

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
    """Calculate con_nit: something nitrogen (m2/d2) #TODO define this

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

    to_return = xr.where(POCdiagenesis_part_option == 2 and t< 1.0E-10, math.sqrt(Dd_tc * SO4 * h2 / (max(roc * JC, 1.0E-10))), 
                         xr.where(POCdiagenesis_part_option==1, math.sqrt(2.0 * Dd_tc * SO4 * h2 / JCc), "NaN"))
    to_return = xr.where(to_return > h2, h2, to_return)
    return to_return

@numba.njit  
def con_sox (
vh2sd_tc: xr.DataArray,
fds1: xr.DataArray,
vh2sp_tc: xr.DataArray,
fps1: xr.DataArray,
DOX: xr.DataArray,
KSh2s: xr.DataArray,

) -> xr.DataArray:
    """Calculate con_sox: #TODO define that this is ()
    
    Args:
      fds1: dissolved fraction for H2S1 and H2S2 in layer 1 (unitless)
      vh2sd_Tc: dissolve sulfide oxidation reaction velocity in sediment layer 1 (m/d)
      vh2sp_tc: particulate sulfide oxidation reaction veolcoity in sediment layer 1 (m/d)
      fps1: particulate fraction for H2S1 and H2S2 in layer 1 (unitless)
      DOX: dissolved oxgyen concentration (mg-O2/L)
      KSh2s: sulfide oxidation normalization constant (mg-O2/L)


    """
    #TODO possibly better way to code

    return xr.where(math.isnan((vh2sd_tc * vh2sd_tc * fds1 + vh2sp_tc * vh2sp_tc * fps1) * DOX / 2.0 / KSh2s), 
              (vh2sd_tc * vh2sd_tc * fds1 + vh2sp_tc * vh2sp_tc * fps1),
              (vh2sd_tc * vh2sd_tc * fds1 + vh2sp_tc * vh2sp_tc * fps1) * DOX / 2.0 / KSh2s)

@numba.njit  
def a12_SO4 (
KL12: xr.DataArray,


) -> xr.DataArray:
    """Calculate a12_SO4: coefficents for implicit finite difference form for SO4 (a11, a12, a21, a22, b1, b2) (m/d)
    
    Args:
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)


    """
    return -KL12

@numba.njit  
def a21_SO4 (
KL12: xr.DataArray,


) -> xr.DataArray:
    """Calculate a12_SO4: coefficents for implicit finite difference form for SO4 (a11, a12, a21, a22, b1, b2) (m/d)
    
    Args:
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)


    """
    return KL12

@numba.njit  
def a22_SO4 (
POCdiagenesis_part_option: xr.DataArray,
KL12: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate a22_SO4: coefficents for implicit finite difference form for SO4 (a11, a12, a21, a22, b1, b2) (m/d)
    
    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless)
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      h2: active sediment layer (m)
      dt: timestep (d)

    """
    return xr.where(POCdiagenesis_part_option==2, 
                    xr.where(SedFlux_solution_option == 1, -KL12, -KL12 - h2/dt), "NaN")

@numba.njit  
def a12_TH2S (
POCdiagenesis_part_option: xr.DataArray,
w12: xr.DataArray,
fps2: xr.DataArray,
KL12: xr.DataArray,
fds2: xr.DataArray,

) -> xr.DataArray:
    """Calculate a12_TH2S: coefficents for implicit finite difference form for TH2S (a11, a12, a21, a22, b1, b2) (m/d)
    
    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless)
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)
      fps2: particulate fraction for H2S1 and H2S2 in layer 2 (unitless)
      fds2: dissolved fraction for H2S1 and H2S2 in layer 2 (unitless)
    """

    return xr.where(POCdiagenesis_part_option==2, -w12 * fps2 - KL12 * fds2, "NaN")

@numba.njit  
def a21_TH2S (
POCdiagenesis_part_option: xr.DataArray,
w12: xr.DataArray,
fps1: xr.DataArray,
KL12: xr.DataArray,
fds1: xr.DataArray,
vb: xr.DataArray,

) -> xr.DataArray:
    """Calculate a21_TH2S: coefficents for implicit finite difference form for TH2S (a11, a12, a21, a22, b1, b2) (m/d)
    
    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless)
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)
      fps1: particulate fraction for H2S1 and H2S2 in layer 1 (unitless)
      fds1: dissolved fraction for H2S1 and H2S2 in layer 1 (unitless)
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
    """
    
    return xr.where(POCdiagenesis_part_option==2, -w12 * fps1 - KL12 * fds1 - vb, "NaN")

@numba.njit  
def a22_TH2S (
SedFlux_solution_option: xr.DataArray,
POCdiagenesis_part_option: xr.DataArray,
vb: xr.DataArray,
a12_TH2S: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate a21_TH2S: coefficents for implicit finite difference form for TH2S (a11, a12, a21, a22, b1, b2) (m/d)
    
    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless)
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      a12_TH2S: coefficents for implicit finite difference form for TH2S (a11, a12, a21, a22, b1, b2) (m/d)
      h2: active sediment layer (m)
      dt: timestep (d)
    """
    
    return xr.where(POCdiagenesis_part_option==2, 
                    xr.where(SedFlux_solution_option==1, -a12_TH2S + vb, -a12_TH2S + vb + h2 / dt), "NaN")


"""
#TODO do I need preserve? 
TH2S1_prev = TH2S1
HSO4_prev = HSO4
SO42_prev = SO42
CH42_prev = CH42
"""
@numba.njit  
def CH4sat(
depth: xr.DataArray,
TsedC: xr.DataArray

) -> xr.DataArray:
    """Calculate CH4sat: saturated concentration of methane in oxygen equivalents (mg-O2/L)
    
    Args:
      depth: water depth (m)
      TsedC: temperature sediment (C)
    """
    
    return 100.0 * (1.0 + depth / 10.0) * 1.024**(20.0 - TsedC)

@numba.njit  
def FOxch(
Methane_solution_option: xr.DataArray,
DOX: xr.DataArray,
KsOxch: xr.DataArray,

) -> xr.DataArray:
    """Calculate FOxch: methane oxidation attenuation due to low oxygen in layer 1 (unitless)
    
    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical) (unitless)
      DOX: dissolved oxgyen concentration (mg-O2/L)
      KsOxch: half-saturation coefficient for oxygen in oxidation of methane (mg-O2/L)

    """
    
    return xr.where(Methane_solution_option ==2,
           xr.where(math.nan(DOX / (KsOxch * 2.0 + DOX)), 0.0, DOX / (KsOxch * 2.0 + DOX)))

@numba.njit  
def con_cox(
Methane_solution_option: xr.DataArray,
vch41_tc: xr.DataArray,
FOxch: xr.DataArray,

) -> xr.DataArray:
    """Calculate con_cox: #TODO define that this is ()
    
    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical) (unitless)
      DOX: dissolved oxgyen concentration (mg-O2/L)
      FOxch: methane oxidation attenuation due to low oxygen in layer 1 (unitless)
      vch41_tc: methane oxidation reaction velocity in sediment layer 1 temperature corrected (m/d)
    """
    
    return xr.where(Methane_solution_option ==2,vch41_tc * vch41_tc * FOxch, "NaN")

@numba.njit  
def a12_CH4(
Methane_solution_option: xr.DataArray,
KL12: xr.DataArray,

) -> xr.DataArray:
    """Calculate a12_CH4: coefficents for implicit finite difference form for CH4 (a11, a12, a21, a22, b1, b2) (m/d)

    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical) (unitless)
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
    """
    
    return xr.where(Methane_solution_option ==2,-KL12, "NaN")

@numba.njit  
def a21_CH4(
Methane_solution_option: xr.DataArray,
KL12: xr.DataArray,

) -> xr.DataArray:
    """Calculate a21_CH4: coefficents for implicit finite difference form for CH4 (a11, a12, a21, a22, b1, b2) (m/d)

    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical) (unitless)
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
    """
    
    return xr.where(Methane_solution_option ==2,-KL12, "NaN")

@numba.njit  
def a22_CH4(
Methane_solution_option: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
KL12: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate a22_CH4: coefficents for implicit finite difference form for CH4 (a11, a12, a21, a22, b1, b2) (m/d)

    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical) (unitless)
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      h2: active Sediment layer thickness (m)
      dt: time (d)
    """
    
    return xr.where(Methane_solution_option ==2,
                    xr.where(SedFlux_solution_option == 1, KL12, KL12 + h2/dt), "NaN") #TODO make sure the NaN is appropriate since there is not alternative 


"""

    #compute SOD
    for i in (1, int(maxit)) :
      #TNH41 and TNH42
"""

@numba.njit  
def FNH4(
KsNh4: xr.DataArray,
fd1: xr.DataArray,
TNH41: xr.DataArray,

) -> xr.DataArray:
    """Calculate FNH4: modification of nitrification reaction in layer 1 (unitless)

    Args:
      KsNh4: half-saturation ammonia constant for sediment nitrification (mg-N/L)
      THN41: total concentration NH4 dissolved layer 1 (mg-N/L)
      fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)

    """
    return xr.where(KsNh4>0, KsNh4 / (KsNh4 + fd1 * TNH41), 1)

#TODO is there a better way to do this, same calculation
@numba.njit  
def TNH41_new(
TNH41: xr.DataArray,
TNH42: xr.DataArray,
a12_TNH4: xr.DataArray,
a21_TNH4: xr.DataArray,
a22_TNH4: xr.DataArray,
b2_TNH4: xr.DataArray,
con_nit: xr.DataArray,
FNH4: xr.DataArray,
KL01: xr.DataArray,
fd1: xr.DataArray,
NH4: xr.DataArray,

) -> xr.DataArray:
    """Calculate THN41_new: newtotal concentration NH4 dissolved layer 1 (mg-N/L)

    Args:
      THN41: total concentration NH4 dissolved layer 1 (mg-N/L)
      THN42: total concentration NH4 dissolved layer 2 (mg-N/L)
      a12_TNH4: coefficents for implicit finite difference form for TNH4 (m/d)
      a12_TNH4: coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2) (m/d)
      a22_TNH4: coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2) (m/d)
      b2_TNH4: coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2) (g-N/m2/d)
      con_nit: something nitrogen (m2/d2) #TODO define this
      FNH4: modification of nitrification reaction in layer 1 (unitless)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)
      NH4: Ammonium water concentration (mg-N/L)

    """
    a11= a21_TNH4 + con_nit * FNH4 / KL01 + KL01 * fd1
    b1  = KL01 * NH4
    TNH41_new, TNH42_new = MatrixSolution(TNH41, TNH42, a11, a12_TNH4, b1, a21_TNH4, a22_TNH4, b2_TNH4)

    return max(TNH41_new,0.00)

@numba.njit  
def TNH42_new(
TNH41: xr.DataArray,
TNH42: xr.DataArray,
a12_TNH4: xr.DataArray,
a21_TNH4: xr.DataArray,
a22_TNH4: xr.DataArray,
b2_TNH4: xr.DataArray,
con_nit: xr.DataArray,
FNH4: xr.DataArray,
KL01: xr.DataArray,
fd1: xr.DataArray,
NH4: xr.DataArray,

) -> xr.DataArray:
    """Calculate THN42_new: total concentration NH4 dissolved layer 2 (mg-N/L)

    Args:
      THN41: total concentration NH4 dissolved layer 1 (mg-N/L)
      THN42: total concentration NH4 dissolved layer 2 (mg-N/L)
      a12_TNH4: coefficents for implicit finite difference form for TNH4 (m/d)
      a12_TNH4: coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2) (m/d)
      a22_TNH4: coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2) (m/d)
      b2_TNH4: coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2) (g-N/m2/d)
      con_nit: something nitrogen (m2/d2) #TODO define this
      FNH4: modification of nitrification reaction in layer 1 (unitless)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)
      NH4: Ammonium water concentration (mg-N/L)

    """
    a11= a21_TNH4 + con_nit * FNH4 / KL01 + KL01 * fd1
    b1  = KL01 * NH4
    TNH41_new, TNH42_new = MatrixSolution(TNH41, TNH42, a11, a12_TNH4, b1, a21_TNH4, a22_TNH4, b2_TNH4)

    return max(TNH42_new,0.00)

@numba.njit  
def NO31_new(
NO31: xr.DataArray,
NO32: xr.DataArray,
TNH41: xr.DataArray,
a12_NO3: xr.DataArray,
a21_NO3: xr.DataArray,
a22_NO3: xr.DataArray,
b2_NO3: xr.DataArray,
con_nit: xr.DataArray,
FNH4: xr.DataArray,
KL01: xr.DataArray,
NO3: xr.DataArray,
vno31_tc: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO31_new: new NO3 sediment layer 1 (mg-N/L)
    Args:
      THN41: total concentration NH4 dissolved layer 1 (mg-N/L)
      a12_NO3: coefficents for implicit finite difference form for NO3 (m/d)
      a12_NO3: coefficents for implicit finite difference form for NO3 (a11, a12, a21, a22, b1, b2) (m/d)
      a22_NO3: coefficents for implicit finite difference form for NO3 (a11, a12, a21, a22, b1, b2) (m/d)
      b2_NO3: coefficents for implicit finite difference form for NO3 (a11, a12, a21, a22, b1, b2) (g-N/m2/d)
      con_nit: something nitrogen (m2/d2) #TODO define this
      FNH4: modification of nitrification reaction in layer 1 (unitless)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)
      NO3: Nitrate water concentration (mg-N/L)
      vno31_tc: denitrification reaction velocity in sediment layer 1 temperature corrected (m/d)
      NO31: NO3 sediment layer 1 (mg-N/L)
      NO32: NO3 sediment layer 2 (mg-N/L)

    """
    a11 = -a21_NO3 + vno31_tc * vno31_tc / KL01 + KL01
    b1  = con_nit * FNH4 / KL01 * TNH41 + KL01 * NO3
    NO31_new, NO32_new = MatrixSolution(NO31, NO32, a11, a12_NO3, b1, a21_NO3, a22_NO3, b2_NO3)
    
    return max(NO31_new, 0.0)

@numba.njit  
def NO32_new(
NO31: xr.DataArray,
NO32: xr.DataArray,
TNH41: xr.DataArray,
a12_NO3: xr.DataArray,
a21_NO3: xr.DataArray,
a22_NO3: xr.DataArray,
b2_NO3: xr.DataArray,
con_nit: xr.DataArray,
FNH4: xr.DataArray,
KL01: xr.DataArray,
NO3: xr.DataArray,
vno31_tc: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO32_new: new NO3 sediment layer 2 (mg-N/L)

    Args:
      THN41: total concentration NH4 dissolved layer 1 (mg-N/L)
      a12_NO3: coefficents for implicit finite difference form for NO3 (m/d)
      a12_NO3: coefficents for implicit finite difference form for NO3 (a11, a12, a21, a22, b1, b2) (m/d)
      a22_NO3: coefficents for implicit finite difference form for NO3 (a11, a12, a21, a22, b1, b2) (m/d)
      b2_NO3: coefficents for implicit finite difference form for NO3 (a11, a12, a21, a22, b1, b2) (g-N/m2/d)
      con_nit: something nitrogen (m2/d2) #TODO define this
      FNH4: modification of nitrification reaction in layer 1 (unitless)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)
      NO3: Nitrate water concentration (mg-N/L)
      vno31_tc: denitrification reaction velocity in sediment layer 1 temperature corrected (m/d)
      NO31: NO3 sediment layer 1 (mg-N/L)
      NO32: NO3 sediment layer 2 (mg-N/L)

    """
    a11 = -a21_NO3 + vno31_tc * vno31_tc / KL01 + KL01
    b1  = con_nit * FNH4 / KL01 * TNH41 + KL01 * NO3
    NO31_new, NO32_new = MatrixSolution(NO31, NO32, a11, a12_NO3, b1, a21_NO3, a22_NO3, b2_NO3)
    
    return max(NO32_new, 0.0)

@numba.njit  
def JC_dn(
NO31: xr.DataArray,
rondn: xr.DataArray,
vno31_tc: xr.DataArray,
KL01: xr.DataArray,
vno32_tc: xr.DataArray,
NO32: xr.DataArray,

) -> xr.DataArray:
    """Calculate JC_dn: carbon (oxygen equivalents) consumbed by denitrification (g-O2/m2/d)

    Args:
      NO31: NO3 sediment layer 1 (mg-N/L)
      NO32: NO3 sediment layer 2 (mg-N/L)
      vno31_tc: denitrification reaction velocity in sediment layer 1 temperature corrected (m/d)
      vno32_tc: denitrification reaction velocity in sediment layer 1 temperature corrected (m/d)
      rondn: oxygen stoichiometric coeff for denitrification (g-O2/g-N)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
    """
    
    return rondn * (vno31_tc * vno31_tc / KL01 * NO31 + vno32_tc * NO32)

@numba.njit  
def JCc(
roc: xr.DataArray,
JC: xr.DataArray,
JC_dn: xr.DataArray,

) -> xr.DataArray:
    """Calculate JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)

    Args:
      roc: oxygen stoichiometric coeff for organic carbon decay (g-O2/g-C)
      JC: total sediment diagenesis flux of POC (g-C/m2/d)
      JC_dn: carbon (oxygen equivalents) consumbed by denitrification (g-O2/m2/d)
    """
    
    return max(roc * JC - JC_dn, 1.0E-10)


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

"""      
      # NO31 and NO32
      
      # SO41, SO42 and TH2S1, TH2S2
      #coefficients for SO41, SO42, H2S1 and H2S2 equations

      bx=[0]*4
      ad=[(0,0,0,0),(0,0,0,0), (0,0,0,0), (0,0,0,0)] # 4x4
      
      #coefficients for SO41 and SO42 equations    
      g=[0]*2
      h=[(0,0),(0,0)] #2x2          
"""
@numba.njit  
def HSO4_new(
POCdiagenesis_part_option: xr.DataArray, 
Dd_tc: xr.DataArray,
SO4: xr.DataArray,
h2: xr.DataArray,
JCc: xr.DataArray,
ra1: xr.DataArray,
ra2: xr.DataArray,
ra0: xr.DataArray,

) -> xr.DataArray:
    """Calculate HSO4_new: hydrogen sulfate water concentration (TODO units)

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      Dd_tc: pore-water diffusion coefficient between layer 1 and 2 temperature corrected (m2/d)
      SO4: sulfate water concentration (mg-O2/L)
      h2: h2: active Sediment layer thickness (m)
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)
      ra1: first order quadratic equation coefficent for when POCdiagenesis_part_option==2 
      ra2: second order quadratic equation coefficent for when POCdiagenesis_part_option==2
      ra0: quadratic equation coefficent for when POCdiagenesis_part_option==2
    """
    HSO4_new=xr.where(POCdiagenesis_part_option == 1, math.sqrt(2.0 * Dd_tc * SO4 * h2 / JCc), 
             xr.where(SO4 <= 0.1, 0.0, (- ra1 + math.sqrt(ra1 * ra1 - 4.0 * ra2 * ra0)) / 2.0 / ra2))
    HSO4_new=xr.where(HSO4_new > h2, h2, HSO4_new)

    return HSO4_new

def KL12SO4(
POCdiagenesis_part_option: xr.DataArray, 
Dd_tc: xr.DataArray,
SO4: xr.DataArray,
h2: xr.DataArray,
HSO4_new: xr.DataArray,
KL12: xr.DataArray,

) -> xr.DataArray:
    """Calculate KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      Dd_tc: pore-water diffusion coefficient between layer 1 and 2 temperature corrected (m2/d)
      SO4: sulfate water concentration (mg-O2/L)
      h2: active Sediment layer thickness (m)
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      HSO4_new: TODO define this

    """

    return xr.where (POCdiagenesis_part_option==1, 
                     xr.where(HSO4_new ==0.0, 1, Dd_tc / (0.5 * HSO4_new)),
                     xr.where(SO4 <= 0.1,1,
                     xr.where(HSO4_new==h2, KL12, Dd_tc / (0.5 * HSO4_new))))
"""
four equations, only solved when POCdiagenesis_part_option==1
0 = bx_1 + ad_1_1 * SO41 + ad_1_2 * SO42 + ad_1_3 * TH2S1
0 = bx_2 + ad_2_1 * SO41 + ad_2_2 * SO42 - JCc * SO42 / (SO42 + KsSO4)
0 = bx_3 + ad_3_3 * TH2S1 + ad_3_4 * TH2S2
0 = bx_4 + ad_4_3 * TH2S1 + ad_4_4 * TH2S2 + JCc * SO42 / (SO42 + KsSO4)

"""

@numba.njit  
def bx_1(
KL01: xr.DataArray, 
SO4: xr.DataArray, 
POCdiagenesis_part_option: xr.DataArray, 

) -> xr.DataArray:
    """Calculate bx_1: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      SO4: sulfate water concentration (mg-O2/L)
    """

    return xr.where(POCdiagenesis_part_option==1,KL01 * SO4,"NaN")

@numba.njit  
def bx_3(
KL01: xr.DataArray, 
H2S: xr.DataArray, 
POCdiagenesis_part_option: xr.DataArray, 

) -> xr.DataArray:
    """Calculate bx_3: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      H2S: sulfide water concentration (mg-O2/L)
    """

    return xr.where(POCdiagenesis_part_option==1,KL01 * H2S,"NaN")

@numba.njit  
def ad_1_1(
KL01: xr.DataArray, 
KL12SO4: xr.DataArray, 
POCdiagenesis_part_option: xr.DataArray, 

) -> xr.DataArray:
    """Calculate ad_1_1: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)
    """

    return xr.where(POCdiagenesis_part_option==1,KL01 - KL12SO4,"NaN")

@numba.njit  
def ad_1_2(
KL12SO4: xr.DataArray, 
POCdiagenesis_part_option: xr.DataArray, 

) -> xr.DataArray:
    """Calculate ad_1_2: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)
    """

    return xr.where(POCdiagenesis_part_option==1,KL12SO4,"NaN")

@numba.njit  
def ad_1_3(
con_sox: xr.DataArray, 
POCdiagenesis_part_option: xr.DataArray, 
KL01: xr.DataArray, 

) -> xr.DataArray:
    """Calculate ad_1_3: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      con_sox: #TODO define that this is ()
    """

    return xr.where(POCdiagenesis_part_option==1,con_sox / KL01,"NaN")

@numba.njit  
def ad_3_3(
con_sox: xr.DataArray, 
POCdiagenesis_part_option: xr.DataArray, 
KL01: xr.DataArray, 
vb: xr.DataArray,
fps1: xr.DataArray,
w12: xr.DataArray,
fds1: xr.DataArray,
KL12SO4: xr.DataArray,
 
) -> xr.DataArray:
    """Calculate ad_3_3: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      con_sox: #TODO define that this is ()
      KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      fps1: particulate fraction for H2S1 and H2S2 in layer 1 (unitless)
      fds1: dissolved fraction for H2S1 and H2S2 in layer 1 (unitless)
      w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)
    """

    return xr.where(POCdiagenesis_part_option==1,- vb - fps1 * w12 - fds1 * KL01 - fds1 * KL12SO4 - con_sox / KL01,"NaN")

@numba.njit  
def ad_3_4( 
POCdiagenesis_part_option: xr.DataArray, 
fps2: xr.DataArray,
w12: xr.DataArray,
fds2: xr.DataArray,
KL12SO4: xr.DataArray,
 
) -> xr.DataArray:
    """Calculate ad_3_4: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)
      fps2: particulate fraction for H2S1 and H2S2 in layer 2 (unitless)
      fds2: dissolved fraction for H2S1 and H2S2 in layer 2 (unitless)
      w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)
    """

    return xr.where(POCdiagenesis_part_option==1,w12 * fps2 + KL12SO4 * fds2,"NaN")

@numba.njit  
def ad_2_1( 
POCdiagenesis_part_option: xr.DataArray, 
KL12SO4: xr.DataArray,
 
) -> xr.DataArray:
    """Calculate ad_2_1: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)

    """

    return xr.where(POCdiagenesis_part_option==1, KL12SO4,"NaN")

@numba.njit  
def ad_4_3( 
POCdiagenesis_part_option: xr.DataArray, 
vb: xr.DataArray,
fps1: xr.DataArray,
w12: xr.DataArray,
fds1: xr.DataArray,
KL12SO4: xr.DataArray,
 
) -> xr.DataArray:
    """Calculate ad_4_3: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      fps1: particulate fraction for H2S1 and H2S2 in layer 1 (unitless)
      fds1: dissolved fraction for H2S1 and H2S2 in layer 1 (unitless)
      w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)

    """

    return xr.where(POCdiagenesis_part_option==1, vb + fps1 * w12 + fds1 * KL12SO4,"NaN")

@numba.njit  
def ad_2_2( 
POCdiagenesis_part_option: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,
KL12SO4: xr.DataArray,
 
) -> xr.DataArray:
    """Calculate ad_2_2: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)
      h2: active Sediment layer thickness (m)
      dt: time (d)

    """

    return xr.where(POCdiagenesis_part_option==1, 
                    xr.where(SedFlux_solution_option==1, - KL12SO4, - KL12SO4 - h2 / dt),"NaN")

@numba.njit  
def bx_2( 
POCdiagenesis_part_option: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,
SO42: xr.DataArray,
 
) -> xr.DataArray:
    """Calculate bx_2: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady) 
      SO42: SO4 in sediment in layers 2 (mg-O/L)
      h2: active Sediment layer thickness (m)
      dt: time (d)

    """

    return xr.where(POCdiagenesis_part_option==1, 
                    xr.where(SedFlux_solution_option==1, 0.0, h2 * SO42 / dt),"NaN")

@numba.njit  
def ad_4_4( 
POCdiagenesis_part_option: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,
fps2: xr.DataArray,
w12: xr.DataArray,
fds2: xr.DataArray,
KL12SO4: xr.DataArray,
vb: xr.DataArray,
 
) -> xr.DataArray:
    """Calculate ad_4_4: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady) 
      KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)
      fps2: particulate fraction for H2S1 and H2S2 in layer 2 (unitless)
      fds2: dissolved fraction for H2S1 and H2S2 in layer 2 (unitless)
      w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)
      h2: active Sediment layer thickness (m)
      dt: time (d)
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units

    """

    return xr.where(POCdiagenesis_part_option==1, 
                    xr.where(SedFlux_solution_option==1, - vb - KL12SO4 * fds2 - w12 * fps2, - vb - KL12SO4 * fds2 - w12 * fps2 - h2 / dt),"NaN")

@numba.njit  
def bx_4( 
POCdiagenesis_part_option: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,
TH2S2: xr.DataArray,
 
) -> xr.DataArray:
    """Calculate bx_2: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady) 
      TH2S2: TH2S sediment layer 2 (mg-O/L)
      h2: active Sediment layer thickness (m)
      dt: time (d)

    """

    return xr.where(POCdiagenesis_part_option==1, 
                    xr.where(SedFlux_solution_option==1, 0.0, h2 * TH2S2 / dt),"NaN")

"""
eliminate H2ST1 and H2ST2 from above equation sets, get two equations
0 = g_1 + h_1_1 * SO41 + h_1_2_ * SO42
0 = g_2 + h_2_1 * SO41 + h_2_2 * SO42 + JCc * SO42 / (SO42 + KsSO4)
"""

@numba.njit  
def g_2( 
POCdiagenesis_part_option: xr.DataArray,
bx_1: xr.DataArray,
ad_3_3: xr.DataArray,
ad_1_3: xr.DataArray,
bx_3: xr.DataArray,
ad_4_4: xr.DataArray,
ad_3_4: xr.DataArray,
ad_4_3: xr.DataArray,
bx_4: xr.DataArray,
 
) -> xr.DataArray:
    """Calculate g_2: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      bx_1, ad_3_3, ad_1_3, bx_3,ad_4_4, ad_3_4, ad_4_3, bx_4: coefficents for system of equation

    """

    return xr.where(POCdiagenesis_part_option==1, ((bx_1 * ad_3_3 - ad_1_3 * bx_3) * ad_4_4 - bx_1 * ad_3_4 * ad_4_3) / (ad_1_3 * ad_3_4) + bx_4 ,"NaN")

@numba.njit  
def g_1( 
POCdiagenesis_part_option: xr.DataArray,
g_2: xr.DataArray,
bx_2: xr.DataArray,

) -> xr.DataArray:
    """Calculate g_1: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      bx_2, g_2: coefficents for system of equation

    """

    return xr.where(POCdiagenesis_part_option==1, g_2+bx_2 ,"NaN")

@numba.njit  
def h_2_1( 
POCdiagenesis_part_option: xr.DataArray,
ad_1_1: xr.DataArray,
ad_3_3: xr.DataArray,
ad_4_4: xr.DataArray,
ad_3_4: xr.DataArray,
ad_4_3: xr.DataArray,
ad_1_3: xr.DataArray,

) -> xr.DataArray:
    """Calculate h_2_1: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      ad_1_1, ad_3_3, ad_4_4, ad_3_4, ad_4_3, ad_1_3: coefficents for system of equation

    """
    return xr.where(POCdiagenesis_part_option==1, ad_1_1 * (ad_3_3 * ad_4_4 - ad_3_4 * ad_4_3) / (ad_1_3 * ad_3_4) ,"NaN")

@numba.njit  
def h_2_2( 
POCdiagenesis_part_option: xr.DataArray,
ad_3_3: xr.DataArray,
ad_4_4: xr.DataArray,
ad_3_4: xr.DataArray,
ad_4_3: xr.DataArray,
ad_1_3: xr.DataArray,
ad_1_2: xr.DataArray,

) -> xr.DataArray:
    """Calculate h_2_2: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      ad_1_2, ad_3_3, ad_4_4, ad_3_4, ad_4_3, ad_1_3: coefficents for system of equation

    """
    return xr.where(POCdiagenesis_part_option==1, ad_1_2 * (ad_3_3 * ad_4_4 - ad_3_4 * ad_4_3) / (ad_1_3 * ad_3_4) ,"NaN")

@numba.njit  
def h_1_1( 
POCdiagenesis_part_option: xr.DataArray,
h_2_1: xr.DataArray,
ad_2_1: xr.DataArray,

) -> xr.DataArray:
    """Calculate h_1_1: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      h_2_1, ad_2_1: coefficents for system of equation

    """
    return xr.where(POCdiagenesis_part_option==1, h_2_1 + ad_2_1 ,"NaN")

@numba.njit  
def h_1_2( 
POCdiagenesis_part_option: xr.DataArray,
h_2_2: xr.DataArray,
ad_2_2: xr.DataArray,

) -> xr.DataArray:
    """Calculate h_1_1: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      h_2_2, ad_2_2: coefficents for system of equation

    """
    return xr.where(POCdiagenesis_part_option==1, h_2_2 + ad_2_2 ,"NaN")

"""
eliminate SO41 and get a quadratic equation of SO42
ra2 * SO42 * SO42 + ra1 * SO42 + ra0 = 0.0 
"""
@numba.njit
def ra0_POC_1( 
POCdiagenesis_part_option: xr.DataArray,
h_1_1: xr.DataArray,
g_2: xr.DataArray,
h_2_1: xr.DataArray, 
g_1: xr.DataArray, 
KsSO4: xr.DataArray,

) -> xr.DataArray:
    """Calculate ra0_POC_1: quadratic equation coefficent (constant)

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      h_1_1, g_2, h_2_1, g_1: coefficents for system of equation
      KsSO4: half-saturation constant for sulfate in sulfate reduction (mg-O2/L)

    """
    return xr.where(POCdiagenesis_part_option==1, (h_1_1 * g_2 - h_2_1 * g_1) * KsSO4 ,"NaN")

@numba.njit
def ra1_POC_1( 
POCdiagenesis_part_option: xr.DataArray,
h_1_1: xr.DataArray,
g_2: xr.DataArray,
h_2_1: xr.DataArray, 
g_1: xr.DataArray, 
KsSO4: xr.DataArray,
JCc: xr.DataArray,

) -> xr.DataArray:
    """Calculate ra1_POC_1: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      h_1_1, g_2, h_2_1, g_1: coefficents for system of equation
      KsSO4: half-saturation constant for sulfate in sulfate reduction (mg-O2/L)
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)
    """
    return xr.where(POCdiagenesis_part_option==1, h_1_1 * g_2 - h_2_1 * g_1 + (h_1_1 * h_2_2 - h_1_2 * h_2_1) * KsSO4 + h_1_1 * JCc,"NaN")

@numba.njit
def ra2_POC_1( 
POCdiagenesis_part_option: xr.DataArray,
h_1_1: xr.DataArray,
h_2_1: xr.DataArray, 
h_2_2: xr.DataArray, 
h_1_2: xr.DataArray, 

) -> xr.DataArray:
    """Calculate ra2_POC_1: equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      h_1_1, h_2_1, h_1_2, h_2_2: coefficents for system of equation
    """
    return xr.where(POCdiagenesis_part_option==1, h_1_1 * h_2_2 - h_1_2 * h_2_1,"NaN")

@numba.njit
def sn1( 
POCdiagenesis_part_option: xr.DataArray,
ra1_POC_1: xr.DataArray,

) -> xr.DataArray:
    """Calculate sn1: roots to quadratic equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      ra1_POC_1: equation coefficent
    """
    return xr.where(POCdiagenesis_part_option==1, 
                    xr.where(ra1_POC_1 <= 0.00, -1,1))

@numba.njit
def disc( 
POCdiagenesis_part_option: xr.DataArray,
ra0_POC_1: xr.DataArray,
ra1_POC_1: xr.DataArray,
ra2_POC_1: xr.DataArray,

) -> xr.DataArray:
    """Calculate disc: roots to quadratic equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      ra0_POC_1, ra1_POC_1, ra2_POC_1 : equation coefficent
    """
    return xr.where(POCdiagenesis_part_option==1, (- ra1_POC_1 - sn1 * math.sqrt(ra1_POC_1 * ra1_POC_1 - 4.0 * ra2_POC_1 * ra0_POC_1)) / 2.0, "NaN")

@numba.njit
def r1( 
POCdiagenesis_part_option: xr.DataArray,
ra0_POC_1: xr.DataArray,
ra1_POC_1: xr.DataArray,
ra2_POC_1: xr.DataArray,
disc: xr.DataArray,

) -> xr.DataArray:
    """Calculate r1: roots to quadratic equation 

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      ra0_POC_1, ra1_POC_1, ra2_POC_1 : equation coefficent
      disc: roots to quadratic equation
    """
    return xr.where(POCdiagenesis_part_option==1, 
                    xr.where(disc != 0.0, disc/ra2_POC_1, 
                    xr.where(ra2_POC_1==0.0, -ra0_POC_1/ra1_POC_1, -ra1_POC_1/ra2_POC_1)), "NaN")

@numba.njit
def r2( 
POCdiagenesis_part_option: xr.DataArray,
ra0_POC_1: xr.DataArray,
r1: xr.DataArray,
ra2_POC_1: xr.DataArray,
disc: xr.DataArray,

) -> xr.DataArray:
    """Calculate r2: roots to quadratic equation 

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      ra0_POC_1, ra1_POC_1, ra2_POC_1 : equation coefficent
      disc,r1: roots to quadratic equation
    """
    return xr.where(POCdiagenesis_part_option==1, 
                    xr.where(disc != 0.0, ra0_POC_1/disc, 
                    xr.where(ra2_POC_1==0.0, -r1, 0.0)), "NaN")

@numba.njit
def SO42_new( 
POCdiagenesis_part_option: xr.DataArray,
r1: xr.DataArray,
r2: xr.DataArray,
h2: xr.DataArray,
KL01: xr.DataArray,
KL12: xr.DataArray,
SO4: xr.DataArray,
CSOD_H2S: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
JCc: xr.DataArray,
SO41: xr.DataArray,
SO42: xr.DataArray,
HSO4: xr.DataArray,
dt: xr.DataArray,
a12_SO4: xr.DataArray,
a21_SO4: xr.DataArray,
a22_SO4: xr.DataArray,
HSO4_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate SO42_new: new SO4 concentration in sediment in layers 2 (mg-O/L)

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      disc,r1,r2: roots to quadratic equation
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady) 
      h2: active Sediment layer thickness (m)
      dt: time (d)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      SO4: sulfate water concentration (mg-O2/L)
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      HSO4_new: hydrogen sulfate water concentration (TODO units)
      CSOD_H2S: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)
      SO41: SO4 concentration in sediment in layers 1 (mg-O/L)
      SO42: SO4 concentration in sediment in layers (mg-O/L)
      HSO4: hydrogen sulfate water concentration (TODO units)
      a12_SO4: coefficents for implicit finite difference form for SO4 (a11, a12, a21, a22, b1, b2) (m/d)
      a21_SO4: coefficents for implicit finite difference form for SO4 (a11, a12, a21, a22, b1, b2) (m/d)
      a22_SO4: coefficents for implicit finite difference form for SO4 (a11, a12, a21, a22, b1, b2) (m/d)
    """
    a11=KL01 + KL12
    b1 = KL01 * SO4 + CSOD_H2S
    b2 = xr.where(SedFlux_solution_option == 1,JCc, JCc - SO42* HSO4 / dt)
    hold, SO42_new = MatrixSolution(SO41, SO42, a11, a12_SO4, b1, a21_SO4, a22_SO4, b2)
    SO42_new=xr.where(HSO4_new==h2, max(SO42_new,0.0), SO41 / 2.0) 
    
    return xr.where(POCdiagenesis_part_option==1, 
                    xr.where(r1<0.0, r2,r1), 
                    xr.where(SO4<=0.1, SO4, SO42_new))

@numba.njit
def SO41_new( 
POCdiagenesis_part_option: xr.DataArray,
g_1: xr.DataArray,
h_1_2: xr.DataArray,
SO42_new: xr.DataArray,
h_1_1: xr.DataArray,
h2: xr.DataArray,
KL01: xr.DataArray,
KL12: xr.DataArray,
SO4: xr.DataArray,
CSOD_H2S: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
JCc: xr.DataArray,
SO41: xr.DataArray,
SO42: xr.DataArray,
HSO4: xr.DataArray,
dt: xr.DataArray,
a12_SO4: xr.DataArray,
a21_SO4: xr.DataArray,
a22_SO4: xr.DataArray,
KL12SO4: xr.DataArray,
HSO4_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate SO41_new: new SO4 concentration in sediment in layers 1 (mg-O/L)

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      h_1_1, h_1_2, g_1: coefficents for system of equation
      SO42_new: new SO4 concentration in sediment in layers 2 (mg-O/L)
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady) 
      TH2S2: TH2S sediment layer 2 (mg-O/L)
      h2: active Sediment layer thickness (m)
      dt: time (d)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)
      SO4: sulfate water concentration (mg-O2/L)
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      HSO4_new: hydrogen sulfate water concentration (TODO units)
      CSOD_H2S: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)
      SO41: SO4 concentration in sediment in layers 1 (mg-O/L)
      SO42: SO4 concentration in sediment in layers (mg-O/L)
      HSO4: hydrogen sulfate water concentration (TODO units)
      a12_SO4: coefficents for implicit finite difference form for SO4 (a11, a12, a21, a22, b1, b2) (m/d)
      a21_SO4: coefficents for implicit finite difference form for SO4 (a11, a12, a21, a22, b1, b2) (m/d)
      a22_SO4: coefficents for implicit finite difference form for SO4 (a11, a12, a21, a22, b1, b2) (m/d)

    """
    a11=KL01 + KL12
    b1 = KL01 * SO4 + CSOD_H2S
    b2 = xr.where(SedFlux_solution_option == 1,JCc, JCc - SO42* HSO4 / dt)
    SO41_new, hold = MatrixSolution(SO41, SO42, a11, a12_SO4, b1, a21_SO4, a22_SO4, b2)
    SO41_new=xr.where(HSO4_new==h2, max(SO41_new,0.0),(KL01 * SO4 + CSOD_H2S) / (KL01 + KL12SO4 * 0.5))
   

    return xr.where(POCdiagenesis_part_option==1, - (g_1 + h_1_2 * SO42_new) / h_1_1, 
                    xr.where(SO4<=0.1, SO4, SO41_new))

@numba.njit
def TH2S1_new( 
POCdiagenesis_part_option: xr.DataArray,
bx_1: xr.DataArray,
ad_1_1: xr.DataArray,
SO41_new: xr.DataArray,
ad_1_2: xr.DataArray,
SO42_new: xr.DataArray,
ad_1_3: xr.DataArray,
a21_TH2S: xr.DataArray,
con_sox: xr.DataArray,
KL01: xr.DataArray,
fds1: xr.DataArray,
H2S: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
JCc: xr.DataArray,
HSO4_new: xr.DataArray,
h2: xr.DataArray,
TH2S2: xr.DataArray,
dt: xr.DataArray,
TH2S1: xr.DataArray,
a12_TH2S: xr.DataArray,
a22_TH2S: xr.DataArray,

) -> xr.DataArray:
    """Calculate TH2S1_new: new TH2S concentration in sediment in layers 1 (mg-O/L)

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      bx_1, ad_1_1, ad_1_2, ad_1_3,: coefficents for system of equation
      SO41_new: new SO4 concentration in sediment in layers 1 (mg-O/L)
      SO42_new: new SO4 concentration in sediment in layers 2 (mg-O/L)
      a21_TH2S: a21_TH2S: coefficents for implicit finite difference form for TH2S (a11, a12, a21, a22, b1, b2) (m/d)
      con_sox: #TODO define that this is ()
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      fds1: dissolved fraction for H2S1 and H2S2 in layer 1 (unitless)
      H2S: sulfide water concentration (mg-O2/L)
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)
      HSO4_new: hydrogen sulfate water concentration (TODO units)
      h2: active sediment layer (m)
      TH2S2: TH2S sediment layer 2 (mg-O/L)
      dt: timestep (d)
      TH2S1: TH2S sediment layer 1 (mg-O/L)
      a12_TH2S: coefficents for implicit finite difference form for TH2S (a11, a12, a21, a22, b1, b2) (m/d)
      a22_TH2S: coefficents for implicit finite difference form for TH2S (a11, a12, a21, a22, b1, b2) (m/d)
    """
    a11=-a21_TH2S + con_sox / KL01 + KL01 * fds1
    b1  = KL01 * H2S
    b2=xr.where(SedFlux_solution_option == 1, JCc * HSO4_new / h2, JCc * HSO4_new / h2 + TH2S2 * h2 / dt)
    TH2S1_new, hold = MatrixSolution(TH2S1, TH2S2, a11, a12_TH2S, b1, a21_TH2S, a22_TH2S, b2)
    TH2S1_new=max(TH2S1_new,0.0)
    return xr.where(POCdiagenesis_part_option==1, - (bx_1 + ad_1_1 * SO41_new + ad_1_2 * SO42_new) / ad_1_3, TH2S1_new)

@numba.njit
def TH2S2_new( 
POCdiagenesis_part_option: xr.DataArray,
bx_3: xr.DataArray,
ad_3_3: xr.DataArray,
TH2S1_new: xr.DataArray,
ad_3_4: xr.DataArray,
a21_TH2S: xr.DataArray,
con_sox: xr.DataArray,
KL01: xr.DataArray,
fds1: xr.DataArray,
H2S: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
JCc: xr.DataArray,
HSO4_new: xr.DataArray,
h2: xr.DataArray,
TH2S2: xr.DataArray,
dt: xr.DataArray,
TH2S1: xr.DataArray,
a12_TH2S: xr.DataArray,
a22_TH2S: xr.DataArray,

) -> xr.DataArray:
    """Calculate TH2S2_new: new TH2S concentration in sediment in layers 2 (mg-O/L)

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      bx_3, ad_3_3, ad_3_4: coefficents for system of equation
      TH2S1_new: new TH2S concentration in sediment in layers 1 (mg-O/L)
      a21_TH2S: a21_TH2S: coefficents for implicit finite difference form for TH2S (a11, a12, a21, a22, b1, b2) (m/d)
      con_sox: #TODO define that this is ()
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      fds1: dissolved fraction for H2S1 and H2S2 in layer 1 (unitless)
      H2S: sulfide water concentration (mg-O2/L)
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)
      HSO4_new: hydrogen sulfate water concentration (TODO units)
      h2: active sediment layer (m)
      TH2S2: TH2S sediment layer 2 (mg-O/L)
      dt: timestep (d)
      TH2S1: TH2S sediment layer 1 (mg-O/L)
      a12_TH2S: coefficents for implicit finite difference form for TH2S (a11, a12, a21, a22, b1, b2) (m/d)
      a22_TH2S: coefficents for implicit finite difference form for TH2S (a11, a12, a21, a22, b1, b2) (m/d)
    """
    a11=-a21_TH2S + con_sox / KL01 + KL01 * fds1
    b1 = KL01 * H2S
    b2=xr.where(SedFlux_solution_option == 1, JCc * HSO4_new / h2, JCc * HSO4_new / h2 + TH2S2 * h2 / dt)
    hold,TH2S2_new = MatrixSolution(TH2S1, TH2S2, a11, a12_TH2S, b1, a21_TH2S, a22_TH2S, b2)
    TH2S2_new=max(TH2S2_new,0.0)

    return xr.where(POCdiagenesis_part_option==1, - (bx_3 + ad_3_3 * TH2S1_new) / ad_3_4, TH2S2_new)

@numba.njit
def CSOD_H2S( 
POCdiagenesis_part_option: xr.DataArray,
con_sox: xr.DataArray,
KL01: xr.DataArray,
TH2S1_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate CSOD_H2S: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      con_sox: #TODO define that this is ()
      TH2S1_new: new TH2S concentration in sediment in layers 1 (mg-O/L)
    """
    
    return xr.where(POCdiagenesis_part_option==1, con_sox / KL01 * TH2S1_new, con_sox / KL01 * TH2S1_new)

@numba.njit
def JCc_CH4( 
POCdiagenesis_part_option: xr.DataArray,
JCc: xr.DataArray,
KsSO4: xr.DataArray,
SO42_new: xr.DataArray,
h2: xr.DataArray,
HSO4_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate JCc_CH4: carbon diagenesis flux consumed for methane formation (g-O2/m2/d)

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      KsSO4: half-saturation constant for sulfate in sulfate reduction (mg-O2/L)
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)
      SO42_new: new SO4 concentration in sediment in layers 2 (mg-O/L)
      HSO4_new: hydrogen sulfate water concentration (TODO units)
      h2: active sediment layer (m)
    """
    
    return xr.where(POCdiagenesis_part_option==1, JCc * KsSO4 / (SO42_new + KsSO4), JCc * (h2 - HSO4_new) / h2)  #TODO find this formula. See page 170 is SO42 supposed to be on top as well?

@numba.njit
def ra2( 
POCdiagenesis_part_option: xr.DataArray,
JCc: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
KL01: xr.DataArray,
h2: xr.DataArray,
SO4: xr.DataArray,
con_sox: xr.DataArray,
TH2S1: xr.DataArray,
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate ra2:  equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady) 
      h2: active Sediment layer thickness (m)
      dt: time (d)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      SO4: sulfate water concentration (mg-O2/L)
      con_sox: #TODO define that this is ()
      TH2S1: TH2S sediment layer 1 (mg-O/L)
    """
    
    return xr.where(POCdiagenesis_part_option==2, 
                    xr.where(SO4>0.1, 
                    xr.where(SedFlux_solution_option == 1, 2.0 * KL01 * JCc / h2, 2.0 * KL01 * JCc / h2 + (KL01 * SO4 + con_sox / KL01 * TH2S1) / dt), "NaN"), "NaN")  


@numba.njit
def ra1( 
POCdiagenesis_part_option: xr.DataArray,
JCc: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
KL01: xr.DataArray,
h2: xr.DataArray,
SO4: xr.DataArray,
dt: xr.DataArray,
HSO4: xr.DataArray,
SO42: xr.DataArray,
Dd_tc: xr.DataArray,

) -> xr.DataArray:
    """Calculate ra1:  equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady) 
      h2: active Sediment layer thickness (m)
      dt: time (d)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      SO4: sulfate water concentration (mg-O2/L)
      HSO4: hydrogen sulfate water concentration (TODO units)
      SO42: SO4 concentration in sediment in layers (mg-O/L)
      Dd_tc: pore-water diffusion coefficient between layer 1 and 2 temperature corrected (m2/d)
    """
    
    return xr.where(POCdiagenesis_part_option==2, 
                    xr.where(SO4>0.1, 
                    xr.where(SedFlux_solution_option == 1, 2.0 * Dd_tc * JCc / h2, 2.0 * Dd_tc * JCc / h2 - 2.0 * KL01 * HSO4 * SO42 / dt), "NaN"), "NaN")  


@numba.njit
def ra0( 
POCdiagenesis_part_option: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
KL01: xr.DataArray,
SO4: xr.DataArray,
con_sox: xr.DataArray,
TH2S1: xr.DataArray,
Dd_tc: xr.DataArray,

) -> xr.DataArray:
    """Calculate ra0:  equation coefficent

    Args:
      POCdiagenesis_part_option: method for partitioin carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady) 
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      SO4: sulfate water concentration (mg-O2/L)
      Dd_tc: pore-water diffusion coefficient between layer 1 and 2 temperature corrected (m2/d)
      con_sox: #TODO define that this is ()
      TH2S1: TH2S sediment layer 1 (mg-O/L)
    """
    
    return xr.where(POCdiagenesis_part_option==2, xr.where(SO4>0.1, xr.where(SedFlux_solution_option == 1, - 2.0 * Dd_tc * (KL01 * SO4 + con_sox / KL01 * TH2S1), - 2.0 * Dd_tc * (KL01 * SO4 + con_sox / KL01 * TH2S1)), "NaN"), "NaN")  


@numba.njit
def CSODmax( 
Methane_solution_option: xr.DataArray,
KL12: xr.DataArray,
CH4sat: xr.DataArray,
JCc_CH4: xr.DataArray,

) -> xr.DataArray:
    """Calculate CSODmax: used for analytical soluton of methane (g-O2/m2/d)

    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical)
      KL12: mass transfer velocity between the two sediment layers (m/d)
      CH4sat: saturated concentration of methane in oxygen equivalents (mg-O2/L)
      JCc_CH4: carbon diagenesis flux consumed for methane formation (g-O2/m2/d)
    """
    
    return xr.where(Methane_solution_option == 1, min(math.sqrt(2.0 * KL12 * CH4sat * JCc_CH4), JCc_CH4), "NaN")

@numba.njit
def CSOD_CH4( 
Methane_solution_option: xr.DataArray,
vch41_tc: xr.DataArray,
KL01: xr.DataArray,
CSODmax: xr.DataArray,
CH41_new: xr.DataArray,
con_cox: xr.DataArray,

) -> xr.DataArray:
    """Calculate CSOD_CH4: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)

    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical)
      vch41_tc: methane oxidation reaction velocity in sediment layer 1 temperature corrected (m/d)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      CSODmax: used for analytical soluton of methane (g-O2/m2/d)
      CH41_new: CH4 sediment layer 1 (mg-O/L)
      con_cox: #TODO define that this is ()

    """
    
    return xr.where(Methane_solution_option == 1, 
                    xr.where(vch41_tc/KL01 <100.0, CSODmax * (1.0 - 2.0 / (math.exp(vch41_tc / KL01) + math.exp(-vch41_tc / KL01))), CSODmax), con_cox / KL01 * CH41_new)

@numba.njit
def CH41_new( 
Methane_solution_option: xr.DataArray,
vch41_tc: xr.DataArray,
KL01: xr.DataArray,
con_cox: xr.DataArray,
KL12: xr.DataArray,
CH4: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
JCc_CH4: xr.DataArray,
CH41: xr.DataArray,
CH42: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,
a12_CH4: xr.DataArray,
a21_CH4: xr.DataArray,
a22_CH4: xr.DataArray,
CH42_new: xr.DataArray,
CH4sat: xr.DataArray,
CSOD_CH4: xr.DataArray

) -> xr.DataArray:
    """Calculate CH41_new: new CH4 sediment layer 1 (mg-O/L)

    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical)
      vch41_tc: methane oxidation reaction velocity in sediment layer 1 temperature corrected (m/d)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)

      con_cox: #TODO define that this is ()
      KL12: mass transfer velocity between the two sediment layers (m/d)
      CH4: methane concentration (mg-o)/L
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JCc_CH4: carbon diagenesis flux consumed for methane formation (g-O2/m2/d)
      CH41: CH4 sediment layer 1 (mg-O/L)
      CH42: CH4 sediment layer 2 (mg-O/L)
      h2: active Sediment layer thickness (m)
      dt: time (d)
      a12_CH4: coefficents for implicit finite difference form for CH4 (a11, a12, a21, a22, b1, b2) (m/d)
      a21_CH4: coefficents for implicit finite difference form for CH4 (a11, a12, a21, a22, b1, b2) (m/d)
      a22_CH4: coefficents for implicit finite difference form for CH4 (a11, a12, a21, a22, b1, b2) (m/d)
      CH42_new: new CH4 sediment layer 2 (mg-O/L)
      CH4sat: saturated concentration of methane in oxygen equivalents (mg-O2/L)

    """
    a11 = KL12 + con_cox / KL01 + KL01
    b1  = KL01 * CH4
    b2 = xr.where(SedFlux_solution_option == 1, KL01 * CH4, JCc_CH4 + CH42 * h2 / dt)
    CH41_new, hold = MatrixSolution(CH41, CH42, a11, a12_CH4, b1, a21_CH4, a22_CH4, b2)
    CH41_new=xr.where(CH42_new > CH4sat, max((b1 - a12_CH4 * CH42) / a11,0.0), max(CH41_new,0.0))
    
    return xr.where(Methane_solution_option == 2, CH41_new, 
                    xr.where( vch41_tc <= 0, 0.0, CSOD_CH4 / (vch41_tc * vch41_tc / KL01)))

@numba.njit
def CH42_new( 
Methane_solution_option: xr.DataArray,
KL01: xr.DataArray,
con_cox: xr.DataArray,
KL12: xr.DataArray,
CH4: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
JCc_CH4: xr.DataArray,
CH41: xr.DataArray,
CH42: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,
a12_CH4: xr.DataArray,
a21_CH4: xr.DataArray,
a22_CH4: xr.DataArray,
CH42_new: xr.DataArray,
CH4sat: xr.DataArray

) -> xr.DataArray:
    """Calculate CH42_new: new CH4 sediment layer 2 (mg-O/L)

    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)

      con_cox: #TODO define that this is ()
      KL12: mass transfer velocity between the two sediment layers (m/d)
      CH4: methane concentration (mg-o)/L
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JCc_CH4: carbon diagenesis flux consumed for methane formation (g-O2/m2/d)
      CH41: CH4 sediment layer 1 (mg-O/L)
      CH42: CH4 sediment layer 2 (mg-O/L)
      h2: active Sediment layer thickness (m)
      dt: time (d)
      a12_CH4: coefficents for implicit finite difference form for CH4 (a11, a12, a21, a22, b1, b2) (m/d)
      a21_CH4: coefficents for implicit finite difference form for CH4 (a11, a12, a21, a22, b1, b2) (m/d)
      a22_CH4: coefficents for implicit finite difference form for CH4 (a11, a12, a21, a22, b1, b2) (m/d)
      CH42_new: new CH4 sediment layer 2 (mg-O/L)
      CH4sat: saturated concentration of methane in oxygen equivalents (mg-O2/L)

    """
    a11 = KL12 + con_cox / KL01 + KL01
    b1  = KL01 * CH4
    b2 = xr.where(SedFlux_solution_option == 1, KL01 * CH4, JCc_CH4 + CH42 * h2 / dt)
    hold, CH42_new = MatrixSolution(CH41, CH42, a11, a12_CH4, b1, a21_CH4, a22_CH4, b2)
    CH42_new=xr.where(CH42_new > CH4sat, max(CH4sat,0.0), max(CH42_new,0.0))
    
    return xr.where(Methane_solution_option == 2, CH42_new, 0.0)

@numba.njit
def NSOD( 
KL01: xr.DataArray,
ron: xr.DataArray,
con_nit: xr.DataArray,
FNH4: xr.DataArray,
TNH41_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate NSOD: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      ron: oxygen stoichiometric coeff for nitrification        (g-O2/g-N)
      con_nit: #TODO define that this is ()
      FNH4: modification of nitrification reaction in layer 1 (unitless)
      THN41_new: newtotal concentration NH4 dissolved layer 1 (mg-N/L)


    """
    
    return ron * con_nit * FNH4 / KL01 * TNH41_new

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
def SOD_Bed( 
CSOD_CH4: xr.DataArray,
res: xr.DataArray,
DOX: xr.DataArray,
CSOD_H2S: xr.DataArray,
NSOD: xr.DataArray,
maxit: xr.DataArray,
JC: xr.DataArray,
roc: xr.DataArray,
JN: xr.DataArray,

) -> xr.DataArray:
    """Calculate SOD_Bed: SedFlux sediment oxygen demand (g-O2/m2/d)

    Args:
      NSOD: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)
      CSOD_CH4: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)
      CSOD_H2S: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)
      NSOD: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)
      JC: total sediment diagenesis flux of POC (g-C/m2/d)
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)
      JN: total sediment diagenesis flux of PON (g-N/m2/d)
    """


    for i in range(1, maxit) :
        if i == 1:
            SOD_Bed= JC * roc + 1.714 * JN
  
        KL01 = KL01(SOD_Bed,DOX) #TODO need this to update for each loop

        if(math.nan(KL01) or KL01==0.0):
            KL01 = 0.00000001
            
        SOD_Bed_old = SOD_Bed
        SOD_Bed = (CSOD_CH4 + CSOD_H2S + NSOD + SOD_Bed_old) / 2.0
        if (abs(SOD_Bed - SOD_Bed_old) / SOD_Bed * 100.0 < res): 
          exit
    
        return SOD_Bed
    return SOD_Bed

"""
Output pathways
"""

@numba.njit  
def NH41_new(
fd1: xr.DataArray,
TNH41_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH41: NH4 sediment layer 1 (mg-N/L)

    Args:
      fd1: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 1 (unitless)
      THN41_new: newtotal concentration NH4 dissolved layer 1 (mg-N/L)

    """

    return fd1 * TNH41_new

@numba.njit  
def NH42_new(
fd2: xr.DataArray,
TNH42_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH42: NH4 sediment layer 2 (mg-N/L)

    Args:
      fd2: fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer 2 (unitless)
      THN42_new: newtotal concentration NH4 dissolved layer 2 (mg-N/L)

    """

    return fd2 * TNH42_new

@numba.njit  
def JNH4(
KL01: xr.DataArray, 
NH41_new: xr.DataArray,
NH4: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH42: NH4 sediment layer 2 (mg-N/L)

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d) 
      NH41: NH4 sediment layer 1 (mg-N/L)
      NH4: Ammonia water concentration (mg-N/L)

    """

    return KL01 * (NH41_new - NH4)

@numba.njit  
def TNH41_Burial( 
vb: xr.DataArray,
TNH41_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate TNH41_Burial: burial of TNH41 in sediment layer 1 (g-N/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      THN41_new: newtotal concentration NH4 dissolved layer 1 (mg-N/L)

    """

    return vb * TNH41_new

@numba.njit  
def NH41_Nitrification( 
con_nit: xr.DataArray, 
FNH4: xr.DataArray,
KL01: xr.DataArray, 
TNH41_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH41_Nitrification: nitrification of TNH41 in sediment layer 1 (g-N/m2/d) 

    Args:
      con_nit: #TODO define that this is ()
      FNH4: modification of nitrification reaction in layer 1 (unitless)
      THN41_new: newtotal concentration NH4 dissolved layer 1 (mg-N/L)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)

    """

    return con_nit * FNH4 / KL01 * TNH41_new

@numba.njit  
def NH41_NH42( 
con_nit: xr.DataArray, 
FNH4: xr.DataArray,
KL01: xr.DataArray, 
TNH41_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH41_NH42: mass transfer between TNH41 and TNH42 in dissolved form  (g-N/m2/d)

    Args:
      con_nit: #TODO define that this is ()
      FNH4: modification of nitrification reaction in layer 1 (unitless)
      THN41_new: newtotal concentration NH4 dissolved layer 1 (mg-N/L)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)

    """

    return con_nit * FNH4 / KL01 * TNH41_new

@numba.njit  
def TNH42_Burial( 
vb: xr.DataArray,
TNH42_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate TNH42_Burial: burial of TNH41 in sediment layer 2 (g-N/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      THN42_new: newtotal concentration NH4 dissolved layer 2 (mg-N/L)

    """

    return vb * TNH42_new

@numba.njit  
def JNO3( 
KL01: xr.DataArray,
NO31_new: xr.DataArray,
NO3: xr.DataArray,

) -> xr.DataArray:
    """Calculate JNO3: sediment-water flux of nitrate (g-N/m2/d)

    Args:
      NO31_new: new NO3 sediment layer 1 (mg-N/L)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      NO3: nitrate concentation water (mg-N/L)
    """

    return KL01 * (NO31_new - NO3)

@numba.njit  
def NO31_Denit( 
KL01: xr.DataArray,
vno31_tc: xr.DataArray,
NO31_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO31_Denit: denitrification of NO31 in sediment layer 1 (g-N/m2/d)

    Args:
      NO31_new: new NO3 sediment layer 1 (mg-N/L)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      vno31_tc: denitrification reaction velocity in sediment layer 1 temperature corrected (m/d)

    """

    return vno31_tc * vno31_tc / KL01 * NO31_new

@numba.njit  
def NO31_NO32( 
KL12: xr.DataArray, 
NO32_new: xr.DataArray,
NO31_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO31_NO32: mass transfer between NO31 and NO32 (g-N/m2/d)

    Args:
      NO31_new: new NO3 sediment layer 1 (mg-N/L)
      NO32_new: new NO3 sediment layer 2 (mg-N/L)
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      
    """

    return KL12 * (NO32_new - NO31_new)

@numba.njit  
def NO32_Denit( 
vno32_tc: xr.DataArray,
NO32_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO32_Denit: denitrification of NO31 in sediment layer 2 (g-N/m2/d)

    Args:
      NO32_new: new NO3 sediment layer 2 (mg-N/L)
      vno32_tc: denitrification reaction velocity in sediment layer 2 temperature corrected (m/d)

    """

    return vno32_tc * NO32_new

@numba.njit  
def JCc_SO4( 
POCdiagenesis_part_option: xr.DataArray,
JCc: xr.DataArray,
SO42_new: xr.DataArray,
KsSO4: xr.DataArray,
HSO4_new: xr.DataArray,
h2: xr.DataArray,

) -> xr.DataArray:
    """Calculate JCc_SO4: carbon diagenesis flux consumed for sulfate reduction (g-O2/m2/d)

    Args:
      POCdiagenesis_part_option: method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth) (unitless) 
      JCc: carbon diagenesis flux corrected for denitrification (g-O2/m2/d)
      SO42_new: new SO4 concentration in sediment in layers 2 (mg-O/L)
      KsSO4: half-saturation constant for sulfate in sulfate reduction (mg-O2/L)
      HSO4_new: hydrogen sulfate water concentration (TODO units)
      h2: active sediment layer (m)

    """

    return xr.where(POCdiagenesis_part_option==1, JCc * SO42_new / (SO42_new + KsSO4), JCc * HSO4_new / h2)

@numba.njit 
def JSO4( 
KL01: xr.DataArray, 
SO41_new: xr.DataArray, 
SO4: xr.DataArray,

) -> xr.DataArray:
    """Calculate JSO4: sediment-water flux of sulfate (g-O2/m2/d)

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      SO41_new: new SO4 concentration in sediment in layers 1 (mg-O/L)
      SO4:  SO4 concentration in water (mg-O/L)
    """

    return KL01 * (SO41_new - SO4)

@numba.njit 
def SO41_SO42( 
SO41_new: xr.DataArray, 
KL12SO4: xr.DataArray, 
SO42_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate SO41_SO42: mass transfer between SO41 and SO42 (g-O2/m2/d) 
  !

    Args:
      KL12SO4: dissolved mass transfer velocity of sulfate between two layers (m/d)
      SO41_new: new SO4 concentration in sediment in layers 1 (mg-O/L)
      SO42_new: new SO4 concentration in sediment in layers 2 (mg-O/L)
    """

    return KL12SO4 * (SO42_new - SO41_new)

@numba.njit 
def H2S1_new( 
fds1: xr.DataArray,
TH2S1_new: xr.DataArray, 

) -> xr.DataArray:
    """Calculate H2S1_new: new H2S1 concentration in sediment layer 1 (mg-O/L) 
  !

    Args:
      fds1: dissolved fraction for H2S1 and H2S2 in layer 1 (unitless)
      TH2S1_new: new TH2S concentration in sediment in layers 1 (mg-O/L)

    """

    return fds1 * TH2S1_new

@numba.njit 
def H2S2_new( 
fds2: xr.DataArray,
TH2S2_new: xr.DataArray, 

) -> xr.DataArray:
    """Calculate H2S2_new: new H2S2 concentration in sediment layer 2 (mg-O/L) 
  !

    Args:
      fds2: dissolved fraction for H2S1 and H2S2 in layer 2 (unitless)
      TH2S2_new: new TH2S concentration in sediment in layers 2 (mg-O/L)

    """

    return fds2 * TH2S2_new

@numba.njit 
def JH2S( 
KL01: xr.DataArray, 
H2S1_new: xr.DataArray, 
H2S: xr.DataArray,

) -> xr.DataArray:
    """Calculate JH2S: sediment-water flux of sulfide (g-O2/m2/d)
  !

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      H2S1_new: new H2S1 concentration in sediment layer 1 (mg-O/L) 
      H2S: H2S water concentration (mg-O/L)

    """

    return KL01 * (H2S1_new - H2S)

@numba.njit 
def H2S1_Oxidation( 
KL01: xr.DataArray, 
TH2S1_new: xr.DataArray, 
con_sox: xr.DataArray,

) -> xr.DataArray:
    """Calculate H2S1_Oxidation: sulfide oxidation in sediment layer 1 (g-O2/m2/d)
  !

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      TH2S1_new: new TH2S concentration in sediment in layers 1 (mg-O/L)
      con_sox: #TODO define that this is ()

    """

    return con_sox / KL01 * TH2S1_new

@numba.njit 
def TH2S1_Burial ( 
vb: xr.DataArray, 
TH2S1_new: xr.DataArray, 

) -> xr.DataArray:
    """Calculate TH2S1_Burial: burial of H2S1 in sediment layer 1 (g-O2/m2/d)

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      TH2S1_new: new TH2S concentration in sediment in layers 1 (mg-O/L)
    """

    return vb * TH2S1_new

@numba.njit 
def H2S1_H2S2( 
KL12: xr.DataArray, 
H2S2_new: xr.DataArray, 
H2S1_new: xr.DataArray,

) -> xr.DataArray:
    """Calculate H2S1_H2S2: mass transfer between H2S1 and H2S2 in dissolved form (g-O2/m2/d)

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      H2S2_new: new H2S2 concentration in sediment layer 2 (mg-O/L) 
      H2S1_new: new H2S1 concentration in sediment layer 1 (mg-O/L) 
    """

    return KL12 * (H2S2_new - H2S1_new)

@numba.njit 
def TH2S2_Burial( 
vb: xr.DataArray,  
TH2S2_new: xr.DataArray, 

) -> xr.DataArray:
    """Calculate H2S1_H2S2: mass transfer between H2S1 and H2S2 in dissolved form (g-O2/m2/d)

    Args:
      vb: burial velocity of POM2 in bed sediment (m/d) #TODO double check units
      TH2S2_new: new TH2S concentration in sediment in layers 2 (mg-O/L)

    """

    return vb * TH2S2_new

@numba.njit 
def CH41_Oxidation( 
Methane_solution_option: xr.DataArray,  
CSOD_CH4: xr.DataArray,
con_cox: xr.DataArray,  
KL01: xr.DataArray,   
CH41_new: xr.DataArray,  

) -> xr.DataArray:
    """Calculate CH41_Oxidation: methane oxidation in sediment layer 1 (g-O2/m2/d)

    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical)  
      CSOD_CH4: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)
      con_cox: #TODO define that this is () 
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d) 
      CH41_new: new CH4 sediment layer 1 (mg-O/L)

    """

    return xr.where(Methane_solution_option==1, CSOD_CH4, con_cox / KL01 * CH41_new)

@numba.njit 
def JCH4( 
Methane_solution_option: xr.DataArray,  
CSOD_CH4: xr.DataArray,  
KL01: xr.DataArray,   
CH41_new: xr.DataArray, 
CSODmax: xr.DataArray,
CH4: xr.DataArray,

) -> xr.DataArray:
    """Calculate JCH4: sediment-water flux of methane (g-O2/m2/d)

    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical)  
      CSOD_CH4: carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d) 
      CH41_new: new CH4 sediment layer 1 (mg-O/L)
      CSODmax: used for analytical soluton of methane (g-O2/m2/d)
      CH4: CH4 concentration water (mg-O/L)

    """

    return xr.where(Methane_solution_option==1, CSODmax - CSOD_CH4, KL01 * (CH41_new - CH4))

@numba.njit 
def JCH4g( 
Methane_solution_option: xr.DataArray,  
CSODmax: xr.DataArray,
CH42_new: xr.DataArray,
CH4sat: xr.DataArray,
JCc_CH4: xr.DataArray,
JCH4: xr.DataArray,
CH41_Oxidation: xr.DataArray,
CH42: xr.DataArray,
dt: xr.DataArray,
h2: xr.DataArray,

) -> xr.DataArray:
    """Calculate JCH4g: methane loss as bubbles from sediment (g-O2/m2/d)

    Args:
      Methane_solution_option: method for solving methane concentration (1 analytical and 2 numerical)  
      CSODmax: used for analytical soluton of methane (g-O2/m2/d)
      CH42_new: new CH4 sediment layer 2 (mg-O/L)
      CH4sat: saturated concentration of methane in oxygen equivalents (mg-O2/L)
      JCc_CH4: carbon diagenesis flux consumed for methane formation (g-O2/m2/d)
      JCH4: sediment-water flux of methane (g-O2/m2/d)
      CH41_Oxidation: methane oxidation in sediment layer 1 (g-O2/m2/d)
      CH42_new: old CH4 sediment layer 2 (mg-O/L),
      h2: active Sediment layer thickness (m)
      dt: time (d)

    """

    return xr.where(Methane_solution_option==1, JCc_CH4 - CSODmax, 
                    xr.where(CH42_new == CH4sat, JCc_CH4 - JCH4 - CH41_Oxidation - (CH42_new - CH42) / dt * h2, 0.0))

@numba.njit 
def DIC1_new( 
KL01: xr.DataArray,  
DIC: xr.DataArray,  
CH41_Oxidation: xr.DataArray,  
roc: xr.DataArray,   
rcdn: xr.DataArray,   
NO31_Denit: xr.DataArray,  
KL12: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
JCc_CH4: xr.DataArray,
JCc_SO4: xr.DataArray,
NO32_Denit: xr.DataArray,
DIC1: xr.DataArray,
DIC2: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate DIC1_new: new DIC sediment layer 1 (mg-C/L)

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)  
      DIC: DIC water concentration (mg-C/L) 
      CH41_Oxidation: methane oxidation in sediment layer 1 (g-O2/m2/d) 
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)  
      rcdn: carbon stoichiometric coeff for denitrification (g-C/g-N)   
      NO31_Denit: NO31_Denit: denitrification of NO31 in sediment layer 1 (g-N/m2/d)  
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JCc_CH4: carbon diagenesis flux consumed for methane formation (g-O2/m2/d)
      JCc_SO4: carbon diagenesis flux consumed for sulfate reduction (g-O2/m2/d)
      NO32_Denit:denitrification of NO31 in sediment layer 2 (g-N/m2/d)
      DIC1: DIC sediment layer 1 (mg-C/L)
      DIC2: DIC sediment layer 1 (mg-C/L)
      h2: active sediment layer (m)
      dt: time (d)

    """
    a11= KL01 + KL12
    a12= -KL12
    b1= KL01 * DIC * 12000.0 + CH41_Oxidation / 2.0 / roc + rcdn * NO31_Denit
    a21= -KL12

    b2=xr.where(SedFlux_solution_option == 1, (JCc_CH4 / 2.0 + JCc_SO4) / roc  + rcdn * NO32_Denit, (JCc_CH4 / 2.0 + JCc_SO4) / roc + rcdn * NO32_Denit + DIC2 * h2 / dt)
    a22=xr.where(SedFlux_solution_option == 1, KL12, KL12 + h2 / dt)
    DIC1_new, hold = MatrixSolution(DIC1, DIC2, a11, a12, b1, a21, a22, b2)

    return max(DIC1_new,0.0)

@numba.njit 
def DIC2_new( 
KL01: xr.DataArray,  
DIC: xr.DataArray,  
CH41_Oxidation: xr.DataArray,  
roc: xr.DataArray,   
rcdn: xr.DataArray,   
NO31_Denit: xr.DataArray,  
KL12: xr.DataArray,
SedFlux_solution_option: xr.DataArray,
JCc_CH4: xr.DataArray,
JCc_SO4: xr.DataArray,
NO32_Denit: xr.DataArray,
DIC1: xr.DataArray,
DIC2: xr.DataArray,
h2: xr.DataArray,
dt: xr.DataArray,

) -> xr.DataArray:
    """Calculate DIC2_new: new DIC sediment layer 2 (mg-C/L)

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)  
      DIC: DIC water concentration (mg-C/L) 
      CH41_Oxidation: methane oxidation in sediment layer 1 (g-O2/m2/d) 
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)  
      rcdn: carbon stoichiometric coeff for denitrification (g-C/g-N)   
      NO31_Denit: NO31_Denit: denitrification of NO31 in sediment layer 1 (g-N/m2/d)  
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JCc_CH4: carbon diagenesis flux consumed for methane formation (g-O2/m2/d)
      JCc_SO4: carbon diagenesis flux consumed for sulfate reduction (g-O2/m2/d)
      NO32_Denit:denitrification of NO31 in sediment layer 2 (g-N/m2/d)
      DIC1: DIC sediment layer 1 (mg-C/L)
      DIC2: DIC sediment layer 1 (mg-C/L)
      h2: active sediment layer (m)
      dt: time (d)
      
    """
    a11= KL01 + KL12
    a12= -KL12
    b1= KL01 * DIC * 12000.0 + CH41_Oxidation / 2.0 / roc + rcdn * NO31_Denit
    a21= -KL12

    b2=xr.where(SedFlux_solution_option == 1, (JCc_CH4 / 2.0 + JCc_SO4) / roc  + rcdn * NO32_Denit, (JCc_CH4 / 2.0 + JCc_SO4) / roc + rcdn * NO32_Denit + DIC2 * h2 / dt)
    a22=xr.where(SedFlux_solution_option == 1, KL12, KL12 + h2 / dt)
    hold, DIC2_new = MatrixSolution(DIC1, DIC2, a11, a12, b1, a21, a22, b2)

    return max(DIC2_new,0.0)

@numba.njit 
def JDIC( 
KL01: xr.DataArray,  
DIC: xr.DataArray,  
DIC1_new: xr.DataArray,


) -> xr.DataArray:
    """Calculate JDIC: sediment-water flux of dissolved inorganic carbon (g-C/m2/d)

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)  
      DIC: DIC water concentration (mg-C/L) 
      DIC1_new: new DIC sediment layer 1 (mg-C/L)
    """

    return KL01 * (DIC1_new - DIC * 12000.0)

@numba.njit 
def DIC1_CH41_Oxidation( 
CH41_Oxidation: xr.DataArray,  
roc: xr.DataArray,  

) -> xr.DataArray:
    """Calculate DIC1_CH41_Oxidation: DIC1 produced by CH41 oxidation in sediment layer 1 (g-C/m2/d)

    Args:
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)  
      CH41_Oxidation: methane oxidation in sediment layer 1 (g-O2/m2/d)
    """

    return CH41_Oxidation / 2.0 / roc

@numba.njit 
def DIC1_NO31_Denit( 
rcdn: xr.DataArray,  
NO31_Denit: xr.DataArray,  

) -> xr.DataArray:
    """Calculate DIC1_NO31_Denit: DIC1 produced by NO31 denitrification in sediment layer 1 (g-C/m2/d) 

    Args:
      rcdn: carbon stoichiometric coeff for denitrification (g-C/g-N)   
      NO31_Denit: denitrification of NO31 in sediment layer 1 (g-N/m2/d)

    """

    return rcdn * NO31_Denit

@numba.njit 
def DIC1_DIC2(  
KL12: xr.DataArray,
DIC2_new: xr.DataArray,   
DIC1_new: xr.DataArray,  

) -> xr.DataArray:
    """Calculate DIC1_DIC2: mass transfer between DIC1 and DIC2 in dissolved form (g-C/m2/d)

    Args:
      KL12: dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)
      DIC1_new: new DIC sediment layer 1 (mg-C/L)
      DIC2_new: new DIC sediment layer 2 (mg-C/L)

    """

    return KL12 * (DIC2_new - DIC1_new)

@numba.njit 
def DIC2_POC2_SO42(  
JCc_SO4: xr.DataArray, 
roc: xr.DataArray,

) -> xr.DataArray:
    """Calculate DIC2_POC2_SO42: DIC2 produced by sulfate reduction in sediment layer 2 (g-C/m2/d)

    Args:
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)  
      JCc_SO4: carbon diagenesis flux consumed for sulfate reduction (g-O2/m2/d)

    """

    return JCc_SO4 / roc

@numba.njit 
def DIC2_CH42(  
JCc_CH4: xr.DataArray, 
roc: xr.DataArray,

) -> xr.DataArray:
    """Calculate DIC2_CH42: DIC2 produced by mathene formation in sediment layer 2 (g-C/m2/d)

    Args:
      roc: oxygen stoichiometric coefficent for organic carbon decay (g-O2/g-C)  
      JCc_CH4: carbon diagenesis flux consumed for methane formation (g-O2/m2/d)

    """

    return JCc_CH4 / 2.0 / roc

@numba.njit 
def DIC2_NO32_Denit(  
rcdn: xr.DataArray, 
NO32_Denit: xr.DataArray,

) -> xr.DataArray:
    """Calculate DIC2_NO32_Denit: DIC2 produced by mathene formation in sediment layer 2 (g-C/m2/d)

    Args:
      rcdn: carbon stoichiometric coeff for denitrification (g-C/g-N) 
      NO32_Denit : denitrification of NO31 in sediment layer 2 (g-N/m2/d)

    """

    return rcdn * NO32_Denit

@numba.njit 
def kdpo41(  
DOX: xr.DataArray,
DOcr: xr.DataArray, 
kdpo42: xr.DataArray, 
d_kpo41: xr.DataArray, 

) -> xr.DataArray:
    """Calculate kdpo41: partition coefficient for inorganic P in sediment layer 1 (L/kg)

    Args:
      DOX: dissolved oxygen concentration (mg-O/L)
      DOcr: critical oxygen concentration for incremental phosphate sorption  (mg-O2/L) [avoid to repeat with DOC] 
      kdpo42: partition coefficient for inorganic P in sediment layer 2 (L/kg)
      d_kpo41: factor that increases the aerobic layer phosphate partition coefficient (unitless) 

    """

    return xr.where(DOX >= DOcr, kdpo42 * d_kpo41, kdpo42 * d_kpo41**(DOX / DOcr))

@numba.njit 
def TIP1_new(  
SedFlux_solution_option: xr.DataArray,
JP: xr.DataArray, 
TIP: xr.DataArray,
fdp: xr.DataArray,
vs: xr.DataArray,
vb: xr.DataArray,
h2: xr.DataArray,
TIP1: xr.DataArray,
TIP2: xr.DataArray,
dt: xr.DataArray,
KL01: xr.DataArray,
Css2: xr.DataArray,
kdpo42: xr.DataArray,
Css1: xr.DataArray,
kdpo41: xr.DataArray,
w12: xr.DataArray,
KL12: xr.DataArray,

) -> xr.DataArray:
    """Calculate TIP1_new: new TIP sediment layer 1 (mg-P/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JP:  total sediment diagenesis flux of POP  (g-P/m2/d) 
      TIP: TIP water concentration (mg-P/L)
      fdp: fraction of dissolved phosphrous (unitless)
      vs: sediment settling velocity (m/d)
      vb: burial velocity of POM2 in bed sediment (m/d)
      h2: active sediment layer (m)
      TIP1: TIP sediment layer 1 (mg-P/L)
      TIP2: TIP sediment layer 2 (mg-P/L)
      dt: time step (d)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      Css2: solids concentration in sediment layer 2 (kg/L)
      kdpo42: partition coefficient for inorganic P in sediment layer 2 (L/kg)
      Css1: solids concentration in sediment layer 1 (kg/L)
      kdpo41: kdpo41: partition coefficient for inorganic P in sediment layer 1 (L/kg)
      w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)
      KL12: Dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)

    """
    fd1 = 1.0 / (1.0 + Css1 * kdpo41)
    fd2 = 1.0 / (1.0 + Css2 * kdpo42)
    fp2 = 1.0 - fd2
    fp1 = 1.0 - fd1

    a21 = -w12 * fp1 - KL12 * fd1 - vb
    a11 = -a21 + KL01 * fd1
    a12 = -w12 * fp2 - KL12 * fd2
    b1  = KL01 * fdp * TIP

    a22 = xr.where (SedFlux_solution_option == 1,-a12 + vb,-a12 + vb + h2 / dt)
    b2= xr.where (SedFlux_solution_option == 1,JP + TIP * (1.0 - fdp) * vs,JP + TIP * (1.0 - fdp) * vs + h2 * TIP2 / dt)
    TIP1_new, hold = MatrixSolution(TIP1, TIP2, a11, a12, b1, a21, a22, b2)

    return max(TIP1_new,0.0)

@numba.njit 
def TIP2_new(  
SedFlux_solution_option: xr.DataArray,
JP: xr.DataArray, 
TIP: xr.DataArray,
fdp: xr.DataArray,
vs: xr.DataArray,
vb: xr.DataArray,
h2: xr.DataArray,
TIP1: xr.DataArray,
TIP2: xr.DataArray,
dt: xr.DataArray,
KL01: xr.DataArray,
Css2: xr.DataArray,
kdpo42: xr.DataArray,
Css1: xr.DataArray,
kdpo41: xr.DataArray,
w12: xr.DataArray,
KL12: xr.DataArray,

) -> xr.DataArray:
    """Calculate TIP2_new: new TIP sediment layer 2 (mg-P/L)

    Args:
      SedFlux_solution_option: numerical method (1 steady, 2 unsteady)
      JP:  total sediment diagenesis flux of POP  (g-P/m2/d) 
      TIP: TIP water concentration (mg-P/L)
      fdp: fraction of dissolved phosphrous (unitless)
      vs: sediment settling velocity (m/d)
      vb: burial velocity of POM2 in bed sediment (m/d)
      h2: active sediment layer (m)
      TIP1: TIP sediment layer 1 (mg-P/L)
      TIP2: TIP sediment layer 2 (mg-P/L)
      dt: time step (d)
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      Css2: solids concentration in sediment layer 2 (kg/L)
      kdpo42: partition coefficient for inorganic P in sediment layer 2 (L/kg)
      Css1: solids concentration in sediment layer 1 (kg/L)
      kdpo41: partition coefficient for inorganic P in sediment layer 1 (L/kg)
      w12: Partical mixing transfer velocity: transfer for NH4, H2S, and PIP between layer 1 and 2 (m/d)
      KL12: Dissolved and particulate phase mixing coefficient between layer 1 and layer 2 (m/d)

    """
    fd1 = 1.0 / (1.0 + Css1 * kdpo41)
    fd2 = 1.0 / (1.0 + Css2 * kdpo42)
    fp2 = 1.0 - fd2
    fp1 = 1.0 - fd1

    a21 = -w12 * fp1 - KL12 * fd1 - vb
    a11 = -a21 + KL01 * fd1
    a12 = -w12 * fp2 - KL12 * fd2
    b1  = KL01 * fdp * TIP

    a22 = xr.where (SedFlux_solution_option == 1,-a12 + vb,-a12 + vb + h2 / dt)
    b2= xr.where (SedFlux_solution_option == 1,JP + TIP * (1.0 - fdp) * vs,JP + TIP * (1.0 - fdp) * vs + h2 * TIP2 / dt)
    hold, TIP2_new = MatrixSolution(TIP1, TIP2, a11, a12, b1, a21, a22, b2)

    return max(TIP2_new,0.0)

@numba.njit 
def DIP1_new(  
TIP1_new: xr.DataArray,
Css1: xr.DataArray,
kdpo41: xr.DataArray,

) -> xr.DataArray:
    """Calculate DIP1_new: new DIP sediment layer 1 (mg-P/L)

    Args:
      TIP1_new: new TIP sediment layer 1 (mg-P/L)
      Css1: solids concentration in sediment layer 1 (kg/L)
      kdpo41: partition coefficient for inorganic P in sediment layer 1 (L/kg)

    """

    return TIP1_new * (1.0 / (1.0 + Css1 * kdpo41))

@numba.njit 
def DIP2_new(  
TIP2_new: xr.DataArray,
Css2: xr.DataArray,
kdpo42: xr.DataArray,

) -> xr.DataArray:
    """Calculate DIP2_new: new DIP sediment layer 2 (mg-P/L)

    Args:
      TIP2_new: new TIP sediment layer 2 (mg-P/L)
      Css2: solids concentration in sediment layer 2 (kg/L)
      kdpo42: partition coefficient for inorganic P in sediment layer 2 (L/kg)

    """

    return TIP2_new * (1.0 / (1.0 + Css2 * kdpo42))

@numba.njit 
def JDIP(  
KL01: xr.DataArray,
DIP1_new: xr.DataArray,
fdp: xr.DataArray, 
TIP: xr.DataArray,

) -> xr.DataArray:
    """Calculate JDIP: sediment-water flux of phosphate (g-P/m2/d)

    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      DIP1_new: new DIP sediment layer 1 (mg-P/L)
      fdp: fraction of dissolved phosphrous (unitless)
      TIP: TIP water concentration (mg-P/L)

    """

    return KL01 * (DIP1_new - fdp * TIP)

@numba.njit 
def TIP_TIP2(  
vs: xr.DataArray, 
fdp: xr.DataArray, 
TIP: xr.DataArray,

) -> xr.DataArray:
    """Calculate TIP_TIP2: settling of PIP of water column into PIP2 in layer 2 (g-P/m2/d)

    Args:
      vs: sediment settling velocity (m/d)
      fdp: fraction of dissolved phosphrous (unitless)
      TIP: TIP water concentration (mg-P/L)

    """

    return TIP * (1.0 - fdp) * vs

@numba.njit 
def TIP1_Burial(  
vb: xr.DataArray, 
TIP1_new: xr.DataArray, 

) -> xr.DataArray:
    """Calculate TIP1_Burial: burial of TIP1 in sediment layer 1 (g-P/m2/d)
    Args:
      vb: burial velocity of POM2 in bed sediment (m/d)
      TIP1_new: new TIP sediment layer 1 (mg-P/L)
    """

    return vb * TIP1_new

@numba.njit 
def DIP1_DIP2(  
KL12: xr.DataArray, 
DIP2_new: xr.DataArray, 
DIP1_new: xr.DataArray, 

) -> xr.DataArray:
    """Calculate DIP1_DIP2: mass transfer between TIP1 and TPO42 in dissolved form   (g-P/m2/d)
    Args:
      KL01: mass transfer velocity between overlying water and the aerobic layer  (m/d)
      DIP1_new: new DIP sediment layer 1 (mg-P/L)
      DIP2_new: new DIP sediment layer 2 (mg-P/L)

    """

    return KL12 * (DIP2_new - DIP1_new)

@numba.njit 
def TIP2_Burial(  
vb: xr.DataArray, 
TIP2_new: xr.DataArray, 

) -> xr.DataArray:
    """Calculate TIP2_Burial: burial of TIP2 in sediment layer 2 (g-P/m2/d)
    Args:
      vb: burial velocity of POM2 in bed sediment (m/d)
      TIP2_new: new TIP sediment layer 2 (mg-P/L)
    """

    return vb * TIP2_new

@numba.njit 
def TPOC2(  
POC2_1_new: xr.DataArray, 
POC2_2_new: xr.DataArray, 
POC2_3_new: xr.DataArray, 

) -> xr.DataArray:
    """Calculate TPOC2: total sediment POC (mg-C/L)

    Args:
      POC2_1_new: POC G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)
      POC2_2_new: POC G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)
      POC2_3_new: POC G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-C/L)
    """

    return POC2_1_new + POC2_2_new + POC2_3_new

@numba.njit 
def TPON2(  
PON2_1_new: xr.DataArray, 
PON2_2_new: xr.DataArray, 
PON2_3_new: xr.DataArray, 

) -> xr.DataArray:
    """Calculate TPON2: total sediment PON (mg-N/L)

    Args:
      PON2_1_new: PON G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)
      PON2_2_new: PON G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)
      PON2_3_new: PON G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-N/L)
    """

    return PON2_1_new + PON2_2_new + PON2_3_new

@numba.njit 
def TPOP2(  
POP2_1_new: xr.DataArray, 
POP2_2_new: xr.DataArray, 
POP2_3_new: xr.DataArray, 

) -> xr.DataArray:
    """Calculate TPOP2: total sediment POP (mg-P/L)

    Args:
      POP2_1_new: POP G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)
      POP2_2_new: POP G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)
      POP2_3_new: POP G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)
    """

    return POP2_1_new + POP2_2_new + POP2_3_new

@numba.njit 
def POM2(  
TPOC2_new: xr.DataArray, 
focm2: xr.DataArray, 

) -> xr.DataArray:
    """Calculate POM2: particulate organic matter (mg/L)

    Args:
      POP2_1_new: POP G1 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)
      POP2_2_new: POP G2 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)
      POP2_3_new: POP G3 in the second layer. Able for diagenesis and depends on depositional flux above (mg-P/L)
    """

    return TPOC2_new / focm2
"""  
    
    # TIP1 and TIP2

    fd1 = 1.0 / (1.0 + Css1 * kdpo41)
    fd2 = 1.0 / (1.0 + Css2 * kdpo42)
    fp1 = 1.0 - fd1
    fp2 = 1.0 - fd2

  
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
"""