from email.feedparser import FeedParser
import math
from numba import types
from _temp_correction import TempCorrection

class ApGrowth:
    def __init__(self, Ap, TwaterC, PAR, lambda0, depth, NH4, NO3, TIP, fdp, mu_max, KsN, KsP, KL, use_NH4, use_NO3, use_TIP, light_limitation_option, growth_rate_option) :
        self.Ap = Ap
        self.TwaterC = TwaterC
        self.mu_max = mu_max
        self.PAR=PAR
        self.lambda0=lambda0
        self.depth=depth
        self.KL=KL
        self.NH4=NH4
        self.NO3=NO3
        self.KsN=KsN
        self.KsP=KsP
        self.fdp=fdp
        self.TIP=TIP
        self.use_NO3=use_NO3
        self.use_NH4=use_NH4
        self.use_TIP=use_TIP
        self.light_limitation_option=light_limitation_option
        self.growth_rate_option=growth_rate_option
    
    def Calculation(self) :
        mu_max_tc = TempCorrection(self.mu_max, 1.047).arrhenius_correction(self.TwaterC)           # Maximum algal growth rate [1/d]
        
        sqrt1 = 0.0
        sqrt2 = 0.0

        # Depth averaged light function
        KEXT = (self.lambda0 * self.depth)    # lambda is light attenuation coefficient (1/m). depth is depth from water surface (m) [unitless] TODO: other depth not initalized method

        # (1) Algal light limitation (FL)

        if (self.Ap <= 0.0 or KEXT <= 0.0 or self.PAR <= 0.0):
            # After sunset or if there is no algae present
            FL = 0.0                                                                        # light limiting factor for algal growth [unitless]
        elif self.light_limitation_option == 1:
            # Half-saturation formulation
            FL = (1.0 / KEXT) * math.log((self.KL + self.PAR) / (self.KL + self.PAR * math.exp(-KEXT)))         
        elif self.light_limitation_option == 2:
            # Smith's model
            if abs(self.KL) < 1.0E-10:
                FL = 1.0                                                                    
            else:
                sqrt1 = (1.0 + (self.PAR / self.KL)**2.0)**0.5                                        
                sqrt2 = (1.0 + (self.PAR * math.exp(-KEXT) / self.KL)**2.0)**0.5
                FL = (1.0 / KEXT) * math.log((self.PAR / self.KL + sqrt1) / (self.PAR * math.exp(-KEXT) / self.KL + sqrt2))
        elif self.light_limitation_option == 3:
            # Steele's model
            if abs(self.KL) < 1.0E-10:
                FL = 0.0
            else:
                FL = (2.718/KEXT) * (math.exp(-self.PAR/self.KL * math.exp(-KEXT)) - math.exp(-self.PAR/self.KL))

        # Limit factor to between 0.0 and 1.0.
        # This should never happen, but it would be a mess if it did.
        if FL > 1.0:
            FL = 1.0
        if FL < 0.0:
            FL = 0.0
        # (2) Algal nitrogen limitation (FN)
        # KsN = Michaelis-Menton half-saturation constant (mg N/L) relating inorganic N to algal growth
        if self.use_NH4 or self.use_NO3:

            FN = (self.NH4 + self.NO3) / (self.KsN + self.NH4 + self.NO3)                  # [unitless]
            if math.isnan(FN):
                FN = 0.0
            if FN > 1.0:
                FN = 1.0
        else:
            FN = 1.0
        # (3) Algal phosphorous limitation (FP)
        # PO4 = Dissolved (inorganic) phosphorous (mg-P/L)
        # KsP = Michaelis-Menton half-saturation constant (mg-P/L) relating inorganic P to algal growth

        if self.use_TIP:
            FP = self.fdp * self.TIP / (self.KsP + self.fdp * self.TIP)              # [unitless]
            if math.isnan(FP):
                FP = 0.0
            if FP > 1.0:
                FP = 1.0
        else:
            FP = 1.0

        # Algal growth rate with three options
        # (a) Multiplicative (b) Limiting nutrient (c) Harmonic Mean
        if self.growth_rate_option == 1:
            # (a) Multiplicative (day-1)
            mu = mu_max_tc * FL * FP * FN                   # [1/d]
        elif self.growth_rate_option == 2:
            # (b) Limiting nutrient (day-1)
            mu = mu_max_tc * FL * min(FP, FN)
        elif self.growth_rate_option == 3:
            # (c) Harmonic Mean Option (day-1)
            if FN == 0.0 or FP == 0.0:
                mu = 0.0
            else:
                mu = mu_max_tc * FL * 2.0 / (1.0 / FN + 1.0 / FP)

        # Algal growth
        self.ApGrowth = mu * self.Ap                      # [ug-Chla/L/d]
    
        return self.ApGrowth
