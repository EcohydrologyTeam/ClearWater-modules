import math
from numba import types
from _temp_correction import TempCorrection

class ApRespiration:
    def __init__(self, krp, TwaterC, Ap) :
        self.krp=krp
        self.TwaterC = TwaterC
        self.Ap = Ap
    
    def Calculation(self) :
        krp_tc = TempCorrection(self.krp, 1.047).arrhenius_correction(self.TwaterC)                 # algae respiration rate [1/d]
        self.ApRespiration = krp_tc * self.Ap             # [ug-Chla/L/d]
                                                        # [ug-Chla/L/d]
    
        return self.ApRespiration
