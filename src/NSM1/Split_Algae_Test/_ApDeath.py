import math
from numba import types
from _temp_correction import TempCorrection

class ApDeath:
    def __init__(self, Ap, TwaterC, kdp=0.15) :
        self.kdp=kdp
        self.TwaterC = TwaterC
        self.Ap = Ap
    
    def Calculation(self) :
        kdp_tc = TempCorrection(self.kdp, 1.047).arrhenius_correction(self.TwaterC)                 # algae mortality rate [1/d]
        self.ApDeath = kdp_tc * self.Ap                                                             # [ug-Chla/L/d]
    
        return self.ApDeath
    
        print (self.kdp)