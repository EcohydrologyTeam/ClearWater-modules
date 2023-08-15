'''
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Alkalinity Kinetics
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

Initial Version: June 12, 2021
'''

import math
from _temp_correction import TempCorrection
from collections import OrderedDict

class Alkalinity:

    def __init__(self, global_module_choices, global_vars, algae_pathways, balgae_pathways, nitrogen_pathways, alkalinity_constant_changes, ):
        self.global_module_choices = global_module_choices
        self.global_vars=global_vars
        self.algae_pathways = algae_pathways
        self.balgae_pathways = balgae_pathways
        self.nitrogen_pathways = nitrogen_pathways
        self.alkalinity_constant_changes = alkalinity_constant_changes

        self.alkalinity_constant=OrderedDict()
        self.alkalinity_constant = {
            'ralkca' : 14/106/12/1000, #translating algal and balgal growht into Alk if NH2 is the N source (eq/ug-chla)
            'ralkcn' : 18/106/12/1000, #ratio translating algal and balgal growth into Alk if NO3 is the N source (eq/ug-Chla)
            'ralkn' : 2/14/1000, #nitrification
            'ralkden' : 4/14/1000, # denitrification
            'pH_solution' : 2,
            'imax' : 13, # maximum iteration number for computing pH
            'es' : 0.003, # maximum relative error for computing pH
        }
        
        for key in self.alkalinity_constant_changes.keys() :
            if key in self.alkalinity_constant:
                self.alkalinity_constant[key] = self.alkalinity_constant_changes[key]
        
        #ralkaa ratio translating algal growth into Alk if NH4 is the N source (eq ug-Chla -1)
        #ralkan # ratio translating algal growth into Alk if NO3 is the N source (eq ug-Chla -1)
        #ralkba # ratio translating benthic algae growth into Alk if NH4 is the N source (eq mg-D -1)
        #ralkbn # ratio translating benthic algae growth into Alk if NO3 is the N source (eq mg-D -1)
        #ralkn # ratio translating NH4 nitrification into Alk (eq mg-N -1)
        #ralkden # ratio translating NO3 denitrification into Alk (eq mg-N -1)

    def Calculations_Alk (self) :

        # 50000.0 is used to convert eq H+/L to mgCaCO3/L

        print("Calculating change in Alkalinity")
        #Algae
        if self.global_module_choices['use_Algae']:
            Alk_ApGrowth= (self.alkalinity_constant['ralkca'] * self.nitrogen_pathways['ApUptakeFr_NH4'] - self.alkalinity_constant['ralkcn'] * (1.0 - self.nitrogen_pathways['ApUptakeFr_NH4'])) * self.algae_pathways['rca'] * self.algae_pathways['ApGrowth'] * 50000.0
            Alk_ApRespiration = self.alkalinity_constant['ralkca'] * self.algae_pathways['rca'] * self.algae_pathways['ApRespiration'] * 50000.0 
        else :
            Alk_ApGrowth = 0.0
            Alk_ApRespiration = 0.0

        #Benthic Algae
        if self.global_module_choices['use_BAlgae'] :
            Alk_AbGrowth = self.balgae_pathways['Fb'] * (self.alkalinity_constant['ralkca'] * self.nitrogen_pathways['AbUptakeFr_NH4'] - self.alkalinity_constant['ralkcn'] * (1.0 - self.nitrogen_pathways['AbUptakeFr_NH4'])) * self.balgae_pathways['rcb'] * self.balgae_pathways['AbGrowth'] / self.global_vars['depth'] * 50000.0
            Alk_AbRespiration = (self.balgae_pathways['Fb'] * self.alkalinity_constant['ralkca'] * self.balgae_pathways['rcb'] * self.balgae_pathways['AbRespiration'] / self.global_vars['depth']) * 50000.0
        else :
            Alk_AbGrowth = 0.0
            Alk_AbRespiration = 0.0

        # Nitrification TODO Check not the same formula as in the text
        if self.global_module_choices['use_NH4'] :
            Alk_Nitrification   = self.alkalinity_constant['ralkn'] * self.nitrogen_pathways['NH4_Nitrification'] * 50000.0
        else :
            Alk_Nitrification   = 0.0

        # Denitrification 
        if self.global_module_choices['use_NO3'] :
            Alk_Denit = self.alkalinity_constant['ralkden'] * self.nitrogen_pathways['NO3_Denit'] * 50000.0 #mg-CaCO3/L/d
        else :
            Alk_Denit = 0.0

        # Kinetic rate, mgCaCO3/L/day
        dAlkdt = Alk_ApGrowth + Alk_ApRespiration - Alk_Nitrification + Alk_Denit - Alk_AbGrowth + Alk_AbRespiration

        print ("dAlkdt", dAlkdt)
    
    #def Calculation_pH (self) :
        TwaterK = self.global_vars['TwaterC'] + 273.15
        self.KW = 10**((-4787.3 / TwaterK) - 7.1321 * math.log10(TwaterK) - 0.010365 * TwaterK + 22.80)
        self.K1 = 10**(-356.3094 - 0.06091964 * TwaterK + (21834.37 / TwaterK) + 126.8339 * math.log10(TwaterK) - 1684915.0 / (TwaterK**2))
        self.K2 = 10**(-107.8871 - 0.03252849 * TwaterK + 5151.79 / TwaterK + 38.92561 * math.log10(TwaterK) - 563713.9 / (TwaterK**2))
        self.hh = 1

        if self.alkalinity_constant['pH_solution'] == 1:
            self.pH=self.NewtonRaphson()

        elif self.alkalinity_constant['pH_solution'] == 2:
            self.pH=self.Bisection()
        
        print("pH",self.pH)

        return dAlkdt, self.pH
    
    def Function (self,x) :
        hh = 10**(-x)
        f = ((self.K1 * hh + 2.0 * self.K1 * self.K2) / (hh * hh + self.K1 * hh + self.K1 * self.K2)) * self.global_vars['DIC'] + (self.KW / hh) - hh - self.global_vars['Alk'] / 50000.0
        return f
    
    def NewtonRaphson (self) :
        pH_new = 1
        ea = 1
    
        # Theory of this method xn+1 = xn - f(xn)/f(xn)'
    
        pH = self.global_vars['pH']
        i = 0   # iteration count
    
        while ea > self.alkalinity_constant['es'] : 
            i = i + 1
            hh = 10**(-pH)
            original = self.Function(pH)
            der = (math.log(10.0) * (self.K1 * hh * self.global_vars['DIC'] * (hh * hh + 4 * self.K2 * hh + self.K1 * self.K2) / (hh * hh + self.K1 * hh + self.K1 * self.K2)**2 + self.KW / hh + hh))
            pH_new = pH - self.Function(pH) / (math.log(10.0) * (self.K1 * hh * self.global_vars['DIC'] * (hh * hh + 4 * self.K2 * hh + self.K1 * self.K2) / (hh * hh + self.K1 * hh + self.K1 * self.K2)**2 + self.KW / hh + hh)) 
 
            if (pH_new != 0.0):
                ea = abs(pH_new - pH) / pH_new
                pH = pH_new
                
            if pH > 14.0 or pH < 0.0 or i >= self.alkalinity_constant['imax'] :
                raise ValueError ("The model is unstable or not converging on a solution for pH. Try using the Bisection method to solve pH or reduce the time step")

        return pH

    def Bisection (self):
        pH = self.global_vars['pH']
        xlow = 1.0
        xup  = 13.0 
        flow = self.Function(xlow)
        fup  = self.Function(xup)
        if flow*fup > 0.0 :
            raise ValueError("Bad pH guesses & Try a smaller time step")

        for i in range(0,self.alkalinity_constant['imax']) :
            xr = (xlow + xup) / 2.0
            fr = self.Function(xr)

            if (flow * fr < 0.0):
                xup = xr
            elif (flow * fr > 0.0):
                xlow = xr
            else:
                pH = xr   # find the exact pH value
                break

            if (abs(xup - xlow) / xr < self.alkalinity_constant['es']):
                pH = xr
                break

        return pH