"""
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
import math
from collections import OrderedDict

from src import shared_functions as wqf
from typing import Union, Optional
import sys
import os
import numpy as np

src_path = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(src_path))

from clearwater_modules_python import water_quality_functions as wqf

"""
Globals
    GC
    TwaterC
    depth

GCM_Constants
    k_rc20 = 1
    k_theta = 1.047

    rgc_rc20 = 1
    rgc_theta = 1.047

"""


class GeneralConstituent:
    def __init__(self, global_vars, gsm_constant_changes):
        self.global_vars = global_vars
        self.gsm_constant_changes = gsm_constant_changes

        self.gsm_constant = OrderedDict()
        self.gsm_constant = {
            "order" : 1,
            "k_rc20": 1,
            "k_theta" : 1.047,

            "rgc_rc20" : 1,
            "rgc_theta" : 1.047,    

            "release" : True,
            "settling" : True,
            "settling_rate" : 0.002,
        }

        for key in self.gsm_constant_changes.keys() :
            if key in self.gsm_constant:
                self.gsm_constant[key] = self.gsm_constant_changes[key]
    
    def Calculations(self) :
        print ("Calculating gsm")

        """
        Compute a single general constituent
        Global
            GC: General constituent concentration of a single cell [mg/L]
            TwaterC: Water temperature in degrees Celsius [C]
            depth : Depth of the bed (m)
 
        Constants
            order: Order of reaction kinetics (0, 1, or 2) [unitless]
            k_rc20: Reaction (decay) rate at 20 degrees Celsius. If order=0 [mg/L/d], order=1 [1/d], order=2 [d/mg/L]
            k_theta: Arrhenius temperature correction factor for decay rate [unitless]
            rgc_rc20: Sediment release rate at 20 degrees Celsius [mg*m/L*d]
            rgc_theta: Arrhenius temperature correction factor for sediment release [unitless]
            release: Compute resuspension, True = on; False = off [unitless]
            settling: Compute setting (i.e., bed loss, on/off) [unitless]
            settling_rate: Settling rate (m/s) 
            

        Returns:
            dGCdt (float or np.array): Rate of change of general constituent concentration
        """

        # Temperature corrections
        # k0_rc20 = k1_rc20 = k2_rc20 = 1.0
        # theta_rc20 = theta_rc20 = theta_rc20 = 1.047

        # Compute concentration changes
        gc_decay: float = 0.0               # Rate of decay
        gc_from_bed: float = 0.0            # Sediment release rate (gain)
        gc_settling: float = 0.0            # Settling rate (loss)

        # Correct reaction rate for current temperature
        k_corr: float  = wqf.ArrheniusCorrection(self.global_vars['TwaterC'], self.gsm_constant['k_rc20'], self.gsm_constant['k_theta'])
        gc_decay: float = k_corr * self.global_vars['GC']**self.gsm_constant['order'] #zero-order (mg/L/d), first order (1/d)
    
        if self.gsm_constant['release']:
            rgc_corr: float  = wqf.ArrheniusCorrection(self.global_vars['TwaterC'], self.gsm_constant['rgc_rc20'], self.gsm_constant['rgc_theta'])
            gc_from_bed = rgc_corr / self.global_vars['depth']
        if self.gsm_constant['settling']:
            gc_settling = self.gsm_constant['settling_rate'] * self.global_vars['GC'] / self.global_vars['depth']

        # Compute net rate of change of general constituent concentration
        dGCdt = - gc_decay + gc_from_bed - gc_settling

        return dGCdt
    
    '''
    def float_test():
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
        
        global_vars = OrderedDict()
        global_vars = {
            "GC": 10,
            "TwaterC": 25,
            "depth": 1.0,
        }

        gsm_constant_changes = OrderedDict()
        gsm_constant_changes = {
            "k0_rc20": 1,
            "k1_rc20" : 0.5,
            "k2_rc20" : 1,
            "k0_theta" : 1.047,
            "k1_theta" : 1.047,
            "k2_theta" : 1.047,

            "rgc_rc20" : 1,
            "rgc_theta" : 1.047,
            "order" : 1,
            "release" : True,
            "settling" : True,
            "settling_rate" : 0.002,
        }

        # Compute change of concentration
        dGCdt = GeneralConstituent(global_vars, gsm_constant_changes).Calculations()
        print (dGCdt)
        print(GC)

        GC_new = GC + dGCdt
        print("============================================")
        print("Float Test:")
        print("-----------")
        print("Initial concentration: %.2f" % GC)
        print("Change rate:           %.2f" % dGCdt)
        print("Final concentration:   %.2f" % GC_new)
        print("============================================")
    '''
    '''
    def array_test():
        GC = np.array(10, 5)                        # Initial concentration
        TwaterC = np.array(25.0, 24.0)              # Water temperature
        order = 1                                   # Compute 1st order kinetics
        k_rc20 = np.array(0.5, 0.25)                # Reaction rate at 20 degrees Celsius, decay
        k_theta = np.array(1.047, 1.0)              # Arrhenius temperature correction factor, decay
        rgc_rc20 = np.array(0.5, 0.4)               # Reaction rate at 20 degrees Celsius, settling
        rgc_theta = np.array(1.047, 1.048)          # Arrhenius temperature correction factor, settling
        release = True                              # Turn suspension on
        settling = True                             # Turn settling on
        depth = np.array(1.0, 2.0)                  # Bed depth
        settling_rate = np.array(0.0002, 0.00025)   # Settling rate

        # Compute change of concentration
        dGCdt = GeneralConstituentKinetics(GC, TwaterC, order, k_rc20 = k_rc20, k_theta = k_theta, 
            rgc_rc20 = rgc_rc20, rgc_theta = rgc_theta, release = release, settling = settling, 
            depth = depth, settling_rate = settling_rate)

        GC_new = GC + dGCdt

        print("============================================")
        print("Array Test:")
        print("-----------")
        print("Initial concentration: %.2f" % GC)
        print("Change rate:           %.2f" % dGCdt)
        print("Final concentration:   %.2f" % GC_new)
        print("============================================")


    def mixed_array_float_test():
        GC = np.array(10, 5)    # Initial concentration
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

        print("============================================")
        print("Mixed Array Float Test:")
        print("-----------")
        print("Initial concentration: %.2f" % GC)
        print("Change rate:           %.2f" % dGCdt)
        print("Final concentration:   %.2f" % GC_new)
        print("============================================")



    '''

if __name__ == '__main__':
    GeneralConstituent.float_test()

