from collections import Orderedtyped.Dict

import pytest
from numba import (
    types,
    typed,
)
import os
import sys

# TODO: there is BenthicAlgae in clearwater_modules.nsm1 and clearwater_modules.nsm2
# I am not sure which one is the correct one to use
from clearwater_modules.nsm1 import BenthicAlgae

@pyest.fixture(scope='module')
def tolerance() -> float:
    return 0.001

@pytest.fixture(scope='module')
def module_choices() -> typed.Dict:
    module_choices = typed.Dict.empty(
        key_type=types.unicode_type, 
        value_type=types.float64,
    )
        
    module_choices = {
        'use_Algae': True,
        'use_NH4': True,
        'use_NO3': True,
        'use_TIP': True,
        'use_POC': True,
        'use_DOC': True,

        'use_BAlgae': True,
        'use_OrgN' : True,
        'use_OrgP' : True,

        'use_SedFlux' : True,
        'use_DOX': True,

        'use_DIC': False,
        'use_N2' : False,
        'use_Pathogen' : False,
        'use_Alk' : False,
        'use_POM2' : False,
    }
    return module_choices

@pytest.fixture(scope='module')
def vars() -> typed.Dict:
    module_choices = typed.Dict.empty(
        key_type=types.unicode_type, 
        value_type=types.float64,
    )
 
    module_choices = {
        'Ap': 100,
        'NH4':100,
        'NO3': 100,
        'TIP':100,
        'TwaterC':20,
        'depth':1,
        'Ab':100,
        'DOX' : 100,
        'OrgN' : 200,
        'vson': 0.01,
        'lambda': 1,
        'fdp': 0.5,
        'PAR': 100
    }
    return module_choices

@pytest.fixture(scope='module')
def constant_changes() -> typed.Dict:
    constant_changes = typed.Dict.empty(
        key_type=types.unicode_type,
        value_type=types.float64,
    )

    constant_changes = {
        'BWd': 100,       
        'BWc': 40,      
        'Bwn' : 7.2,      
        'BWp' : 1,      
        'BWa' : 3500,       
        'KLb': 10,       
        'KsNb' : 0.25,       
        'KsPb' : 0.125,    
        'Ksb' : 10,       
        'mub_max' : 0.4,   
        'krb' : 0.2,     
        'kdb': 0.3,   
        'b_growth_rate_option' : 1,
        'b_light_limitation_option' : 1,
        'Fw' : 0.9,
        'Fb' : 0.9, 
    }
    
    # TODO: Why were these parameter commented out earlier?
    # For now I am remaping to be blank
    constant_changes = {}
    return {}


def test_original (module_choices, vars, constant_changes, tolerance) -> None:
    
    Change, path = BenthicAlgae(
        module_choices,
        vars, 
        constant_changes,
    ).Calculations()
    
    assert pytest.approx(Change, tolerance) == -47.1515


#Change Ab
def test_Ab (module_choices, vars, constant_changes, tolerance) -> None:
    
    vars['Ab'] = 200

    Change, path = BenthicAlgae(
        module_choices,
        vars, 
        constant_changes,
    ).Calculations()
 
    assert pytest.approx(Change, tolerance) == -97.0159


#Change mub_max
def test_mub_max (module_choices, vars, constant_changes, tolerance) -> None:
    constant_changes['mub_max'] = 0.2

    Change, path = BenthicAlgae(
        module_choices,
        vars,
        constant_changes,
    ).Calculations()
   
    assert pytest.approx(Change, tolerance) == -48.57576

#Change growth_rate_option
def test_GRO (module_choices, vars, constant_changes, tolerance) -> None:
    constant_changes['b_growth_rate_option'] = 2

    Change, path = BenthicAlgae(
        module_choices,
        vars, 
        constant_changes,
    ).Calculations()

    assert pytest.approx(Change, tolerance) == -47.14797

#Change light_limitation_option
def test_LLO (module_choices, vars, constant_changes, tolerance) -> None:
    constant_changes['b_light_limitation_option'] = 2

    Change, path = BenthicAlgae(
        module_choices, 
        vars, 
        constant_changes,
    ).Calculations()

    assert pytest.approx(Change, tolerance) == -46.50408888
