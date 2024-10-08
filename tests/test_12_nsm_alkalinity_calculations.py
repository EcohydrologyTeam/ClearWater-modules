from numba import (
    types,
    typed,
)
import pytest

from clearwater_modules.nsm1 import NutrientBudget
from clearwater_modules.nsm1.constants import (
    AlgaeStaticVariables,
    AlkalinityStaticVariables,
    BalgaeStaticVariables,
    NitrogenStaticVariables,
    CarbonStaticVariables,
    CBODStaticVariables,
    DOXStaticVariables,
    N2StaticVariables,
    POMStaticVariables,
    PathogenStaticVariables,
    PhosphorusStaticVariables,
    GlobalParameters,
    GlobalVars
)


@pytest.fixture(scope='function')
def initial_nsm1_state() -> dict[str, float]:
    """Return initial state values for the model."""
    return {            'Ap': 36.77,
                        'Ab': 24,
                        'NH4': .063,
                        'NO3': 5.54,
                        'OrgN': 1.726,
                        'N2': 1,
                        'TIP': 0.071,
                        'OrgP': 0.24,
                        'POC': 4.356,
                        'DOC': 1,
                        'DIC': 1,
                        'POM': 1,
                        'CBOD': 1,
                        'DOX': 8,
                        'PX': 1,
                        'Alk': 1}

@pytest.fixture(scope='module')
def time_steps() -> int:
    return 1

@pytest.fixture(scope='function')
def default_algae_params() -> AlgaeStaticVariables:
    """Returns default algae static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return AlgaeStaticVariables(
        AWd = 100,
        AWc= 40,
        AWn= 7.2,
        AWp= 1,
        AWa= 1000,
        KL= 10,
        KsN= 0.04,
        KsP= 0.0012,
        mu_max_20= 1,
        kdp_20= 0.15,
        krp_20= 0.2,
        mu_max_theta= 1.047,
        kdp_theta= 1.047,
        krp_theta= 1.047,
        vsap= 0.15,
        growth_rate_option = 1,
        light_limitation_option = 1,
    )

@pytest.fixture(scope='function')
def default_alkalinity_params() -> AlkalinityStaticVariables:
    """Returns default alkalinity static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return AlkalinityStaticVariables(
        r_alkaa = 14.0 / 106.0 / 12.0 / 1000.0,
        r_alkan = 18.0 / 106.0 / 12.0 / 1000.0,
        r_alkn = 2.0 / 14.0 / 1000.0,
        r_alkden = 4.0 / 14.0 / 1000.0,
        r_alkba = 14.0 / 106.0 / 12.0 / 1000.0,
        r_alkbn = 18.0 / 106.0 / 12.0 / 1000.0
    )

@pytest.fixture(scope='function')
def default_balgae_params() -> BalgaeStaticVariables:
    """Returns default Benthic Algae static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return BalgaeStaticVariables(
        BWd = 100,
        BWc = 40,
        BWn = 7.2,
        BWp = 1,
        BWa = 5000,
        KLb = 10,
        KsNb= 0.25,
        KsPb=0.125,
        Ksb=10,
        mub_max_20=0.4,
        krb_20=0.2,
        kdb_20=0.3,
        mub_max_theta = 1.047,
        krb_theta = 1.047,
        kdb_theta = 1.047,
        b_growth_rate_option=1,
        b_light_limitation_option=1,
        Fw=0.9,
        Fb=0.9
    )

@pytest.fixture(scope='function')
def default_nitrogen_params() -> NitrogenStaticVariables:
    """Returns default nitrogen static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return NitrogenStaticVariables(
        KNR= 0.6 ,
        knit_20= 0.1,
        kon_20=0.1,
        kdnit_20=0.002,
        rnh4_20=0,
        vno3_20=0,
        knit_theta= 1.047,
        kon_theta= 1.047,
        kdnit_theta= 1.047,
        rnh4_theta= 1.047,
        vno3_theta= 1.047,
        KsOxdn=0.1,
        PN=0.5,
        PNb=0.5
    )

@pytest.fixture(scope='function')
def default_carbon_params() -> CarbonStaticVariables:
    """Returns default carbon static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return CarbonStaticVariables(
        f_pocp = 0.9,
        kdoc_20= 0.01,
        kdoc_theta = 1.047,
        f_pocb=0.9,
        kpoc_20= 0.005,
        kpoc_theta = 1.047,
        KsOxmc=1.0,
        pCO2 = 383.0,
        FCO2 = 0.2,
        roc = 32.0/12.0
    )

@pytest.fixture(scope='function')
def default_CBOD_params() -> CBODStaticVariables:
    """Returns default CBOD static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return CBODStaticVariables(
        KsOxbod = 0.5,
        kbod_20 =  0.12,
        ksbod_20 = 0.0,
        kbod_theta =  1.047,
        ksbod_theta = 1.047
    )

@pytest.fixture(scope='function')
def default_DOX_params() -> DOXStaticVariables:
    """Returns default DOX static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return DOXStaticVariables(
        ron = 2.0 * 32.0 / 14.0,
        KsSOD = 1,
    )

@pytest.fixture(scope='function')
def default_N2_params() -> N2StaticVariables:
    """Returns default N2 static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return N2StaticVariables(

    )

@pytest.fixture(scope='function')
def default_POM_params() -> POMStaticVariables:
    """Returns default POM static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return POMStaticVariables(
        kpom_20 = 0.1,
        kpom_theta = 1.047
    )

@pytest.fixture(scope='function')
def default_pathogen_params() -> PathogenStaticVariables:
    """Returns default Pathogens static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return PathogenStaticVariables(
        kdx_20=0.8,
        kdx_theta = 1.047,
        apx=1,
        vx=1
    )

@pytest.fixture(scope='function')
def default_phosphorus_params() -> PhosphorusStaticVariables:
    """Returns default phosphorus static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return PhosphorusStaticVariables(
        kop_20 = 0.1,
        rpo4_20 = 0,
        kop_theta = 1.047,
        rpo4_theta = 1.047,
        kdpo4 = 100.0
    )

@pytest.fixture(scope='function')
def default_gp_params() -> GlobalParameters:
    """Returns default global parameter static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return GlobalParameters(
        use_NH4= True,
        use_NO3= True, 
        use_OrgN= True,
        use_OrgP = True,
        use_TIP= True,  
        use_SedFlux= False,
        use_DOX= True,
        use_Algae= True,
        use_Balgae= True,
        use_POC = False,
        use_DOC = False,
        use_DIC= False,
        use_N2 = True,
        use_Pathogen = True,
        use_Alk = True,
        use_POM = True 
    )

@pytest.fixture(scope='function')
def default_gvars_params() -> GlobalVars:
    """Returns default global variables static variable values for the model.

    NOTE: As of now (3/18/2022) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return GlobalVars(
        vson = 0.01,
        vsoc = 0.01,
        vsop = 0.01,
        vs = 0.01,
        SOD_20 = 0.5,
        SOD_theta = 1.047,
        vb = 0.01,
        fcom = 0.4,
        kaw_20_user = 0,
        kah_20_user = 1,
        kaw_theta = 1.047,
        kah_theta = 1.047,
        hydraulic_reaeration_option = 1,
        wind_reaeration_option = 1,  
        dt = 1,    #TODO Dynamic or static?
        depth = 1.5,         #TODO Dynamic or static?
        TwaterC = 25,
        theta = 1.047,
        velocity = 1,
        flow = 150,
        topwidth = 100,
        slope = 0.0002,
        shear_velocity = 0.05334,
        pressure_mb = 1013.25,
        wind_speed = 3,
        q_solar = 500,
        Solid = 1,
        lambda0 = 0.02,
        lambda1 = 0.0088,
        lambda2 = 0.054,
        lambdas = 0.056,
        lambdam = 0.174, 
        Fr_PAR = .47  
    )

def get_nutrient_budget_instance(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,

) -> NutrientBudget:
    """Return an instance of the NSM1 class."""
    return NutrientBudget(
        time_steps=time_steps,
        initial_state_values=initial_nsm1_state,
        algae_parameters=default_algae_params,
        alkalinity_parameters=default_alkalinity_params,
        balgae_parameters=default_balgae_params,
        nitrogen_parameters=default_nitrogen_params,
        carbon_parameters=default_carbon_params,
        CBOD_parameters=default_CBOD_params,
        DOX_parameters=default_DOX_params,
        N2_parameters=default_N2_params,
        POM_parameters=default_POM_params,
        pathogen_parameters=default_pathogen_params,
        phosphorus_parameters=default_phosphorus_params,
        global_parameters=default_gp_params,
        global_vars=default_gvars_params,
        time_dim='nsm1_time_step',
    )

@pytest.fixture(scope='module')
def tolerance() -> float:
    """Controls the precision of the pytest.approx() function."""
    return 0.01

def test_defaults(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary

    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.27
    
def test_changed_Ap(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['Ap'] = 10
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.31
    
def test_changed_Ab(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['Ab'] = 10
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.75
    
def test_changed_Alk(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['Alk'] = 5
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 7.27
    
def test_changed_NH4(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['NH4'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.10
    
def test_changed_TIP(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['TIP'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.54
    
def test_changed_DOX(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['DOX'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.32
    
def test_changed_NO3(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['NO3'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.10
    
def test_changed_TwaterC(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_gvars_dict = default_gvars_params
    initial_gvars_dict['TwaterC'] = 15
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.43
    
def test_changed_q_solar(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_gvars_dict = default_gvars_params
    initial_gvars_dict['q_solar'] = 250
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.16
    
def test_changed_KL(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['KL'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.36
    
def test_changed_KLb(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['KLb'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.30
    
def test_changed_lambda0(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_gvars_dict = default_gvars_params
    initial_gvars_dict['lambda0'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.09
    
def test_changed_lambda1(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_gvars_dict = default_gvars_params
    initial_gvars_dict['lambda1'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.02
    
def test_changed_lambda2(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_gvars_dict = default_gvars_params
    initial_gvars_dict['lambda2'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.18
    
def test_changed_KsN(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['KsN'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.10
    
def test_changed_KsP(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['KsP'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.20
    
def test_changed_KsNb(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['KsNb'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.25
    
def test_changed_KsPb(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['KsPb'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.13
    
def test_changed_light_limitation_option(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['light_limitation_option'] = 2
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.37
    
def test_changed_b_light_limitation_option(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['b_light_limitation_option'] = 2
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.30
    
def test_changed_Fb(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['Fb'] = 0.5
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.84
    
def test_changed_Ksb(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['Ksb'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.12
    
def test_changed_AWa(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['AWa'] = 100
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 15.43
    
def test_changed_AWc(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['AWc'] = 20
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.59
    
def test_changed_AWd(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['AWd'] = 50
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.27
    
def test_changed_AWp(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['AWp'] = 0.1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.27
    
def test_changed_AWn(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['AWn'] = 0.1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.27
    
def test_changed_BWd(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['BWd'] = 50
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 4.24
    
def test_changed_BWc(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['BWc'] = 10
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.54
    
def test_changed_BWn(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['BWn'] = 0.1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.27
    
def test_changed_BWp(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['BWp'] = 0.1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.27
    
def test_changed_BWa(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['BWa'] = 1000
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.27
    
def test_changed_depth(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_gvars_dict = default_gvars_params
    initial_gvars_dict['depth'] = 5
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.20
    
def test_changed_mu_max_20(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['mu_max_20'] = 0.2
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.35
    
def test_changed_mub_max_20(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['mub_max_20'] = 1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.53
    
def test_changed_PN(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_nitrogen_dict = default_nitrogen_params
    initial_nitrogen_dict['PN'] = 0.9
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.10
    
def test_changed_PNb(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_nitrogen_dict = default_nitrogen_params
    initial_nitrogen_dict['PNb'] = 0.9
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.24
    
def test_changed_PNb(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['PNb'] = 0.9
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.24
    
def test_changed_knit_20(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_nitrogen_dict = default_nitrogen_params
    initial_nitrogen_dict['knit_20'] = 0.01
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.31
    
def test_changed_KsOxdn(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_nitrogen_dict = default_nitrogen_params
    initial_nitrogen_dict['KsOxdn'] = 0.01
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.27
    
def test_changed_kdnit_20(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_nitrogen_dict = default_nitrogen_params
    initial_nitrogen_dict['kdnit_20'] = 0.1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.40
    
def test_changed_kdnit_20(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_nitrogen_dict = default_nitrogen_params
    initial_nitrogen_dict['kdnit_20'] = 0.1
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.40
    
def test_changed_krp_20(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['krp_20'] = 0.01
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.08
    
def test_changed_krp_20(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['krp_20'] = 0.01
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.08
    
def test_changed_krb_20(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['krb_20'] = 0.01
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 2.51
    
def test_changed_knit_theta(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_nitrogen_dict = default_nitrogen_params
    initial_nitrogen_dict['knit_theta'] = 1.2
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.21
    
def test_changed_kdnit_theta(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_nitrogen_dict = default_nitrogen_params
    initial_nitrogen_dict['kdnit_theta'] = 1.2
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.27
    
def test_changed_mu_max_theta(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['mu_max_theta'] = 1.2
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 4.39
    
def test_changed_mub_max_theta(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['mub_max_theta'] = 1.2
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.44
    
def test_changed_krp_theta(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_algae_dict = default_algae_params
    initial_algae_dict['krp_theta'] = 1.2
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 3.47
    
def test_changed_krb_theta(
    time_steps,
    initial_nsm1_state,
    default_algae_params,
    default_alkalinity_params,
    default_balgae_params,
    default_nitrogen_params,
    default_carbon_params,
    default_CBOD_params,
    default_DOX_params,
    default_N2_params,
    default_POM_params,
    default_pathogen_params,
    default_phosphorus_params,
    default_gp_params,
    default_gvars_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_balgae_dict = default_balgae_params
    initial_balgae_dict['krb_theta'] = 1.2
    
    # instantiate the model
    nsm1: NutrientBudget = get_nutrient_budget_instance(
        time_steps=time_steps,
        initial_nsm1_state=initial_nsm1_state,
        default_algae_params=default_algae_params,
        default_alkalinity_params=default_alkalinity_params,
        default_balgae_params=default_balgae_params,
        default_nitrogen_params=default_nitrogen_params,
        default_carbon_params=default_carbon_params,
        default_CBOD_params=default_CBOD_params,
        default_DOX_params=default_DOX_params,
        default_N2_params=default_N2_params,
        default_POM_params=default_POM_params,
        default_pathogen_params=default_pathogen_params,
        default_phosphorus_params=default_phosphorus_params,
        default_gp_params=default_gp_params,
        default_gvars_params=default_gvars_params
    )

    # Run the model
    nsm1.increment_timestep()
    Alk = nsm1.dataset.isel(nsm1_time_step=-1).Alk.values.item()

    assert isinstance(Alk, float)
    assert pytest.approx(Alk, tolerance) == 4.05