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
        use_POC = True,
        use_DOC = True,
        use_DIC= True,
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
        pressure_atm = 1013.25,
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

def test_defaults_POC(
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
    POC = nsm1.dataset.isel(nsm1_time_step=-1).POC.values.item()

    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    
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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_defaults_DIC(
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
    DIC = nsm1.dataset.isel(nsm1_time_step=-1).DIC.values.item()

    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_POC(
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
    """Test the model with default parameters and changed POC."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['POC'] = 2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 3.98
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.22
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_DOC(
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
    """Test the model with default parameters and changed DOC."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['DOC'] = 2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 2.23
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_DIC(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['DIC'] = 2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 1.54
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['DOX'] = 5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_CBOD(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['CBOD'] = 10000

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.81   
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['Ap'] = 3

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.08
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.21
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['Ab'] = 2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 4.70
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.06
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['NH4'] = 7000

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['NO3'] = 5000

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77

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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['TIP'] = 5000

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['TwaterC'] = 15

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 5.58
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.15
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.85
    
def test_changed_fpocp(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_carbon_dict = default_carbon_params
    default_carbon_dict['f_pocp'] = 0.4

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.17
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.38
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_fpocb(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_carbon_dict = default_carbon_params
    default_carbon_dict['f_pocb'] = 0.4

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 5.33
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 2.22
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_KsOxmc(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_carbon_dict = default_carbon_params
    default_carbon_dict['KsOxmc'] = 3

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77

def test_changed_pCO2(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_carbon_dict = default_carbon_params
    default_carbon_dict['pCO2'] = 425

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_FCO2(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_carbon_dict = default_carbon_params
    default_carbon_dict['FCO2'] = 0.5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.42
    
def test_changed_vsoc(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['vsoc'] = 0.1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.04
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['depth'] = 5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 5.10
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.10
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['light_limitation_option'] = 2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_fcom(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['fcom'] = 0.01

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['Fb'] = 0.01

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 4.56
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.05
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_Fw(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['Fw'] = 0.08

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 4.70
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.06
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_kpoc_20(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_carbon_dict = default_carbon_params
    default_carbon_dict['kpoc_20'] = 0.01

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.27
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.27
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_kdoc_20(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_carbon_dict = default_carbon_params
    default_carbon_dict['kdoc_20'] = 0.1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.14
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_SOD_20(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['SOD_20'] = 1000

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.79
    
def test_changed_kdp_20(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['kdp_20'] = 1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 7.72
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.40
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_kdb_20(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['kdb_20'] = 1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 10.41
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.70
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_kbod_20(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_CBOD_dict = default_CBOD_params
    default_CBOD_dict['kbod_20'] = 1000

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == -45.77
    
def test_changed_kah_20_user(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['kah_20_user'] = 5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == -0.16
    
def test_changed_kaw_20_user(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['kaw_20_user'] = 5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == -0.0063
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['krp_20'] = 10000

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 2.31
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['krb_20'] = 10000

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 6.81
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['AWa'] = 5000

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.10
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.22
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['AWc'] = 100

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.68
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.28
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['AWd'] = 500

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['AWp'] = 5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['AWn'] = 1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['BWd'] = 1000

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 4.72
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.06
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['BWc'] = 400

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 22.15
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 3.00
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['BWn'] = 1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['BWp'] = 10

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['BWa'] = 100

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['q_solar'] = 250

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['mu_max_20'] = 10

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['mub_max_20'] = 10

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['lambda0'] = 0.5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['lambda1'] = 0.1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['lambda2'] = 1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_lambdam(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['lambdam'] = 5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['Ksb'] = 1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_KsOxbod(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_CBOD_dict = default_CBOD_params
    default_CBOD_dict['KsOxbod'] = 10

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['KLb'] = 1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['KsNb'] = 5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['KsPb'] = 5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['KL'] = 0.1

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['KsN'] = 10

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['KsP'] = 10

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_KsSOD(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_DOX_dict = default_DOX_params
    default_DOX_dict['KsSOD'] = 5

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.30
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_kpoc_theta(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_carbon_dict = default_carbon_params
    default_carbon_dict['kpoc_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.28
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.27
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_kdoc_theta(
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
    """Test the model with default parameters and changed DIC."""
    # alter parameters as necessary
    default_carbon_dict = default_carbon_params
    default_carbon_dict['kdoc_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.31
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.23
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_SOD_theta(
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
    """Test the model with default parameters and changed SOD_theta."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['SOD_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.31
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_kdp_theta(
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
    """Test the model with default parameters and changed kdp_theta."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['kdp_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.55
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.27
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_kdb_theta(
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
    """Test the model with default parameters and changed kdb_theta."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['kdb_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 8.03
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.43
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_kbod_theta(
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
    """Test the model with default parameters and changed kbod_theta."""
    # alter parameters as necessary
    default_CBOD_dict = default_CBOD_params
    default_CBOD_dict['kbod_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.31
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
def test_changed_kah_theta(
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
    """Test the model with default parameters and changed kah_theta."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['kah_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.31
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.54
    
def test_changed_kaw_theta(
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
    """Test the model with default parameters and changed kaw_theta."""
    # alter parameters as necessary
    default_gvars_dict = default_gvars_params
    default_gvars_dict['kaw_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.31
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed krp_theta."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['krp_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.31
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77
    
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
    """Test the model with default parameters and changed krb_theta."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['krb_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.31
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77

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
    """Test the model with default parameters and changed mu_max_theta."""
    # alter parameters as necessary
    default_algae_dict = default_algae_params
    default_algae_dict['mu_max_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.31
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77

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
    """Test the model with default parameters and changed mub_max_theta."""
    # alter parameters as necessary
    default_balgae_dict = default_balgae_params
    default_balgae_dict['mub_max_theta'] = 1.2

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
    POC = nsm1.dataset.isel(
        nsm1_time_step=-1).POC.values.item()
    assert isinstance(POC, float)
    assert pytest.approx(POC, tolerance) == 6.31
    DOC = nsm1.dataset.isel(
        nsm1_time_step=-1).DOC.values.item()
    assert isinstance(DOC, float)
    assert pytest.approx(DOC, tolerance) == 1.24
    DIC = nsm1.dataset.isel(
        nsm1_time_step=-1).DIC.values.item()
    assert isinstance(DIC, float)
    assert pytest.approx(DIC, tolerance) == 0.77