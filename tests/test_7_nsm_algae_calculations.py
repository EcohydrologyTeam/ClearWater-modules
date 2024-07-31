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
    return {

        'Ap': 40,
        'Ab': 24,
        'NH4': 0.05,
        'NO3': 5,
        'OrgN': 1.726,
        'N2': 1, 
        'TIP': 0.07,
        'OrgP': 0.24,
        'POC': 4,
        'DOC': 1,
        'DIC': 1,
        'POM': 10, 
        'CBOD': 5,
        'DOX': 8,
        'PX': 1,
        'Alk': 1  

    }

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
        vsap= 0.15,
        growth_rate_option = 1,
        light_limitation_option = 1,
        mu_max_theta= 1.047,
        kdp_theta= 1.047,
        krp_theta= 1.047,
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
        r_alkan= 18.0 / 106.0 / 12.0 / 1000.0,
        r_alkn = 2.0 / 14.0 / 1000.0,
        r_alkden = 4.0 / 14.0 / 1000.0,
        r_alkba = 14.0 / 106.0 / 12.0 / 1000.0,
        r_alkbn =18.0 / 106.0 / 12.0 / 1000.0
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
        BWc= 40,
        BWn=7.2,
        BWp= 1,
        BWa= 5000,
        KLb= 10,
        KsNb= 0.25,
        KsPb=0.125,
        Ksb=10,
        mub_max_20=0.4,
        krb_20=0.2,
        kdb_20=0.3,
        b_growth_rate_option=1,
        b_light_limitation_option=1,
        Fw=0.9,
        Fb=0.9,
        mub_max_theta = 1.047,
        krb_theta = 1.06,
        kdb_theta = 1.047,
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
        KsOxdn=0.1,
        PN=0.5,
        PNb=0.5,
        knit_theta= 1.083,
        kon_theta= 1.074,
        kdnit_theta= 1.08,
        rnh4_theta= 1.047,
        vno3_theta= 1.045,
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
        f_pocb=0.9,
        kpoc_20= 0.005,
        KsOxmc=1.0,
        pCO2 = 383.0,
        FCO2 = 0.2,
        roc = 32.0/12.0,
        kpoc_theta = 1.047,
        kdoc_theta = 1.047,
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
        KsSOD =1,
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
        h2 = 0.1,
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
        apx=1,
        vx=1,
        kdx_theta = 1.07,
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
        rpo4_20 =0,
        kdpo4 = 0.0,
        kop_theta = 1.047,
        rpo4_theta = 1.047,
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
        use_POC = True,
        use_DOC = True,
        use_DOX= True,
        use_DIC= True,
        use_Algae= True,
        use_Balgae= True,
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
        vs = 1,
        SOD_20 = 999,
        SOD_theta = 999,
        theta=1.047,
        vb = 0.01,
        fcom = 0.4,
        kaw_20_user = 0,
        kah_20_user = 1,
        hydraulic_reaeration_option = 1,
        wind_reaeration_option = 1,  
        timestep = 1,    #TODO Dynamic or static?
        depth = 1.5,         #TODO Dynamic or static?
        TwaterC = 25,
        kaw_theta = 1.024,
        kah_theta = 1.024,
        velocity = 1,
        flow = 150,
        topwidth = 100,
        slope = 0.0002,
        shear_velocity = 0.05334,
        pressure_atm = 1,
        wind_speed = 3,
        q_solar = 500,
        Solid = 1,
        lambda0 = .02,
        lambda1 = .0088,
        lambda2 = .054,
        lambdas = .052,
        lambdam = .174, 
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
    return 0.000001

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
    Ap = nsm1.dataset.isel(nsm1_time_step=-1).Ap.values.item()

    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 52.668069

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
    initial_state_dict['Ap'] = 60.0

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 74.849998

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
    initial_state_dict['NH4'] = 0.3

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 52.680781

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
    initial_state_dict = initial_nsm1_state
    default_algae_params['KL'] = 15

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 49.229049

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
    initial_state_dict = initial_nsm1_state
    default_algae_params['KsN'] = 0.02

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 52.803304

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
    initial_state_dict = initial_nsm1_state
    default_algae_params['KsP'] = 0.005

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 50.931103

def test_changed_mu_max(
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
    default_algae_params['mu_max_20'] = 1.5

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 69.809173

def test_changed_kdp(
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
    default_algae_params['kdp_20'] = 0.09

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 55.687635

def test_changed_krp(
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
    default_algae_params['krp_20'] = 0.09

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 58.203941

def test_changed_vsap(
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
    default_algae_params['vsap'] = 0.2

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 51.334735

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
    initial_state_dict['NO3'] = 10

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 52.802629    

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
    initial_state_dict['TIP'] = 0.1

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 52.842286  

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
    initial_state_dict = initial_nsm1_state
    default_gvars_params['TwaterC'] = 35

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 62.384696  

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
    initial_state_dict = initial_nsm1_state
    default_gvars_params['depth'] = 2

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 47.171170

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
    initial_state_dict = initial_nsm1_state
    default_gvars_params['q_solar'] = 250

def test_changed_Fr_PAR(
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
    default_gvars_params['Fr_PAR'] = 0.2

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 45.185221
  
def test_changed_use_NH4(
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
    default_gp_params['use_NH4'] = False

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 52.665375     

def test_changed_use_NO3(
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
    default_gp_params['use_NO3'] = False

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 37.582388         

def test_changed_use_TIP(
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
    default_gp_params['use_TIP'] = False

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 53.255764    

def test_changed_g2_l1(
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
    default_algae_params['growth_rate_option'] = 2

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 52.939611    

def test_changed_g3_l1(
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
    default_algae_params['growth_rate_option'] = 3

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 53.096967    

def test_changed_g1_l2(
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
    #default_algae_params['growth_rate_option'] = 3
    default_algae_params['light_limitation_option'] = 2

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 59.847275    

def test_changed_g2_l2(
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
    default_algae_params['growth_rate_option'] = 2
    default_algae_params['light_limitation_option'] = 2

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 60.175682   

def test_changed_g3_l2(
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
    default_algae_params['growth_rate_option'] = 3
    default_algae_params['light_limitation_option'] = 2

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 60.365992  

def test_changed_g1_l3(
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
    default_algae_params['growth_rate_option'] = 1
    default_algae_params['light_limitation_option'] = 3

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 40.479472

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
    initial_state_dict = initial_nsm1_state
    default_gvars_params['lambda0'] = 5

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 31.938218

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
    initial_state_dict = initial_nsm1_state
    default_gvars_params['lambda1'] = 0.03

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 46.601473

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
    initial_state_dict = initial_nsm1_state
    default_gvars_params['lambda2'] = 0.02

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 55.635473

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
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    default_gvars_params['lambdam'] = 0.009

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 62.784764

def test_changed_use_POC(
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
    default_gp_params['use_POC'] = False

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 63.113783

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
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    initial_state_dict['POC'] = 6

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 46.458270

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
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    default_gvars_params['fcom'] = 0.6

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
    Ap = nsm1.dataset.isel(
        nsm1_time_step=-1).Ap.values.item()
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 56.937679