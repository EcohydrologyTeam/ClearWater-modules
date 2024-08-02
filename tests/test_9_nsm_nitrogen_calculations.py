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

        'Ap': 36.77,
        'Ab': 24,
        'NH4': 0.063,
        'NO3': 5.54,
        'OrgN': 1.726,
        'N2': 1, 
        'TIP': 0.071,
        'OrgP': 0.25,
        'POC': 4.356,
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
        kon_theta= 1.047,
        kdnit_theta= 1.045,
        rnh4_theta= 1.074,
        vno3_theta= 1.08,
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
        use_TIP= False,  
        use_SedFlux= False,
        use_POC = False,
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

    nsm1.increment_timestep()

 

    krb_tc = nsm1.dataset.isel(nsm1_time_step=-1).krb_tc.item()
    assert isinstance(krb_tc, float)
    print("krb_tc",krb_tc)

    AbRespiration = nsm1.dataset.isel(nsm1_time_step=-1).AbRespiration.item()
    assert isinstance(AbRespiration, float)
    print("AbRespiration",AbRespiration)


    NH4_AbRespiration = nsm1.dataset.isel(nsm1_time_step=-1).NH4_AbRespiration.item()
    assert isinstance(NH4_AbRespiration, float)
    print("NH4_AbRespiration",NH4_AbRespiration)


    # Run the model
    NH4 = nsm1.dataset.isel(nsm1_time_step=-1).NH4.values.item()

    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.6101494

    NO3 = nsm1.dataset.isel(nsm1_time_step=-1).NO3.values.item()

    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.1260636

    OrgN = nsm1.dataset.isel(nsm1_time_step=-1).OrgN.values.item()

    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.899502

def test_changed_knit(
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
    default_nitrogen_params['knit_20'] = 5

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.154017957758354

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.58219498187554

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.89950169149861

def test_changed_kon(
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
    default_nitrogen_params['kon_20'] = 0.5

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 1.47877811233963

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.12606356028487

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.03087295850801

def test_changed_kdnit(
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
    default_nitrogen_params['kdnit_20'] = 2

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.610149

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 4.955769

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.8995017

def test_changed_rnh4(
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
    default_nitrogen_params['rnh4_20'] = -0.5

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.133827915386814

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.12606356028487

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.89950169

def test_changed_rnh4_KsOxdn(
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
    default_nitrogen_params['rnh4_20'] = 1
    default_nitrogen_params['KsOxdn'] = 10

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 1.5627923

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.11856308

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.89950169

def test_changed_vno3(
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
    default_nitrogen_params['vno3_20'] = 0.9

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.610149379

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 0.242017033

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.89950169

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
    default_gvars_params['depth'] = 3

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.472669296

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.23749182

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.729153886

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
    default_gvars_params['TwaterC'] = 10

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.31604458

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.330296324

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.8073881

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
    initial_nsm1_state['NH4'] = 0.01

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.569020898

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.114263537

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.89950169

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
    initial_nsm1_state['NO3'] = 20

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.61359955

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 19.57664053

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.89950169

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
    initial_nsm1_state['Ap'] = 5

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.5553608647

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.3696889

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.85633245

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
    initial_nsm1_state['DOX'] = 1

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.6152233076

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.1199048523

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.899501691

def test_changed_OrgN(
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
    initial_nsm1_state['OrgN'] = 7

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 1.2736991965

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.12606356

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 6.47479187

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.063

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.1120288

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.89950169

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.40213002

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.54

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.89950169

def test_changed_use_OrgN(
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
    default_gp_params['use_OrgN'] = False

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.392992196

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.126063560

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.726

def test_changed_use_DOX(
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
    default_gp_params['use_DOX'] = False

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.61007213

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.12631127

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.89950169

def test_changed_use_Algae(
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
    default_gp_params['use_Algae'] = False

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.54676037

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.409982005

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.849538428

def test_changed_use_Balgae(
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
    default_gp_params['use_Balgae'] = False

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.33406138337

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.24974272

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.547299413

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
    initial_state_dict = initial_nsm1_state
    default_nitrogen_params['PN'] = 0.3

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.6120854746

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.124127465

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.899501691

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
    initial_state_dict = initial_nsm1_state
    default_nitrogen_params['PNb'] = 0.7

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.60832257266

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.127890367

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.89950169

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
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_nsm1_state
    default_balgae_params['Fw'] = 0.5

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.610149379

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.12606356028

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.7429673455

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
    initial_state_dict = initial_nsm1_state
    default_balgae_params['Fb'] = 0.4

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.456767159

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.194774206

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.703833759

def test_changed_use_Algae_use_Balgae(
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
    default_gp_params['use_Algae'] = False
    default_gp_params['use_Balgae'] = False    

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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.2708483787

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.549138339

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.49733615

def test_changed_use_Algae_use_NH4(
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
    default_gp_params['use_Algae'] = False
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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.063

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.39915906

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.849538428

def test_changed_use_Algae_use_NH4_use_NO3(
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
    default_gp_params['use_Algae'] = False
    default_gp_params['use_NH4'] = False    
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
    NH4 = nsm1.dataset.isel(
        nsm1_time_step=-1).NH4.values.item()
    assert isinstance(NH4, float)
    assert pytest.approx(NH4, tolerance) == 0.063

    NO3 = nsm1.dataset.isel(
        nsm1_time_step=-1).NO3.values.item()
    assert isinstance(NO3, float)
    assert pytest.approx(NO3, tolerance) == 5.54

    OrgN = nsm1.dataset.isel(
        nsm1_time_step=-1).OrgN.values.item()
    assert isinstance(OrgN, float)
    assert pytest.approx(OrgN, tolerance) == 1.849538428