from numba import (
    types,
    typed,
)
import pytest

from clearwater_modules_python.tsm import EnergyBudget


@pytest.fixture(scope='module')
def use_sed_temp() -> bool:
    return True


@pytest.fixture(scope='module')
def vars() -> typed.Dict:
    vars = typed.Dict.empty(
        key_type=types.unicode_type,
        value_type=types.float64,
    )
    vars = {
        'TwaterC': 20,
        'surface_area': 1,
        'volume': 1
    }
    return vars


@pytest.fixture(scope='module')
def t_changes() -> typed.Dict:
    return typed.Dict.empty(
        key_type=types.unicode_type,
        value_type=types.float64,
    )


@pytest.fixture(scope='module')
def met_changes() -> typed.Dict:
    return typed.Dict.empty(
        key_type=types.unicode_type,
        value_type=types.float64,
    )


@pytest.fixture(scope='module')
def tolerance() -> float:
    """Controls the precision of the pytest.approx() function."""
    return 0.0001


@pytest.fixture(scope='module')
def answers(function_name: str) -> float:
    answers = {
        'test_TwaterC_20': 19.997,
        'test_TwaterC_40': 39.964,
        'test_surface_area_2': 19.994,
        'test_surface_area_4': 19.987,
        'test_volume_2': 19.998,
        'test_volume_4': 19.999,
        'test_TairC_30': 20,
        'test_TairC_40': 20.001,
        'test_TsedC_10': 19.999,
        'test_TsedC_15': 20.001,
        'test_q_solar_450': 19.997,
        'test_q_solar_350': 19.996,
        'test_wind_kh_kw_0_5': 19.997,
        'test_wind_kh_kw_1_5': 19.997,
        'test_eair_mb_2': 19.997,
        'test_eair_mb_5': 19.997,
        'test_pressure_mb_970': 19.997,
        'test_pressure_mb_1050': 19.997,
        'test_cloudiness_0': 19.997,
        'test_cloudiness_0_5': 19.997,
        'test_wind_speed_5': 19.997,
        'test_wind_speed_30': 19.983,
        'test_wind_a_1En7': 19.997,
        'test_wind_a_7En7': 19.997,
        'test_wind_b_1En6': 19.997,
        'test_wind_b_2En6': 19.996,
        'test_wind_c_0_5': 19.998,
        'test_wind_c_3': 19.980,
    }
    return answers[function_name]


def test_TwaterC_20(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    vars['TwaterC'] = 20

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_TwaterC_20')


def test_TwaterC_40(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    vars['TwaterC'] = 40

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_TwaterC_40')


def test_surface_area_2(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    vars['surface_area'] = 2

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_surface_area_2')


def test_surface_area_4(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    vars['surface_area'] = 4

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_surface_area_4')


def test_volume_2(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    vars['volume'] = 2

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_volume_2')


def test_volume_4(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    vars['volume'] = 4

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_volume_4')


def test_TairC_30(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['TairC'] = 30

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_TairC_30')


def test_TairC_40(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['TairC'] = 40

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_TairC_40')


def test_TsedC_10(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['TsedC'] = 10

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_TsedC_10')


def test_TsedC_15(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['TsedC'] = 15

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_TsedC_15')


def test_q_solar_450(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['q_solar'] = 450

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_q_solar_450')


def test_q_solar_350(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['q_solar'] = 350

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_q_solar_350')


def test_wind_kh_kw_0_5(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['wind_kh_kw'] = .5

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_wind_kh_kw_0_5')


def test_wind_kh_kw_1_5(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['wind_kh_kw'] = 1.5

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_wind_kh_kw_1_5')


def test_eair_mb_2(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['eair_mb'] = 2

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_eair_mb_2')


def test_eair_mb_5(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['eair_mb'] = 5

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_eair_mb_5')


def test_pressure_mb_970(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['pressure_mb'] = 970

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_pressure_mb_970')


def test_pressure_mb_1050(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['pressure_mb'] = 1050

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers(
        'test_pressure_mb_1050')


def test_cloudiness_0(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['cloudiness'] = 0

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_cloudiness_0')


def test_cloudiness_0_5(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['cloudiness'] = .5

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_cloudiness_0_5')


def test_wind_speed_5(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['wind_speed'] = 5

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_wind_speed_5')


def test_wind_speed_30(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['wind_speed'] = 30

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_wind_speed_30')


def test_wind_a_1En7(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['wind_a'] = 1*10**-7

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_wind_a_1En7')


def test_wind_a_7En7(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['wind_a'] = 7*10**-7

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_wind_a_7En7')


def test_wind_b_1En6(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['wind_b'] = 1*10**-6

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_wind_b_1En6')


def test_wind_b_2En6(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['wind_b'] = 2*10**-6

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_wind_b_2En6')


def test_wind_c_0_5(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['wind_c'] = .5

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_wind_c_0_5')


def test_wind_c_3(use_sed_temp, vars, met_changes, t_changes, tolerance, answers) -> None:

    met_changes['wind_c'] = 3

    dTwaterCdt = EnergyBudget(
        met_changes,
        t_changes,
        use_sed_temp=use_sed_temp,
    ).run(variables=vars)
    TwaterC = vars['TwaterC'] + dTwaterCdt * 60

    assert pytest.approx(TwaterC, tolerance) == answers('test_wind_c_3')
