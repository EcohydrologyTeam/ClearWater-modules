from typing import (
    TypedDict,
)


class Temperature(TypedDict):
    stephan_boltzmann: float
    Cp_air: float
    emissivity_water: float
    gravity: float
    a0: float
    a1: float
    a2: float
    a3: float
    a4: float
    a5: float
    a6: float
    pb: float
    Cps: float
    h2: float
    alphas: float
    richardson_option: bool


class Meteorological(TypedDict):
    TairC: float
    Q_solar: float
    TsedC: float
    eair_mb: float
    pressure_mb: float
    cloudiness: float
    wind_speed: float
    wind_a: float
    wind_b: float
    wind_c: float
    wind_kh_kw: float


DEFAULT_METEOROLOGICAL = Meteorological(
    TairC=20,
    Q_solar=400,
    TsedC=5.0,
    eair_mb=1.0,
    pressure_mb=1013.0,
    cloudiness=0.1,
    wind_speed=3.0,
    wind_a=0.3,
    wind_b=1.5,
    wind_c=1.0,
    wind_kh_kw=1.0,
)

DEFAULT_TEMPERATURE = Temperature(
    stephan_boltzmann=5.67e-8,
    Cp_air=1005,
    emissivity_water=0.97,
    gravity=-9.806,
    a0=6984.505294,
    a1=-188.903931,
    a2=2.133357675,
    a3=-1.288580973E-2,
    a4=4.393587233E-5,
    a5=-8.023923082E-8,
    a6=6.136820929E-11,
    pb=1600.0,
    Cps=1673.0,
    h2=0.1,
    alphas=0.0432,
    richardson_option=True,
)
