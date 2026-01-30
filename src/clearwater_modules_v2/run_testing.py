from datetime import datetime, timedelta
from clearwater_modules_v2.processes.temperature import Temperature
from clearwater_data.variables.registry import VariableRegistry
from clearwater_data.variables.float import FloatVariable

# Instantiate the Variable Registry
# This is were you will define data the process will access
variable_registry = VariableRegistry()
# To add data to the registry you need a variable
water_temp = FloatVariable(20)
surf_area = FloatVariable(1)
vol = FloatVariable(1)

air_temp_c = FloatVariable(50.0)  # was 20.0 in v1 example
q_solar = FloatVariable(900.0)  # was 400.0 in v1 example
sed_temp_c = FloatVariable(5.0)

eair_mb = FloatVariable(
    1.0
)  # v1 example uses 1.0. This is the atmospheric vapor pressure.
# This value looks like a low value; a google search shows typical values around 10mb

pressure_mb = FloatVariable(1013.0)
cloudiness_frct = FloatVariable(0.1)
wind_spd = FloatVariable(3.0)

sed_temp_c = FloatVariable(50.0)  # was 5.0 in v1 example
sed_thick_m = FloatVariable(
    1.0
)  # looks like the the thickness was 0.1m in the v1 example

# define wind function parameters as floats
wind_a_userInput = 0.3
wind_b_userInput = 1.5
wind_c_userInput = 1.0
# wind_kh_kw = FloatVariable(1.0) #this is the air diffusivity ratio in v1 example
# this parameter is used in v2 temperature process and it defaults to 1.0 so no need to define it here

sediment_density_userInput = 1.6  # g/cm3 (1600 kg/m3 was used in v1 example);
# default value in v2 temperature process is 1.67 g/cm3

sediment_specific_heat_userInput = 1673.0  # J/kg/K was used in v1 example
# default value in v2 temperature process is 1000.0 J/kg/K

sediment_diffusivity_userInput = 0.0432  # m^2/s #used in v1 example
# default value in v2 temperature process is 0.0061 m^2/s

# define time step for the process
time_step_userInput = timedelta(seconds=1)  # also tested with 300sec (5min)
# Now you can register this
variable_registry.register("water_temperature", water_temp)
variable_registry.register("wetted_surface_area", surf_area)
variable_registry.register("volume", vol)
variable_registry.register("cloudiness", cloudiness_frct)
variable_registry.register("air_temperature", air_temp_c)
variable_registry.register("solar_radiation", q_solar)
variable_registry.register("wind_speed", wind_spd)
variable_registry.register("atmospheric_pressure", pressure_mb)
variable_registry.register("atmospheric_vapor_pressure", eair_mb)
variable_registry.register("sediment_temperature", sed_temp_c)
variable_registry.register("sediment_thickness", sed_thick_m)
process = Temperature(
    wind_a=wind_a_userInput,
    wind_b=wind_b_userInput,
    wind_c=wind_c_userInput,
    sediment_density=sediment_density_userInput,
    sediment_specific_heat=sediment_specific_heat_userInput,
    sediment_diffusivity=sediment_diffusivity_userInput,
    time_step=time_step_userInput,
)
date_time = datetime(2026, 1, 1, 0, 0, 0)
variable_registry.get("water_temperature")
process.run(date_time, variable_registry)
variable_registry.get("water_temperature")
temp_time_series = []

for _ in range(84600):
    temp_at_start = variable_registry.get("water_temperature")
    process.run(date_time, variable_registry)
    temp_at_end = variable_registry.get("water_temperature")
    temp_pair = (temp_at_start, temp_at_end)
    temp_time_series.append(temp_pair)

temp_time_series
