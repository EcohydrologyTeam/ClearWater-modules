import xarray as xr


def celcius_to_kelvin(tempc: xr.DataArray) -> xr.DataArray:
    return tempc + 273.16


def density_water(temperature: xr.Dataset) -> xr.Dataset:
        """Compute the density of water (kg/m3) as a function of water temperature (Celsius)
        Inputs: water_temperature [xr.Dataset] Water temperature in units of Celsius
        Outputs: density_water [xr.Dataset] Density of water in units of kg/m3
        """
        return 999.973 * (
            1.0
            - (
                ((temperature - 3.9863) ** 2 * (temperature + 288.9414))
                / (508929.2 * (temperature + 68.12963))
            )
        )