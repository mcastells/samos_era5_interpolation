#https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

from datetime import timedelta

import cdsapi
import numpy as np
import pygrib
from metpy.units import masked_array, units
from numpy import ma

c = cdsapi.Client()

class Era5:
    def __init__(self, gribs, dt_start, dt_end) -> None:
        self.timesteps = []
        self.lons = None
        self.lats = None
        self.lons_ocean = None
        self.lats_ocean = None

        #self.wind_10m = []
        self.u_wind_10m = []
        self.v_wind_10m = []
        self.mslp = []
        self.temp_2m = []
        self.dewpoint_2m = []
        self.sst = []

        for grb in gribs:
            if grb.validDate >= dt_start and grb.validDate <= dt_end:
                if grb.validDate not in self.timesteps:
                    self.timesteps.append(grb.validDate)
                if '10 metre U wind component' in grb.parameterName: 
                    self.u_wind_10m.append(grb.values)
                elif '10 metre V wind component' in grb.parameterName: 
                    self.v_wind_10m.append(grb.values)
                elif 'Mean sea level pressure' in grb.parameterName: 
                    self.mslp.append(grb.values)
                elif '2 metre temperature' in grb.parameterName: 
                    self.temp_2m.append(grb.values)
                elif '2 metre dewpoint' in grb.parameterName: 
                    self.dewpoint_2m.append(grb.values)
                elif grb.parameterName == 'Sea surface temperature': 
                    self.sst.append(grb.values)
                    if self.lats_ocean is None and self.lons_ocean is None:
                        self.lats_ocean, self.lons_ocean = grb.latlons()
                else:
                    raise(Exception(f'Match not found for parameter name: {grb.parameterName}'))
        self.timesteps = np.array(self.timesteps)
        ntimesteps = len(self.timesteps)
        
        self.lats, self.lons = grb.latlons()  # get the lats and lons for the grid.
        nlats_atmos = self.lats.shape[0]
        nlons_atmos = self.lats.shape[1]
        if self.lats_ocean is not None:
            nlats_ocean = self.lats_ocean[0]
            nlons_ocean = self.lats_ocean[1]
        else:
            nlats_ocean = 0
            nlons_ocean = 0

        def invalidate_data_missing_timesteps(data,grid):
            return data if len(data) == ntimesteps else ma.masked_all(shape=(ntimesteps,nlats_atmos,nlons_atmos) if grid=='atmos' else (ntimesteps,nlats_ocean,nlons_ocean))

        def add_units(data, unit_name, grid='atmos'):
            return masked_array(invalidate_data_missing_timesteps(data,grid),units(unit_name))
        
        self.u_wind_10m = add_units(self.u_wind_10m, 'm/s')
        self.v_wind_10m = add_units(self.v_wind_10m, 'm/s')
        self.mslp = add_units(self.mslp, 'Pa')
        self.temp_2m = add_units(self.temp_2m, 'K')
        self.dewpoint_2m = add_units(self.dewpoint_2m, 'K')
        self.sst = add_units(self.sst, 'K')

        #self.wind_10m = (self.u_wind_10m**2 + self.v_wind_10m**2)**0.5


def download(dt_start, dt_end, data_extent, filename):
    years = []
    months = []
    days = []
    date_iter = dt_start.replace(hour=0)
    while date_iter <= dt_end + timedelta(days=1):
        if str(date_iter.year) not in years:
            years.append(str(date_iter.year))
        if str(date_iter.month) not in months:
            months.append(str(date_iter.month))
        if str(date_iter.day) not in days:
            days.append(str(date_iter.day))
        date_iter += timedelta(days=1)

    era5_variables = ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature', '2m_dewpoint', 'mean_sea_level_pressure', 'sea_surface_temperature']

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'grib',
            'variable': era5_variables,
            'year': years,
            'month': months,
            'day': days,
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'area': [
                data_extent[3], data_extent[0], data_extent[2],
                data_extent[1],
            ],
        },
        filename)
    gribs = pygrib.open(filename)

    for grb in gribs[1:len(era5_variables)]:
        print(grb.parameterName)
        print(f'\t{grb.units}')
    gribs.rewind()

    return gribs