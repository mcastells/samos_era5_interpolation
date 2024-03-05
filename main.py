import datetime
from pathlib import Path

import numpy as np
import numpy.ma as ma
import pandas
import pygrib
import xarray as xr
from metpy.interpolate import interpolate_to_points
from termcolor import colored
from tqdm import tqdm, trange

import era5

samos_input_dir = Path('/Net/mdc/marineflux/data/SAMOS_Fluxes_2024-02-13/B23')
era5_dir = Path('~/era5_data/')

def interpolate_timestep(previoushour, nexthour, minute):
    weight = minute/60.
    return previoushour * (1-weight) + nexthour * weight

if __name__ == '__main__':
    ncpaths = sorted(list(samos_input_dir.rglob('*.nc')))

    for path in tqdm(ncpaths, desc=f'Processing {samos_input_dir.name} files', leave=False):
        
        if path.is_file() and path.suffix == '.nc':

            tqdm.write(f'Opening {colored(path, "blue")}')
            samos_ds = xr.open_dataset(path)

            timepadding = 1
            synoptic_dt_start = pandas.to_datetime(samos_ds.time.data[0]) - datetime.timedelta(hours=timepadding)
            synoptic_dt_end = pandas.to_datetime(samos_ds.time.data[-1]) + datetime.timedelta(hours=timepadding)

            extent_padding = 1
            extent = [ma.min(samos_ds.lon[:]) - extent_padding, ma.max(samos_ds.lon[:]) + extent_padding, ma.min(samos_ds.lat[:]) - extent_padding, ma.max(samos_ds.lat[:]) + extent_padding]
            
            grib_filepath = f'{era5_dir / path.stem}_{datetime.datetime.strftime(synoptic_dt_start,"%Y%m%d-%H%M")}_{datetime.datetime.strftime(synoptic_dt_end,"%Y%m%d-%H%M")}_era5.grib'
            grib_p = Path(grib_filepath)
            if not grib_p.is_file():
                # if a grib file with the specific datetime range doesn't exist, look for grib files with a wider datetime range
                for p in list(era5_dir.glob(f'{path.stem}*.grib')):
                    p_split = p.name.split('_')
                    grib_start = datetime.datetime.strptime(p_split[-3],"%Y%m%d-%H%M")
                    grib_end = datetime.datetime.strptime(p_split[-2],"%Y%m%d-%H%M")
                    if grib_start <= synoptic_dt_start and grib_end >= synoptic_dt_end:
                        grib_p = p
                        break

            try:
                grbs = pygrib.open(str(grib_p))
            except Exception as error:
                grbs = era5.download(synoptic_dt_start, synoptic_dt_end, extent, grib_filepath)

            era5_ds = era5.Era5(grbs,synoptic_dt_start,synoptic_dt_end)

            era5_interpolated_to_ship = {
                'uwind_10m': np.empty(shape=len(samos_ds.time)),
                'vwind_10m': np.empty(shape=len(samos_ds.time)),
                'mslp': np.empty(shape=len(samos_ds.time)),
                'temp_2m': np.empty(shape=len(samos_ds.time)),
                'dewpoint_2m': np.empty(shape=len(samos_ds.time)),
                'sst': np.empty(shape=len(samos_ds.time))
            }

            syn_timesteps = era5_ds.timesteps

            for era5_var_name in era5_interpolated_to_ship.keys():
                if era5_var_name == 'sst':
                    syn_lons = era5_ds.lons_ocean
                    syn_lats = era5_ds.lats_ocean
                else:
                    syn_lons = era5_ds.lons
                    syn_lats = era5_ds.lats
                
                era5_var = getattr(era5_ds, era5_var_name)

                syn_lons[np.where(syn_lons < 0)] += 360

                points = np.array((syn_lons.flatten(),syn_lats.flatten())).T

                for i in trange(len(samos_ds.time), desc='Interpolating minute ERA5 variables', leave=False):
                    dt = samos_ds.time[i]
                    if dt.minute == 0:
                        timestep = np.where(syn_timesteps == dt)[0][0]
                        era5var_atship = interpolate_to_points(points=points,values=era5_var[timestep].m.flatten(),xi=(samos_ds.lon[i],samos_ds.lat[i]))
                    else:
                        previoushour = np.where(syn_timesteps < dt)[0][-1]
                        era5var_attime = interpolate_timestep(era5_var[previoushour].m,era5_var[previoushour+1].m,dt.minute)
                        era5var_atship = interpolate_to_points(points=points,values=era5var_attime.flatten(),xi=(samos_ds.lon[i],samos_ds.lat[i]))
                    era5_interpolated_to_ship[era5_var_name][i] = float(era5var_atship)

            print(era5_interpolated_to_ship)

            # calculate fluxes using these interpolated values