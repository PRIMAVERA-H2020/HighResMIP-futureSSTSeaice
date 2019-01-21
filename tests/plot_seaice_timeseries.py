import matplotlib
matplotlib.use('Agg')
import iris, os, sys, glob
import iris.coord_categorisation as icc
import iris.analysis
import iris.plot as iplt
from matplotlib.colors import BoundaryNorm  
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

dir_in = '/gws/nopw/j04//hrcm/cache/malcolm/HadISST2/1x1/processing_2018/future2/'

fice = 'future_ice_conc.nc'
pole_lat = {}; pole_lon = {}
pole_lat['npole'] = iris.Constraint(coord_values = {'latitude' : lambda l : 60 < (l.point) <= 90})
pole_lon['npole'] = iris.Constraint(coord_values = {'longitude' : lambda l : 0 < (l.point) <= 360})
pole_lat['spole'] = iris.Constraint(coord_values = {'latitude' : lambda l : -90 < (l.point) <= -50})
pole_lon['spole'] = iris.Constraint(coord_values = {'longitude' : lambda l : 0 < (l.point) <= 360})

CMIP5_ref = {'ACCESS1-0': 'CSIRO-BOM', 'ACCESS1-3': 'CSIRO-BOM', 'GFDL-CM3': 'NOAA-GFDL', 'IPSL-CM5A-LR': 'IPSL', 'IPSL-CM5A-MR': 'IPSL', 'MPI-ESM-MR': 'MPI-M', 'CNRM-CM5': 'CNRM-CERFACS', 'HadGEM2-ES':'MOHC'}
#CMIP5_ref = {'IPSL-CM5A-MR': 'IPSL'}
#CMIP5_ref = {'ACCESS1-0': 'CSIRO-BOM', 'ACCESS1-3': 'CSIRO-BOM', 'GFDL-CM3': 'NOAA-GFDL', 'IPSL-CM5A-MR': 'IPSL', 'CNRM-CM5': 'CNRM-CERFACS', 'HadGEM2-ES':'MOHC'}

dir_base = '/gws/nopw/j04/hrcm/cache/malcolm/HighResMIP/sst_forcing/processing_ocean_v2/'
dir_area = '/badc/cmip5/data/cmip5/output1/{}/{}/historical/fx/ocean/fx/r0i0p0/latest/areacello/areacello_fx_{}_historical_r0i0p0.nc'

def history_period_callback(cube, field, filename):
    # remove attributes preventing cube concatenation
    attrib = ['history','time_period']
    for att in attrib:
        try:
            del cube.attributes[att]
        except:
            pass

def guess_areas(cube):
    coords = ['latitude','longitude']
    for coord in coords:
        cube.coord(coord).bounds = None
        cube.coord(coord).guess_bounds()

    grid_areas = iris.analysis.cartography.area_weights(cube)
    return grid_areas

def history_period_callback(cube, field, filename):
    # remove attributes preventing cube concatenation
    attrib = ['history','time_period']
    for att in attrib:
        try:
            del cube.attributes[att]
        except:
            pass
    coords = ['day_of_year', 'year', 'month']
    for coord in coords:
        try:
            cube.remove_coord(coord)
        except:
            pass
    try:
        cube.coord('time').bounds = None
    except:
        pass

def highresmip_sic(pole, mon):
    dir_in = '/group_workspaces/jasmin2/primavera1/WP6/forcing/HadISST2_submit/v1.2/'
    dir_in_future = '/gws/nopw/j04/hrcm/cache/malcolm/HadISST2/1x1/processing_2018/future2/siconc/'
    dir_in_future = '/group_workspaces/jasmin2/primavera1/WP6/forcing/FutureSSTandSeaice/future3/siconc/'
    ice_ts_filename = os.path.join(dir_in_future, 'HadISST2.2_historic_ice_ts_'+pole+'.nc')
    ice_ts_future_filename = os.path.join(dir_in_future, 'HadISST2.2_future_ice_ts_v1_'+pole+'.nc')

    fice = 'siconc_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_'
    fnames = sorted(glob.glob(os.path.join(dir_in, fice+'*.nc')))
    mask = iris.load_cube(os.path.join(dir_in, 'tos_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_20150101-20151231.nc'))[0]

    if not os.path.exists(ice_ts_filename):
        ice_l = iris.cube.CubeList()
        for fn in fnames:
            ice = iris.load_cube(fn, callback = history_period_callback)
            print ice
            ice_l.append(ice)
        ice = ice_l.concatenate_cube()
        ice_pole = ice.extract(pole_lat[pole])
        mask_pole = mask.extract(pole_lat[pole])
        try:
            icc.add_year(ice_pole, 'time', name = 'year')
            icc.add_month_number(ice_pole, 'time', name = 'month')
        except:
            pass
        area = guess_areas(ice_pole)
        masked = np.where(mask_pole.data.mask == True)
        ice_pole.data[:, masked[0], masked[1]] = 0.0
        #ice_pole.data.mask = mask_pole.data.mask
        ice_ts = ice_pole.collapsed(['longitude','latitude'], iris.analysis.SUM, weights = area)
        iris.save(ice_ts, ice_ts_filename)
    else:
        ice_ts = iris.load_cube(ice_ts_filename)

    if mon != 0:
        month_constraint = iris.Constraint(coord_values = {'month': lambda l : l == mon})
        ice_ts1 = ice_ts.extract(month_constraint)
    else:
        ice_ts1 = ice_ts
    print ice_ts1
    ice_ts2 = ice_ts1.aggregated_by(['year', 'month'], iris.analysis.MEAN)

    fice = 'future_siconc_????_025_daily_v1'
    fnames = sorted(glob.glob(os.path.join(dir_in_future, fice+'.nc')))
    
    if not os.path.exists(ice_ts_future_filename):

        ice_fut_l = iris.cube.CubeList()
        for fn in fnames:
            ice = iris.load_cube(fn, callback = history_period_callback)
            print ice
            ice_fut_l.append(ice)
        ice_fut = ice_fut_l.concatenate_cube()
        ice_fut_pole = ice_fut.extract(pole_lat[pole])
        try:
            icc.add_year(ice_fut_pole, 'time', name = 'year')
            icc.add_month_number(ice_fut_pole, 'time', name = 'month')
        except:
            pass
        area = guess_areas(ice_fut_pole)
        ice_fut_ts = ice_fut_pole.collapsed(['longitude','latitude'], iris.analysis.SUM, weights = area)
        iris.save(ice_fut_ts, ice_ts_future_filename)
    else:
        ice_fut_ts = iris.load_cube(ice_ts_future_filename)

    if mon != 0:
        month_constraint = iris.Constraint(coord_values = {'month': lambda l : l == mon})
        ice_fut_ts1 = ice_fut_ts.extract(month_constraint)
    else:
        ice_fut_ts1 = ice_fut_ts
    print ice_fut_ts1
    ice_fut_ts2 = ice_fut_ts1.aggregated_by(['year', 'month'], iris.analysis.MEAN)
    print ice_fut_ts2
    return ice_ts2, ice_fut_ts2


def reference_ts(pole):
    model = 'HadGEM2-ES'
    fname = os.path.join(dir_base, model+'_monthly_1950_2100_sic.nc')
    fname_area = dir_area.format(CMIP5_ref[model], model, model)
    ice = iris.load_cube(fname)
    ice_pole = ice.extract(pole_lat[pole])
    try:
        icc.add_year(ice_pole, 'time', name = 'year')
        icc.add_month_number(ice_pole, 'time', name = 'month')
    except:
        pass
    area = guess_areas(ice_pole)
    ice_ts = ice_pole.collapsed(['longitude','latitude'], iris.analysis.SUM, weights = area)
    return ice_ts


def cmip5_ice(model, mon, pole = 'npole'):
    ref_ts = reference_ts(pole)
    month_constraint = iris.Constraint(coord_values = {'month': lambda l : l == mon})

    fname = os.path.join(dir_base, model+'_monthly_1950_2100_sic.nc')

    fname_out = fname[:-3]+'_monthly_timeseries_'+pole+'.nc'
    fname_area = dir_area.format(CMIP5_ref[model], model, model)
    print fname_area

    if not os.path.exists(fname_out):
        ice_ts = ref_ts.copy()

        with Dataset(fname, 'r') as fh, Dataset(fname_area, 'r') as fa:
            sic = fh.variables['sic']
            area = fa.variables['areacello']
        
            print sic.shape
            for t in range(sic.shape[0]):
            #miss = np.where(sic[t, :, :] < 15.0)
            #sic[t, miss[0], miss[1]] = 0.0
                sic_area = sic[t, :, :] * area[:, :]
                nhemi = sic.shape[1] // 2
                if pole == 'npole':
                    if 'MPI' in model:
                        sic_total = np.sum(sic_area[:nhemi, :])
                    else:
                        sic_total = np.sum(sic_area[nhemi:, :])
                else:
                    if 'MPI' in model:
                        sic_total = np.sum(sic_area[nhemi:, :])
                    else:
                        sic_total = np.sum(sic_area[:nhemi, :])
                ice_ts.data[t] = sic_total
    
            iris.save(ice_ts, fname_out)
    else:
        ice_ts = iris.load_cube(fname_out)

    if mon != 0:
        ref_ts1 = ref_ts.extract(month_constraint)
        ice_ts1 = ice_ts.extract(month_constraint)
    else:
        ref_ts1 = ref_ts
        ice_ts1 = ice_ts

    return ref_ts1, ice_ts1


if __name__ == '__main__':

    months = range(1,13)
    years = range(2015, 2050)
    year_constraint = iris.Constraint(coord_values = {'year': lambda l : years[0] <= l <= years[-1]})


    for pole in ['npole','spole']:
        for mon in range(1,13):
            if mon == 0:
                time_coord = 'time'
                xlim = [60000, 70000]
            else:
                time_coord = 'year'
                xlim = [1950, 2050]
            hrmip_hist, hrmip_fut = highresmip_sic(pole, mon)
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            for im, model in enumerate(CMIP5_ref):
                ice_ref, ice_fut = cmip5_ice(model, mon, pole = pole)
                #print ice_ref
                #print ice_ref.coord('time')
                if im == 0:
                    plt.plot(ice_ref.coord(time_coord).points, ice_ref.data, label = 'HadGEM2-ES')
                    #iplt.plot(ice_ref)
                plt.plot(ice_fut.coord(time_coord).points, ice_fut.data, label = model)
                #iplt.plot(ice_fut)
                plt.title('Month '+str(mon)+', '+pole)
            plt.plot(hrmip_hist.coord(time_coord).points, hrmip_hist.data, label = 'HighResMIP_past', color = 'black')
            plt.plot(hrmip_fut.coord(time_coord).points, hrmip_fut.data, label = 'HighResMIP_fut', color = 'black', linestyle = '--')
            plt.legend(loc = 'best', ncol=2)
            ax.set_xlim(xlim[0], xlim[1])
            plt.savefig('./cmip5_icearea_'+str(mon)+'_'+pole+'_v3.png')
    plt.show()

