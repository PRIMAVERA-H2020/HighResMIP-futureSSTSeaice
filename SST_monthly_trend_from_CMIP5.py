#!/bin/env python2.7
'''
NAME:
    SST_trend_from_CMIP5

DESCRIPTION:
    Take the monthly surface temperature from each CMIP5 RCP8/5 simulation from 2010 to 2100
    calculate a trend at each point
    ensemble mean the trends from each model
    regrid to the HadISST grid
    use this for the future trend for the SSTs for HighResMIP, rather than extrapolating
        the observed trend

    Anomalise:
        calculate time mean, monthly mean 2D field (month, y, x)
        calculate 2D spatial linear trend by month (year, month, y, x)
    
    Issues: 
        In order to add back monthly anomaly, need to know where ice-edge is
        Antarctic - SST is not warming so much here (and sea-ice is currently increasing a little). 
        Hence need to take care about trend here. May want to extrapolate trend from ~1982 since in obs there is a jump
        lack of data pre-satellite era
        Issues around ice edge

USAGE: 
    Execute as a script from command line

AUTHOR:
    Malcolm Roberts (hadom)
    with various subroutines from Chris Roberts (hadrr)    

LAST MODIFIED:
    2015-06-22

'''

import iris.coord_categorisation as icc
import os
import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.analysis
import iris.experimental as ie
import iris.analysis as ia
import scipy.linalg
import iris.io
import subprocess
import iris.plot as iplt
from scipy.stats import norm
import datetime
import sys
import glob
import cPickle as pickle #cPickle is an optimised version of Pickle and is O(100) times faster
from shutil import move
#import cf_units
from scipy import ndimage as nd
import gridsmoother  # additional package from M.Mizielinski (Met Office)
from multiprocessing import Process, Pool
import sst_future

# need to run from the output directory
# assumes that years in file names are real years

BADC_CMIP5_path1 = '/badc/cmip5/data/cmip5/output1'
BADC_CMIP5_path_hist = {}
BADC_CMIP5_path_hist['tos'] = 'historical/mon/ocean/Omon/r1i1p1/latest/tos'
BADC_CMIP5_path_hist['sic'] = 'historical/mon/seaIce/OImon/r1i1p1/latest/sic'
BADC_CMIP5_path_rcp85 = {}
BADC_CMIP5_path_rcp85['tos'] = 'rcp85/mon/ocean/Omon/r1i1p1/latest/tos'
BADC_CMIP5_path_rcp85['sic'] = 'rcp85/mon/seaIce/OImon/r1i1p1/latest/sic'
BADC_CMIP5_path_mask1 = 'historical/fx/ocean/fx/r0i0p0/latest/sftof/'
BADC_CMIP5_path_mask2 = 'rcp85/fx/ocean/fx/r0i0p0/latest/sftof/'

variable1 = {}; variable2 = {}
variable1['tos'] = 'sea_surface_temperature'
variable2['tos'] = 'surface_temperature'
variable1['sic'] = 'sea_ice_area_fraction'
variable2['sic'] = 'sea_ice_area_fraction'
YEAR_LAST_RCP85 = '2120'

HadISST_ref = '/gws/nopw/j04/hrcm/cache/malcolm/HadISST2.1.0.0/HadISST.2.1.0.0_single_sst.nc'

DATADIR = '/gws/nopw/j04/hrcm/cache/malcolm/HighResMIP/sst_forcing'
savedir = os.path.join(DATADIR, 'processing_ocean_v3')
savedir_v0 = os.path.join(DATADIR, 'processing_ocean_v2')
if not os.path.exists(savedir): os.mkdir(savedir)

YEARS_after_1870 = iris.Constraint(coord_values = {'year' : lambda l : l >= 1870})
YEARS_after_1900 = iris.Constraint(coord_values = {'year' : lambda l : l >= 1900})
YEARS_after_1910 = iris.Constraint(coord_values = {'year' : lambda l : l >= 1910})
YEARS_after_1950 = iris.Constraint(coord_values = {'year' : lambda l : l >= 1950})
YEARS_after_2010 = iris.Constraint(coord_values = {'year' : lambda l : l >= 2010})
YEARS_between_1970_2010 = iris.Constraint(coord_values = {'year' : lambda l : 1970 <= l <= 2010})
YEARS_between_1950_2014 = iris.Constraint(coord_values = {'year' : lambda l : 1950 <= l <= 2014})
YEARS_between_2004_2024 = iris.Constraint(coord_values = {'year' : lambda l : 2004 <= l <= 2024})
YEARS_between_2012_2016 = iris.Constraint(coord_values = {'year' : lambda l : 2012 <= l <= 2016})
YEARS_between_2010_2100 = iris.Constraint(coord_values = {'year' : lambda l : 2010 <= l <= 2100})
YEARS_between_1950_2100 = iris.Constraint(coord_values = {'year' : lambda l : 1950 <= l <= 2100})
YEARS_before_1950 = iris.Constraint(coord_values = {'year' : lambda l : l <= 1949})

verbose = True
# chunks of data to get at once
interval = 30
colourname = ['red','blue','green','orange','aqua','lime','sienna','black']

def CMIP5_datasets():
    # dictionary of CMIP5 centre/model in order to get the path to the data
#    CMIP5_ref = {'CCCma': 'CanESM2', \
#                 'CNRM-CERFACS': 'CNRM-CM5', 'CSIRO-BOM': 'ACCESS1-0', 'CSIRO-BOM': 'ACCESS1-3', \
#                 'IPSL':'IPSL-CM5A-MR', 'MPI-M':'MPI-ESM-MR', 'MRI':'MRI-CGCM3',\
#                 'NCAR':'CCSM4', 'NOAA-GFDL':'GFDL-CM3', 'NSF-DOE-NCAR': 'CESM1-CAM5', 'MOHC': 'HadGEM2-ES'}
    # subset v3
    #CMIP5_ref = {'CNRM-CERFACS': 'CNRM-CM5', \
    #             'IPSL':'IPSL-CM5A-MR', \
    #             'NSF-DOE-NCAR': 'CESM1-CAM5', 'MOHC': 'HadGEM2-ES', 'NOAA-GFDL':'GFDL-CM3'}

    # Massonet 2012 has:
    CMIP5_ref = {'ACCESS1-0': 'CSIRO-BOM', 'ACCESS1-3': 'CSIRO-BOM', 'GFDL-CM3': 'NOAA-GFDL', 'IPSL-CM5A-LR': 'IPSL', 'IPSL-CM5A-MR': 'IPSL', 'MPI-ESM-MR': 'MPI-M', 'CNRM-CM5': 'CNRM-CERFACS', 'HadGEM2-ES': 'MOHC'}
    # find IPSL-CM5A-LR and MPI-ESM-MR are outliers in ice retreat

#    CMIP5_ref = {'ACCESS1-0': 'CSIRO-BOM', 'ACCESS1-3': 'CSIRO-BOM', 'GFDL-CM3': 'NOAA-GFDL', 'IPSL-CM5A-MR': 'IPSL', 'CNRM-CM5': 'CNRM-CERFACS', 'HadGEM2-ES':'MOHC'}
    # failed# 'CMCC': 'CMCC-CM', - different units (days since)
    # failed 'MOHC': 'HadGEM2-ES' - terrible trouble with time coords
    #  made this by hand and put into directory
    return CMIP5_ref

def form_BADC_path_rcp85(runid, var = 'tos'):
    return os.path.join(BADC_CMIP5_path1, runid, BADC_CMIP5_path_rcp85[var])

def form_BADC_path_hist(runid, var = 'tos'):
    return os.path.join(BADC_CMIP5_path1, runid, BADC_CMIP5_path_hist[var])

def form_BADC_path_mask1(runid):
    return os.path.join(BADC_CMIP5_path1, runid, BADC_CMIP5_path_mask1)

def form_BADC_path_mask2(runid):
    return os.path.join(BADC_CMIP5_path1, runid, BADC_CMIP5_path_mask2)

def fill_in_masked_data(cube):
    # now fill the data across the land points by using the nearest value - so no new extrema created
    cube_filled = cube.copy()
    for k in range(cube.shape[0]):
        data = cube.data[k,...]
        data[np.where(cube.data.mask[k,...] == True)] = np.nan
#        y_sampling = cube_filled.coord('latitude')
#        x_sampling = cube_filled.coord('longitude')
        # give some differential between grid spacing in x and y - particularly for Arctic, want
        # to use same latitude more than euclidean distance
        sampling = [1.0, 0.6]
        data_fill = fill(data, sampling=sampling)
        cube_filled.data[k,...] = data_fill
    return cube_filled

def fill(data, invalid=None, sampling = None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid') 
    by the value of the nearest valid data cell

    Input:
        data:    numpy array of any dimension
        invalid: a binary array of same shape as 'data'. True cells set where data
                 value should be replaced.
                 If None (default), use: invalid  = np.isnan(data)

    Output: 
        Return a filled array. 
    Source: 
        http://stackoverflow.com/questions/3662361/fill-in-missing-values-with-nearest-neighbour-in-python-numpy-masked-arrays
    """

    if invalid is None: invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True, sampling = sampling)
    return data[tuple(ind)]

def fix_cubelist(cubelist):
    u = cubelist[1].coord('time').units
    points = cubelist[1].coord('time').points
    #u = time_coord
    #cubelist[0].coord('time').points = cubelist[0].coord('time').points.astype(np.float64)
    for cube in cubelist:
        cube.coord('time').attributes['comment'] = 'overwritten'
        #cube.coord('time').convert_units(u)
        cube.coord('time').units = u
        cube.coord('time').bounds = None
        cube.coord('time').points = points
    return cubelist

def add_mask(cube):
        if not "MaskedArray" in str(type(cube.data)):
                cube.data = np.ma.MaskedArray(cube.data, shrink=False)
        cube.data.mask = False
        return cube

def space_filter_and_save(cube, ic):
    cube_filter = filter_data_boxcar(cube, 20.0, 10.0, 'hadisst2','025x025', return_smooth=True, masked=True)
    iris.save(cube_filter, os.path.join(savedir, 'cmip5_mean_trend_filtered_modelmean_masked_'+str(ic)+'.nc'), unlimited_dimensions = [])

def remove_nemo_wrap(file_in, file_out):
    # remove the extra two columns of NEMO data (which are masked and hence cdo maps the mask)
    if os.path.exists(file_out): os.remove(file_out)
    x_dimension = 'i'
    cmd = 'ncks -d'+x_dimension+',1,180 '+file_in+' '+file_out
    print cmd
    os.system(cmd)

def cmip5_callback(cube, field, filename):
    # add a categorical year coordinate to the cube
    attrib_remove = ['history','creation_date','references','mo_runid','tracking_id','cmor_version',\
                     'table_id','forcing', 'CDO', 'CDI']
    for attrib in attrib_remove:
        try:
            del cube.attributes[attrib]
        except:
            pass

def fix_time_coord(cube_in, cube_time):
    time_coord = cube_time.coord('time')
    u = cube_time.coord('time').units
    #cube_in.coord('time').attributes = cube_time('time').attributes
    cube_in.coord('time').points = time_coord.points
    cube_in.coord('time').units = u
    return cube_in
    
def cmip5_callback_delall(cube, field, filename):
    # add a categorical year coordinate to the cube
    attrib_remove = ['history','creation_date','references','mo_runid','tracking_id','cmor_version',\
                     'table_id','forcing']
    attrib_remove = cube.attributes.copy()
    for attrib in attrib_remove:
        try:
            del cube.attributes[attrib]
        except:
            pass
    cube.coord('time').attributes = {}
    cube.coord('time').bounds = None

def use_ncrcat_to_merge(files, final_file):
    if os.path.exists(final_file):
        os.remove(final_file)
    cmd = 'ncrcat '+' '.join(files)+' '+final_file
    print cmd
    os.system(cmd)

def fix_fill_value_with_nco(fname):
    cmd = 'ncatted -h -O -a _FillValue,unknown,o,d,1.0e20 '+fname
    os.system(cmd)

#def nco_add_dimension(file):
#    cmd = 'ncap2 -s 'defdim("model",1);model[model]=1;model@long_name="Model"''

def filter_data_boxcar(cube_in, lon_width, lat_width, field, ocean_resol, return_smooth=False, masked = False):
        
        gs = gridsmoother.GridSmoother()
        filename= "boxcar_%s_%s_x_%s.json_%s.gz" % (field, lon_width, lat_width, ocean_resol)
        if os.path.exists(filename):
            gs.load(filename)
        else:
            ndims = len(cube_in.shape)
            if ndims == 3:
                cube2d = cube_in[0]
            elif ndims == 2:
                cube2d = cube_in
            else:
                sys.exit('Dims in cube is not 3 or 2 for spatial filtering')
            gs.build_boxcar_lat_lon(cube2d, lat_width, lon_width)
            gs.save(filename)
                
        smoothed_cube = gs.smooth_3d_cube(cube_in, masked=masked)
        if return_smooth:
            return smoothed_cube
        else:
            return cube_in - smoothed_cube

def find_and_merge_historical_rcp85(CMIP5_ref, var = 'tos'):
    for model in CMIP5_ref:
        run = CMIP5_ref[model]
        print run
        runid = run+'/'+model
        datapath_hist = form_BADC_path_hist(runid, var = var)
        datapath_rcp85 = form_BADC_path_rcp85(runid, var = var)
        print runid, datapath_hist, datapath_rcp85
        output_file = os.path.join(savedir, model+'_monthly_1950_'+YEAR_LAST_RCP85+'_'+var+'.nc')
        if os.path.exists(output_file):
            continue
        datafile = glob.glob(os.path.join(datapath_hist, var+'_*mon_'+model+'_historical_r1i1p1_*'))
        print 'datafile ',datafile
        # read all files into cube, 
        cube_tmp = iris.load(datafile, variable1[var], callback = cmip5_callback_delall)
        if len(cube_tmp) == 0:
            cube_tmp = iris.load(datafile, variable2[var], callback = cmip5_callback_delall)
        sst_cube_full = cube_tmp.concatenate_cube()
        print 'sst_cube_full ',sst_cube_full
        tmp_file = os.path.join(savedir, model+'_monthly_hist_tmp_'+var+'.nc')
        iris.save(sst_cube_full, tmp_file, unlimited_dimensions = ['time'], fill_value = 1.0e20)
        
        # add year coordinate
        icc.add_year(sst_cube_full, 'time', name = 'year')
        end_year = sst_cube_full.coord('year').points.max()
        print 'end_year ',end_year

        YEARS_from_1950 = iris.Constraint(coord_values = {'year' : lambda l : l >= 1950})
        # extract years from 1950-last year, noting last year included
        sst_cube_from1950 = sst_cube_full.extract(YEARS_from_1950)
        hist_file = os.path.join(savedir, model+'_monthly_hist_'+var+'.nc')
        iris.save(sst_cube_from1950, hist_file, unlimited_dimensions = ['time'], fill_value = 1.0e20)
            
        # then read the rcp85 directory, 
        datafile = glob.glob(os.path.join(datapath_rcp85, var+'_*mon_'+model+'_rcp85_r1i1p1_*'))
        print 'datafile for rcp85 ',datafile
        # limit the files to those not past YEAR_LAST_RCP85 (2120)
        datafile_in_period = []
        for fname in datafile:
            name = os.path.basename(fname)
            years = name.split('_')[-1]
            if int(years[0:4]) <= int(YEAR_LAST_RCP85):
                datafile_in_period.append(fname)
        datafile_in_period = sorted(datafile_in_period)
        print 'datafile_in_period ',datafile_in_period
        # concatenate files, need to try ncrcat instead
        rcp85_file = os.path.join(savedir, model+'_monthly_rcp85_'+var+'.nc')
        if not 'HadGEM' in model: # need to do HadGEM by hand due to duplicates
            use_ncrcat_to_merge(datafile_in_period, rcp85_file)

        sst_cube_tmp = iris.load(rcp85_file, variable1[var], callback = cmip5_callback_delall)
        if len(sst_cube_tmp) == 0:
            sst_cube_tmp = iris.load(rcp85_file, variable2[var], callback = cmip5_callback_delall)
        sst_cube_full = sst_cube_tmp.concatenate_cube()
        # add year coordinate
        icc.add_year(sst_cube_full, 'time', name = 'year')
        # find first year of this data
        start_year = sst_cube_full.coord('year').points.min()

        # make constraint from the year after the hist cube to YEAR_LAST_RCP85 (2120)
        YEARS_between_histend_2120 = iris.Constraint(coord_values = {'year' : lambda l : int(end_year)+1 <= l <= int(YEAR_LAST_RCP85)})
        sst_cube_to2120 = sst_cube_full.extract(YEARS_between_histend_2120)
        
        # remove the year coordinate (prevents cube merge)
        for c in [sst_cube_from1950, sst_cube_to2120]:
            try:
                c.remove_coord('year')
            except:
                pass
        #print sst_cube_from1950
        #print sst_cube_to2120

        # make a cube list, and then merge these 2 cubes together
        # to obtain a cube 1950-2100
        full_cube_list = iris.cube.CubeList()
        full_cube_list.append(sst_cube_from1950)
        full_cube_list.append(sst_cube_to2120)
        iris.util.unify_time_units(full_cube_list)
        #print 'full_cube_list ',full_cube_list

        cube_files = []
        if not os.path.exists(savedir): os.makedirs(savedir)
        for ic, cubef in enumerate(full_cube_list):
            #print cubef.coord('time')
            fout = os.path.join(savedir,model+'_monthly_1950_'+str(ic)+'_'+var+'.nc')
            if 'HadGEM3' in model and ic == 1:
                print 'skip the duplicate december'
                cubef = cubef[1:]
            iris.save(cubef, fout, unlimited_dimensions = ['time'], fill_value = 1.0e20)
            cube_files.append(fout)
        use_ncrcat_to_merge(cube_files, output_file)
#        if run == 'IPSL':
#            # need to remove wrap rows from NEMO
#            output_tmp = output_file[:-3]+'_tmp.nc'
#            remove_nemo_wrap(output_file, output_tmp)
#            move(output_tmp, output_file)

        #remap_with_cdo(output_file, output_file[:-3]+'_1x1.nc')
  
def remove_monthly_avg(cube):
    '''
    Removes monthly mean from each mon/ann value
    '''
    sst_month = cube.aggregated_by('month', iris.analysis.MEAN)
    print 'sst_month ',sst_month
    print 'cube ',cube
    print 'month number ',cube.coord('month')
    print 'month number ',len(cube.coord('month').points)
    monthly_anom = cube.copy()
    for m in range(0,11+1):
        inds = np.where(cube.coord('month').points == m+1)[0]
        for i in inds:
            monthly_anom.data[i] -= sst_month.data[m]
    
    return monthly_anom,sst_month
                           
def make_trend_files_for_models(model, CMIP5_ref):
    # read the merged monthly SST for 1950-2100 (datafile)
    # extract data between 2012_2016 (sst_cube)
    #   aggregate this by year (sst_cube_year)
    #   and aggregate by month (monthly_mean)
    #   subtract this monthly_mean from each year (yearly_monthly_trend) - should be called delta
    #   save as allyrs_sst.nc
    # read the land sea mask
    # mask the yearly_monthly_trend file with this mask 
    #   save as _allyrs_masked_sst.nc
    # 

    global_temp_change = {}  
    run = CMIP5_ref[model]
    print run
    runid = run+'/'+model
    print runid
    datafile = os.path.join(savedir, model+'_monthly_1950_2120_tos.nc')
    print datafile
    if not os.path.exists(savedir): os.makedirs(savedir)
        
    file_trend = os.path.join(savedir,model+'_month_annual_1950_2120_delta.nc')
    if not os.path.exists(file_trend[:-3]+'_allyrs_masked_sst.nc'):
    
        # Read SST data and calculate annual means
        # the callback method only works well is there are full years of 12 months data
        print 'start analysis'
        cube_tmp = iris.load(datafile, variable1['tos'], callback = cmip5_callback)
        if len(cube_tmp) == 0:
            cube_tmp = iris.load(datafile, variable2['tos'], callback = cmip5_callback)
        for c in cube_tmp:
            print c.coords('time')
            print 
        print cube_tmp

        sst_cube_full = cube_tmp.concatenate_cube()
        icc.add_year(sst_cube_full, 'time', name = 'year')
        icc.add_month_number(sst_cube_full, 'time', name = 'month')
            
        # convert to SST if surface temperature
        #if sst_cube_full.data[0,...].max() >= 100.: sst_cube_full -= 273.15
            
        # now the surface temperature over sea-ice is below -1.8, so reset it to -1.8
        #sst_cube_full.data[sst_cube_full.data < -1.8] = -1.8
        print 'sst_cube_full ',sst_cube_full
            
        # extract all years 10 years either side of 2015 (where we want delta to be about zero)
        year_centre = 2015
        YEARS_between = iris.Constraint(coord_values = {'year' : lambda l : year_centre - 10 <= l <= year_centre + 10})
            
        sst_cube = sst_cube_full.extract(YEARS_between)
        
        # now decompose this cube into its components: spatial time-mean, spatial linear trend, 
        # mean seasonal cycle centred around year_centre
        monthly_mean = sst_cube.aggregated_by('month', iris.analysis.MEAN)
        print monthly_mean
        
        # remove mean seasonal cycle from all months
        yearly_monthly_trend = sst_cube_full.copy()
        for m in range(0,12):
            inds = np.where(sst_cube_full.coord('month').points == m+1)[0]
            for i in inds:
                yearly_monthly_trend.data[i] -= monthly_mean.data[m]
        iris.save(yearly_monthly_trend,file_trend[:-3]+'_allyrs_sst.nc', fill_value = 1.0e20)

        # add the model mask
        mask_path = form_BADC_path_mask1(runid)
        print mask_path
        if not os.path.exists(mask_path):
            print 'mask path not exist'
            mask_path = form_BADC_path_mask2(runid)
        mask = iris.load_cube(os.path.join(mask_path,'*.nc'),'sea_area_fraction')
        if 'IPSL' in run:
            # need to remove wrap rows from NEMO
            remove_nemo_wrap(os.path.join(mask_path,'*.nc'), os.path.join(savedir,'ipsl_mask.nc'))
            mask = iris.load_cube(os.path.join(savedir,'ipsl_mask.nc'))
            land = np.where(mask.data.mask == True)
            mask.data[...] = 1.0
            mask.data[land[0], land[1]] = 0.0
        print mask
        land = np.where(mask.data == 0.0)
        mask_land = yearly_monthly_trend[0].copy()
        add_mask(mask_land)
        add_mask(yearly_monthly_trend)
        print mask_land
        mask_land.data[...] = 1.0
        mask_land.data[land] = 0.0
        #mask_land.data[land] = np.ma.masked
        mask_land.data.mask[land] = True
        yearly_monthly_trend *= mask_land
        yearly_monthly_trend.data.mask[:,...] = mask_land.data.mask[...]

        cube_filled = fill_in_masked_data(yearly_monthly_trend)
        iris.save(cube_filled,file_trend[:-3]+'_allyrs_masked_sst.nc', fill_value = 1.0e20)

def plot_global_temp_changes(cube_dict, CMIP5_ref):
    fig = plt.figure(figsize=(10,6),dpi=100)#dpi=300
    subpl=fig.add_subplot(1,1,1)
    for run in CMIP5_ref:
    
        iplt.plot(cube_dict[run], linewidth=1.)
        plt.title("Timeseries of variability SST")
        #plt.savefig(os.path.join(savedir,'HighResMIP_ISST_v3_variability.png'))
    plt.show()
    

def calculate_global_temp_change(cube):
    grid_areas = guess_areas(cube)
    cube_ts = cube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
    return cube_ts
    
def guess_areas(cube):
    """
    Estimate area of grid cells by guessing bounds of 'latitude'
    and 'longitude' coordinates and return areas as an array.

    Args:

    * cube:
        iris.cube.Cube

    Retuns:
        Numpy array

    """
    try:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
    except:
        pass

    grid_areas = iris.analysis.cartography.area_weights(cube)

    return grid_areas

def remap_with_cdo(file_in, file_out, togrid = '025', rtype = 'bil', smooth = True):
    file_tmp = file_out[:-3]+'_tmp.nc'
    if rtype == 'bil':
        cdo_cmd = 'remapbil'
    else:
        cdo_cmd = 'remapycon'

    if togrid == '025':
        namelist_file = make_cdo_namelist_hadisst()
        cmd = 'cdo '+cdo_cmd+',r1440x720'+' '+file_in+' '+file_tmp
    elif togrid == '1':
        namelist_file = make_cdo_namelist_1x1()
        cmd = 'cdo '+cdo_cmd+',r360x180'+' '+file_in+' '+file_tmp
    elif togrid == '2':
        namelist_file = make_cdo_namelist_2x2()
        cmd = 'cdo '+cdo_cmd+',r180x90'+' '+file_in+' '+file_tmp
    elif togrid == 'hadisst025':
        namelist_file = make_cdo_namelist_hadisst()
        cmd = 'cdo '+cdo_cmd+','+namelist_file+' '+file_in+' '+file_tmp
    else:
        raise Exception('Dest grid for cdo regridding is unknown '+togrid)

#     cmd = 'cdo remapycon,'+namelist_file+' '+file_in+' '+file_tmp
    print ('regrid command ',cmd)
    os.system(cmd)

    if smooth:
        cmd = 'cdo smooth9 '+file_tmp+' '+file_out
    else:
        cmd = 'mv '+file_tmp+' '+file_out
    os.system(cmd)

    if os.path.exists(file_tmp):
        os.system('rm '+file_tmp)

def make_cdo_namelist_hadisst():
    namelist_file = os.path.join(os.getcwd(), 'namelist_grid')
    with open(namelist_file,'w') as outp:
        outp.write('gridtype = lonlat \n')
        outp.write('xsize = 1440 \n')
        outp.write('ysize = 720 \n')
        outp.write('xfirst = 0.12500 \n')
        outp.write('yfirst = -89.87500 \n')
        outp.write('yinc = 0.25000 \n')
        outp.write('\n')
    return namelist_file

def make_cdo_namelist_1x1():
    namelist_file = os.path.join(os.getcwd(), 'namelist_grid1x1')
    with open(namelist_file,'w') as outp:
        outp.write('gridtype = lonlat \n')
        outp.write('xsize = 360 \n')
        outp.write('ysize = 180 \n')
        outp.write('xfirst = 0.5 \n')
        outp.write('yfirst = -89.5 \n')
        outp.write('yinc = 1.0 \n')
        outp.write('\n')
    return namelist_file

def make_cdo_namelist_2x2():
    namelist_file = os.path.join(os.getcwd(), 'namelist_grid2x2')
    with open(namelist_file,'w') as outp:
        outp.write('gridtype = lonlat \n')
        outp.write('xsize = 180 \n')
        outp.write('ysize = 90 \n')
        outp.write('xfirst = 1.0 \n')
        outp.write('yfirst = -89. \n')
        outp.write('yinc = 2 \n')
        outp.write('\n')
    return namelist_file

def gridpoint_zonal_filter(cube, niter = 2):
    '''
    Filter the 2 grid waves on the model grid
    '''
    cube121 = cube.copy()
    for nit in range(niter):
        for i in range(1, cube.shape[-1]-1):
            cube121.data[:, :, i] = (cube121.data[:, :, i-1] + 2.0* cube121.data[:, :, i] + cube121.data[:, :, i+1]) / 4.0
        i = 0
        cube121.data[:, :, i] = (cube121.data[:, :, -1] + 2.0* cube121.data[:, :, i] + cube121.data[:, :, i+1]) / 4.0
        i = cube.shape[-1]-1
        cube121.data[:, :, i] = (cube121.data[:, :, i-1] + 2.0* cube121.data[:, :, i] + cube121.data[:, :, 0]) / 4.0
    return cube121


def remap_trends_to_common_grid(model, CMIP5_ref, do_plot = False):
    # remap each model trend field onto a 1x1 grid, so that they can be meansed together
    # make a simple lonlat grid for cdo to remap to
    global_temp_change = {}
    files_remapped = []
    #pickle_file = os.path.join(savedir, 'trends_1x1_allmodels.pkl')

    file_trend = os.path.join(savedir, model+'_month_annual_1950_2120_delta.nc')
    file_masked = file_trend[:-3]+'_allyrs_masked_sst.nc'
    file_masked_121 = file_trend[:-3]+'_allyrs_masked_sst_121.nc'

    if not os.path.exists(file_masked_121):
        c = iris.load_cube(file_masked)
        c_remove121 = gridpoint_zonal_filter(c, niter=3)
        iris.save(c_remove121, file_masked_121, fill_value = 1.0e20)

    #file_1x1 = file_masked[:-3]+'_1x1.nc'
    file_1x1 = file_masked_121[:-3]+'_1x1.nc'
    files_remapped.append(file_1x1)

    if not os.path.exists(file_1x1):
        #remap_with_cdo(file_masked, file_1x1)
        remap_with_cdo(file_masked_121, file_1x1, togrid = '1', rtype = 'bil', smooth = False)

    if do_plot:
        cube_ts_list = iris.load(file_1x1)
        trend_cube = cube_ts_list.concatenate_cube()
        guess_areas(trend_cube)
        icc.add_year(trend_cube, 'time', name = 'year')
        cube_ts = calculate_global_temp_change(trend_cube)
        global_temp_change[model] = cube_ts

    return files_remapped

def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[1:-1]

def time_filter_by_month(infile):
    '''
    Filter each month individually with a simple rolling window filter
    To reduce the interannual variability magnitude, since this is being provided by HadISST2.2
    Lanczos filter
    Note this removes the first and last n years (1/n filter)
    '''
    c = iris.load_cube(infile)
    icc.add_month_number(c, 'time', name = 'month')
    icc.add_year(c, 'time', name = 'year')

    cube_list_filtered = iris.cube.CubeList()

    # window length in years
    window = 11
    # Construct 5-year low pass filter
    wgts5 = low_pass_weights(window, 1. / 7.)
    fmonths = {}
    for mon in range(1,13):
        fmonths[mon] = []
    for mon in range(1,13):
        fout_month = infile[:-3]+'_filtered_'+str(mon).zfill(2)+'_7_11.nc'
        fmonths[mon].append(fout_month)
        if not os.path.exists(fout_month):
            print 'filter month ',mon
            con_mon = iris.Constraint(coord_values = {'month' : lambda l : l == mon})
            c_mon = c.extract(con_mon)
            c5 = c_mon.rolling_window('time',
                               iris.analysis.SUM,
                               len(wgts5),
                               weights=wgts5)
            cube_list_filtered.append(c5)
            iris.save(c5, fout_month, unlimited_dimensions = ['time'], fill_value = 1.0e20)
            fmonths[mon].append(fout_month)
        else:
            c5 = iris.load_cube(fout_month)

        year_start = int(c5.coord('year').points[0])
        year_end = int(c5.coord('year').points[-1])
    print 'years ',year_start, year_end
    years = range(year_start, year_end+1)

    fout_all_years_filtered = infile[:-3]+'_filtered_'+str(years[0])+'-'+str(years[-1])+'_7_11.nc'

    c_all_yr_l = iris.cube.CubeList()
    fnames = []
    for year in years:
        con_yr = iris.Constraint(coord_values = {'year' : lambda l : l == year})
        c_yr_l = iris.cube.CubeList()
        fyear = []
        for im in range(1,13):
            c = iris.load_cube(fmonths[im])
            c.coord('year').bounds = None
            c_yr = c.extract(con_yr)
            #c_yr.remove_coord('month')
            c_yr_l.append(c_yr)
        c_this_yr = c_yr_l.merge_cube()
        c_this_yr.remove_coord('month')
        c_this_yr.remove_coord('year')
        c_yr_l.append(c_this_yr)
        fout = infile[:-3]+'_filtered_'+str(year)+'_7_11.nc'
        fnames.append(fout)
        iris.save(c_this_yr, fout, unlimited_dimensions = ['time'], fill_value = 1.0e20)
    cmd = 'ncrcat '+' '.join(fnames)+' '+fout_all_years_filtered
    subprocess.call(cmd, shell = True)
#    c_all_yr = c_yr_l.concatenate_cube()
#    iris.save(c_all_yr, fout_all_years_filtered, unlimited_dimensions = ['time'], fill_value = 1.0e20)

    return fout_all_years_filtered

def monthly_running_mean(infile):
    '''
    Try removing mean seasonal cycle from all data
    Then perhaps filtering the rest (or else extracting the trend for each month
    '''
    file_runningmean = infile[:-3]+'_running_monthlymean.nc'
    if os.path.exists(file_runningmean):
        return file_runningmean
    c = iris.load_cube(infile)
    icc.add_month_number(c, 'time', name = 'month')
    icc.add_year(c, 'time', name = 'year')

    #cube_anomaly = sst_future.remove_monthly_mean_time_avg(c)
    cube_mean = sst_future.monthly_mean_running_time_avg(c)

    #iris.save(cube_anomaly, infile[:-3]+'_remove_monthlymean.nc', unlimited_dimensions = ['time'], fill_value = 1.0e20)
    iris.save(cube_mean, file_runningmean, unlimited_dimensions = ['time'], fill_value = 1.0e20)

    return file_runningmean
    
def calculate_mean_delta(CMIP5_ref):
    # calculate the mean of all the remapped trend files
    # should simply mean averaging over all cubes, since should have same dimensions at this point
    remapped_files = {}
    model_list = []
    model_trend_ending = '_month_annual_1950_2120_delta'
    for im, model in enumerate(CMIP5_ref):
        file_trend = os.path.join(savedir, model+model_trend_ending+'.nc')
        file_masked = file_trend[:-3]+'_allyrs_masked_sst_121.nc'
        file_1x1 = file_masked[:-3]+'_1x1.nc'
        remapped_files[model] = file_1x1
        if im == 0:
            cube_ref1 = iris.load_cube(file_1x1, callback = cmip5_callback)

        model_list.append(model.split('_')[0])

    model_subset = '_'.join(model_list)
    print remapped_files

    modelmean_file_out = os.path.join(savedir, 'cmip5_modelmean'+model_trend_ending+'_on_hadisst2_'+model_subset+'.nc')
    cube_list_new = iris.cube.CubeList()
    hadisst2_mask = iris.load_cube('/gws/nopw/j04/hrcm/cache/malcolm/HadISST2.2.2.0/hadisst_0to360_alldays_sst_1963.nc','sea_surface_temperature')[0]

    if not os.path.exists(modelmean_file_out):
        fnames = []
        cube_list_new = iris.cube.CubeList()
        for im, model in enumerate(CMIP5_ref):
            fname = os.path.join(savedir, 'cmip5_modelmean'+model_trend_ending+'_on_hadisst2_'+str(model)+'.nc')
            if not os.path.exists(fname):

                cube_model = iris.load_cube(remapped_files[model], callback = cmip5_callback)
                icc.add_year(cube_model, 'time', name = 'year')
                print 'cube_model ',cube_model
                cube_model = cube_model.extract(YEARS_between_1950_2100)
                cube_model.remove_coord('year')
            # this assumes all the files have the same period, and just a different calendar
            #fix_time_coord(cube_model, cube_ref1)
                model_coord = iris.coords.AuxCoord(im, long_name = 'Model', units = '1')
                cube_model.add_aux_coord(model_coord)
                iris.save(cube_model, fname)
                fix_fill_value_with_nco(fname)  
            fnames.append(fname)

            print 'fnames to read into list ',fname
            c = iris.load_cube(fname, callback = cmip5_callback)
            c.coord('time').bounds = None
            if im == 0:
                cube_ref1 = c.copy()
                c1 = c
            else:
                c1 = fix_time_coord(c, cube_ref1)
                iris.util.describe_diff(c, cube_ref1)
                print c1.coord('time')
            c1.coord('Model').points = im

            cube_list_new.append(c1)
            #print 'cube_list_new ',cube_list_new
        
        cube = cube_list_new.merge_cube()
        if not os.path.exists(modelmean_file_out):
            cm = cube.collapsed('Model', iris.analysis.MEAN)
            iris.save(cm, modelmean_file_out, unlimited_dimensions = ['time'], fill_value = 1.0e20)

        del cube

    return modelmean_file_out

def regrid_mean_to_025(modelmean_file_out):
    '''
    Remap the averaged multi-model mean to the 0.25x0.25 grid
    '''

    hadisst2_grid = iris.load_cube('/group_workspaces/jasmin2/primavera1/WP6/forcing/HadISST2_submit/v1.2/tos_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_20140101-20141231.nc','sea_surface_temperature')[0]

    modelmean_025_file_out = modelmean_file_out[:-3]+'_025.nc'
    modelmean_025_file_out_iris = modelmean_file_out[:-3]+'_025_iris.nc'
    file_to_remap = modelmean_file_out

    if not os.path.exists(modelmean_025_file_out):
        remap_with_cdo(file_to_remap, modelmean_025_file_out, rtype = 'bil', togrid = 'hadisst025')
    if not os.path.exists(modelmean_025_file_out_iris):
        c = iris.load_cube(file_to_remap)
        cregrid = c.regrid(hadisst2_grid, iris.analysis.Linear())
        iris.save(cregrid, modelmean_025_file_out_iris)

if __name__ == '__main__':
    
    CMIP5_ref = CMIP5_datasets()

    cube_1950_2100 = find_and_merge_historical_rcp85(CMIP5_ref)
    #cube_1950_2100 = find_and_merge_historical_rcp85(CMIP5_ref, var = 'sic')
    
    for model in CMIP5_ref:
        make_trend_files_for_models(model, CMIP5_ref)

        remapped_files = remap_trends_to_common_grid(model, CMIP5_ref)

    modelmean_file_out = calculate_mean_delta(CMIP5_ref)
    modelmean_filtered_out = monthly_running_mean(modelmean_file_out)

    regrid_mean_to_025(modelmean_filtered_out)
    
    
