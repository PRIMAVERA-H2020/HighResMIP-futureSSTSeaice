#!/bin/env python2.7
'''
NAME:
    sst_ice_maximum_hadisst2

DESCRIPTION:
    calculate a mask indicating where sea-ice exists (at any time, for given month)
    will use this mask to constrain future sea-ice to be within region defined

USAGE: 
    Execute as a script from command line

AUTHOR:
    Malcolm Roberts (hadom)

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
import sys, copy
import glob
import calendar
from scipy.interpolate import interp1d
import cPickle as pickle #cPickle is an optimised version of Pickle and is O(100) times faster

DATADIR = '/gws/nopw/j04/hrcm/cache/malcolm/HadISST2/1x1/processing_2018/'
savedir = os.path.join(DATADIR, 'future2')
#DATADIR_new = '/group_workspaces/jasmin2/primavera1/WP6/forcing/FutureSSTandSeaice/'
savedir_ice = os.path.join(DATADIR, 'future2')
YEARS = range(1950, 2016)
CODE_DIR = '/home/users/mjrobert/workspace/primavera-code-repository/HighResMIP/'

month_jan = iris.Constraint(coord_values = {'month_number' : lambda l : l == 1})
npole = iris.Constraint(coord_values = {'latitude' : lambda l : 50 < (l.point) <= 90})
spole = iris.Constraint(coord_values = {'latitude' : lambda l : -90 < (l.point) <= -50})
pole_title = ['Arctic, daily HadISST2, 1962-1989','Antarctic, daily HadISST2, 1962-1989']

# need to run from the output directory
# assumes that years in file names are real years
BIN = {}; BIN['SST'] = {}; BIN['SICE'] = {}
BIN['SST']['MAX_MIN'] = [-2.8, 5.0]
BIN['SST']['SIZE'] = 0.02
BIN['SST']['XBINS'] = np.arange(BIN['SST']['MAX_MIN'][0], BIN['SST']['MAX_MIN'][1], BIN['SST']['SIZE'])
BIN['SST']['NBINS'] = len(BIN['SST']['XBINS'])-1

SH_CONC_MIN = 50.0
FIRST_FUTURE_YEAR = 2016

CMIP5_ref = {'GFDL-CM3': 'NOAA-GFDL', 'IPSL-CM5A-MR': 'IPSL', 'CNRM-CM5': 'CNRM-CERFACS', 'HadGEM2-ES':'MOHC'}

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

def interpolate_histogram(year, mean_freq):
    '''
    Interpolate the monthly histogram to generate daily values, for the number
    of days in this year
    '''
    if calendar.isleap(year):
        ndays = 366
    else:
        ndays = 365

    midmonth = np.zeros((14))
    midmonth[0] = 1
    daynumber = 0
    for im in range(0,12):
        monthdays = calendar.monthrange(year, im+1)[1]
        if im > 0:
            monthadd = calendar.monthrange(year, im)[1]
        else:
            monthadd = 0
        daynumber += monthadd
        midmonth[im+1] = monthdays // 2 + daynumber
    midmonth[-1] = ndays
    alldays = np.arange(1, ndays+1)

    mean_freq_tmp = np.ma.zeros((14))
    mean_freq_days = np.ma.zeros((BIN['SST']['NBINS'], 2, ndays))
    
    print 'shapes ',mean_freq[year].shape, mean_freq_tmp.shape, midmonth.shape
    print midmonth

    for ir in range(0,2):
        for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
            mean_freq_tmp[1:13] = mean_freq[year][ib, ir, :]
            # set wrap points for year
            mean_freq_tmp[0] = (mean_freq[year][ib, ir, 0] + mean_freq[year][ib, ir, 11]) / 2.0
            mean_freq_tmp[13] = mean_freq_tmp[0]
            f = interp1d(midmonth, mean_freq_tmp, kind = 'linear')
            
            mean_freq_days[ib, ir, :] = f(alldays)

    return mean_freq_days

def read_histogram(year):
    yearstart = 2015; yearmeanstart = 2010
    mean_freq_future = np.ma.zeros((BIN['SST']['NBINS'], 2, 12, yearstart - yearmeanstart+1))
    mean_freq = {}; mean_freq_days = {}
    fnames = []

    if year < 2016:
        filename_pickle = os.path.join(DATADIR, 'sst_func_siconc_month_relationship_'+str(year)+'.pkl')
        if os.path.exists(filename_pickle):
            fh = open(filename_pickle, 'r')
            mean_freq_yr = pickle.load(fh)
            if year <= yearstart:
                mean_freq[year] = mean_freq_yr
            fh.close()
    else:
        for yr in range(yearmeanstart, yearstart+1):
            filename_pickle = os.path.join(DATADIR, 'sst_func_siconc_month_relationship_'+str(yr)+'.pkl')
            if os.path.exists(filename_pickle):
                fh = open(filename_pickle, 'r')
                mean_freq_yr = pickle.load(fh)
                mean_freq_future[:,:,:,yr-yearmeanstart] = mean_freq_yr
                #mean_freq[year] = mean_freq_yr
                fh.close()

        mean_freq_hadisst2 = np.mean(mean_freq_future, axis=3)

    # read CMIP5 histograms
        mean_freq_cmip5_all = np.ma.zeros((BIN['SST']['NBINS'], 2, 12, len(CMIP5_ref)))
        for im, model in enumerate(CMIP5_ref):
            filename_pickle = os.path.join(DATADIR, 'sst_func_siconc_month_relationship_'+model+'.pkl')
            if os.path.exists(filename_pickle):
                fh = open(filename_pickle, 'r')
                mean_freq_model = pickle.load(fh)
                fh.close()

            mean_freq_cmip5_all[:, :, :, im] = mean_freq_model
        mean_freq_cmip5 = np.mean(mean_freq_cmip5_all, axis = 3)

        # by 2030 we want to be using the CMIP5 relationship completely, at 2015 we want to use HadISST2.2
        year_start_linear = 2016
        year_end_linear = 2050
        if year < year_start_linear:
            mean_freq[year] = mean_freq_hadisst2
        elif year > year_end_linear:
            mean_freq[year] = mean_freq_cmip5
        else:
            weight = float(year - year_start_linear) / float(year_end_linear - year_start_linear)
            year_weight_hadisst = (1.0 - weight)
            year_weight_cmip5 = weight
            mean_freq[year] = year_weight_hadisst * mean_freq_hadisst2 + year_weight_cmip5 * mean_freq_cmip5
            
    mean_freq_days[year] = interpolate_histogram(year, mean_freq)

    return mean_freq_days

def ice_maximum_extent(years):
    '''
    Generate a 0/1 mask for any point which has ice concentration > 0 over time period
    Use this to mask future ice (make sure that you don't have ice where it never was before)
    '''
    dir_in = '/group_workspaces/jasmin2/primavera1/WP6/forcing/HadISST2_submit/v1.2/'
    year = 2015
    ice_hist_file = os.path.join(dir_in, 'siconc_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_'+str(year)+'0101-'+str(year)+'1231.nc')
    ice_ref = iris.load_cube(ice_hist_file)[0]

    ice_extent_filename = dir_in+'ice_max_month_'+str(years[0])+'-'+str(years[-1])+'.nc'
    ice_extent_min_filename = dir_in+'ice_min_month_'+str(years[0])+'-'+str(years[-1])+'.nc'

    if os.path.exists(ice_extent_filename) and os.path.exists(ice_extent_min_filename):
        ice_max_month = iris.load_cube(ice_extent_filename)
        ice_min_month = iris.load_cube(ice_extent_min_filename)
        return ice_max_month, ice_min_month, ice_ref

    ice_hist_l = iris.cube.CubeList()
    for year in years:
        ice_hist_file = os.path.join(dir_in, 'siconc_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_'+str(year)+'0101-'+str(year)+'1231.nc')
        ice = iris.load_cube(ice_hist_file, callback = history_period_callback)
        try:
            icc.add_month_number(ice, 'time', name = 'month')
        except:
            pass
        ice_hist_l.append(ice)
    ice_hist = ice_hist_l.concatenate_cube()

    cube_mask = ice_hist[:12].copy()
    cube_mask_min = ice_hist[:12].copy()
    for im in range(0,12):
        mon = im+1
        print 'process month ',mon
        month_constraint = iris.Constraint(month = mon)

        # calculate maximum ice conc over time period
        ice = ice_hist.extract(month_constraint)
        ice_max = ice.collapsed('time', iris.analysis.MAX)
        ice_max.remove_coord('month')
        ice_max.coord('time').bounds = None

        # calculate minimum ice conc over time period
        ice_min = ice.collapsed('time', iris.analysis.MIN)
        ice_min.remove_coord('month')
        ice_min.coord('time').bounds = None

        # calculate mask of 1 if any point ever has some ice
        cube_mask.data[im, :, :] = 0.0
        some_ice = np.where((ice_max.data > 1.0) & (ice_max.data <= 100.0))
        cube_mask.data[im, some_ice[0], some_ice[1]] = 1.0
        print 'cube_mask ',cube_mask[im].coord('time')

        # calculate mask of 1 if any point always has some ice
        cube_mask_min.data[im, :, :] = 0.0
        some_ice = np.where((ice_min.data > 80.0) & (ice_min.data <= 100.0))
        cube_mask_min.data[im, some_ice[0], some_ice[1]] = 1.0
        print 'cube_mask ',cube_mask_min[im].coord('time')

    iris.save(cube_mask, ice_extent_filename)
    iris.save(cube_mask_min, ice_extent_min_filename)
            
    return cube_mask, cube_mask_min, ice_ref

def add_metadata(cube):
    '''
    Change metadata from tos to siconc
    '''
    cube.var_name = 'siconc'
    cube.long_name = 'HighResMIP Future Sea Ice Concentration'
    cube.standard_name = 'sea_ice_area_fraction'

def get_midmonth_icevalues(sst_this_month, sst_last_month, sst_next_month, mean_freq, year, im, lat_shape, siconc_temp, mm1, mp1):
    '''
    Get the mid-month sea-ice values from the current and surrounding months
'''
    ndays = sst_this_month.shape[0]
    ndays_last = sst_last_month.shape[0]
    ndays_next = sst_next_month.shape[0]

    mid_month = ndays // 2
    mid_month_last = ndays_last // 2
    mid_month_next = ndays_next // 2

    siconc_this = siconc_temp.copy()
    siconc_last = siconc_temp.copy()
    siconc_next = siconc_temp.copy()

    for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
        sst_range = np.where((bin < sst_this_month.data[mid_month, :, :]) & (bin+1 >= sst_this_month.data[mid_month, :, :]) & (ice_mask.data[im, :, :] == 1))
        nhemi = np.where(sst_range[0] > lat_shape / 2)
        shemi = np.where(sst_range[0] < lat_shape / 2)
        siconc_this.data[sst_range[0][nhemi], sst_range[1][nhemi]] = mean_freq[year][ib, 0, im]
        siconc_this.data[sst_range[0][shemi], sst_range[1][shemi]] = mean_freq[year][ib, 1, im]

    for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
        sst_range = np.where((bin < sst_last_month.data[mid_month_last, :, :]) & (bin+1 >= sst_last_month.data[mid_month_last, :, :]) & (ice_mask.data[mm1, :, :] == 1))
        nhemi = np.where(sst_range[0] > lat_shape / 2)
        shemi = np.where(sst_range[0] < lat_shape / 2)
        siconc_last.data[sst_range[0][nhemi], sst_range[1][nhemi]] = mean_freq[year][ib, 0, mm1]
        siconc_last.data[sst_range[0][shemi], sst_range[1][shemi]] = mean_freq[year][ib, 1, mm1]

    for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
        sst_range = np.where((bin < sst_next_month.data[mid_month_next, :, :]) & (bin+1 >= sst_next_month.data[mid_month_next, :, :]) & (ice_mask.data[mp1, :, :] == 1))
        nhemi = np.where(sst_range[0] > lat_shape / 2)
        shemi = np.where(sst_range[0] < lat_shape / 2)
        siconc_next.data[sst_range[0][nhemi], sst_range[1][nhemi]] = mean_freq[year][ib, 0, mp1]
        siconc_next.data[sst_range[0][shemi], sst_range[1][shemi]] = mean_freq[year][ib, 1, mp1]

    return siconc_this, siconc_last, siconc_next

def merge_first_year_to_historic(future_fname):
    '''
    Read last year of real sea-ice data, and first year of new data
    Merge them over the first month
    '''
    dir_in_hist = '/group_workspaces/jasmin2/primavera1/WP6/forcing/HadISST2_submit/v1.2/'
    fice = 'siconc_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_'+str(FIRST_FUTURE_YEAR-1)+'0101-'+str(FIRST_FUTURE_YEAR-1)+'1231.nc'

    sic = iris.load_cube(dir_in_hist+fice)
    icc.add_month_number(sic, 'time', name = 'month')
    month_constraint_dec = iris.Constraint(coord_values = {'month': lambda l : l == 12})
    month_constraint_jan = iris.Constraint(coord_values = {'month': lambda l : l == 1})
    sic_hist = sic.extract(month_constraint_dec)

    sic_future = iris.load_cube(future_fname)
    sic_future_month = sic_future.extract(month_constraint_jan)
    ndays = sic_future_month.shape[0]

    for iday in range(0, ndays):
        weight = float(iday)/float(ndays-1)
        sic_future.data[iday, :, :] = sic_hist.data[-1, :, :] * (1.0 - weight) + sic_future_month.data[iday, :, :] * weight

    iris.save(sic_future, future_fname[:-3]+'_merged.nc')
    cmd = 'mv '+future_fname+' '+future_fname[:-3]+'_unmerged.nc'
    os.system(cmd)
    cmd = 'mv '+future_fname[:-3]+'_merged.nc'+' '+future_fname
    os.system(cmd)

def generate_future_siconc_from_sst(mean_freq_days, year, months, ice_mask, ice_mask_min, ice_ref):
    '''
    Use the maximum ice extent, the SST and the pdf relationship between SST and sea-ice to generate a future sea ice concentration
    '''
    #dir_in = '/group_workspaces/jasmin2/primavera1/WP6/forcing/HadISST2_submit/v1.2/'
    if year >= 2016:
        #future_sst_file = os.path.join(DATADIR, 'future', 'full_sst_1950_2100_025_daily_fixed.nc')
        future_sst_file = os.path.join(savedir, 'sst', 'future_sst_'+str(year)+'_025_daily_v1.nc')
    else:
        if len(months) == 12:
            future_sst_file = os.path.join(DATADIR, 'hadisst2_tos_daily_'+str(year)+'_fixed_day_v1_monthly_under_ice.nc')
        elif len(months) == 1:
            future_sst_file = os.path.join(DATADIR, 'hadisst2_tos_daily_'+str(year)+'_fixed_'+str(months[0]+1).zfill(2)+'.nc')
        else:
            raise Exception('Must be either 12 months or 1 '+str(len(months)))

    future_siconc_file = os.path.join(savedir_ice, 'siconc', 'future_siconc_{}_025_daily_v1.nc')
    print 'read sst ',future_sst_file
    full_sst = iris.load_cube(future_sst_file)

    try:
        #icc.add_month_number(full_siconc, 'time', name = 'month')
        icc.add_month_number(full_sst, 'time', name = 'month')
    except:
        pass

    try:
        #icc.add_year(full_siconc, 'time', name = 'year')
        icc.add_year(full_sst, 'time', name = 'year')
    except:
        pass

    #print 'full cube ',full_siconc.coord('year')
    new_sice_year = iris.cube.CubeList()

    print 'YEAR ',year
    fout_year = future_siconc_file.format(str(year))
    new_sice_month = iris.cube.CubeList()
    year_con = iris.Constraint(coord_values = {'year' : lambda l : l == int(year)})
    sst_this_year = full_sst.extract(year_con)
    full_siconc = sst_this_year.copy()
    print 'cube year ',sst_this_year

    lat_shape = sst_this_year.shape[1]
    for im in months:
        if im > 0 and im < 11:
            mm1 = im - 1 
            mp1 = im + 1
        elif im == 0:
            mm1 = 11
            mp1 = im + 1
        elif im == 11:
            mm1 = im - 1
            mp1 = 0

        month = im+1
        month_con = iris.Constraint(coord_values = {'month' : lambda l : l == int(month)})
        sst_this_month = sst_this_year.extract(month_con)

        if im == 0:
            month_con_m1 = iris.Constraint(coord_values = {'month' : lambda l : l == 1})
        else:
            month_con_m1 = iris.Constraint(coord_values = {'month' : lambda l : l == int(month-1)})
        sst_last_month = sst_this_year.extract(month_con_m1)

        if im ==11:
            month_con_p1 = iris.Constraint(coord_values = {'month' : lambda l : l == 11})
        else:
            month_con_p1 = iris.Constraint(coord_values = {'month' : lambda l : l == int(month+1)})
        sst_next_month = sst_this_year.extract(month_con_p1)


        siconc_template = full_siconc.extract(month_con)
        coord_names = [coord.name() for coord in siconc_template.coords()]
        print coord_names
        if 'day_of_year' not in coord_names:
            icc.add_day_of_year(siconc_template, 'time', name = 'day_of_year')
        siconc_template.data[...] = 0.0

        ndays = sst_this_month.shape[0]
        day_of_year = siconc_template.coord('day_of_year').points

        for day in range(ndays):
            day_number = day_of_year[day] - 1
            print 'process ',day, month, year, day_number
            for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
                sst_range = np.where((bin < sst_this_month.data[day, :, :]) & (bin+1 >= sst_this_month.data[day, :, :]) & (ice_mask.data[im, :, :] == 1))
                nhemi = np.where(sst_range[0] > lat_shape / 2)
                shemi = np.where(sst_range[0] < lat_shape / 2)
                siconc_template.data[day, sst_range[0][nhemi], sst_range[1][nhemi]] = mean_freq_days[year][ib, 0, day_number]
                siconc_template.data[day, sst_range[0][shemi], sst_range[1][shemi]] = mean_freq_days[year][ib, 1, day_number]
            print 'completed ',day, month, year

        siconc_template.data.mask[:, sst_this_month[0].data.mask] = True
        siconc_template.data.mask[:, ice_ref.data.mask] = True
        new_sice_month.append(siconc_template)
    new_sice_cube = new_sice_month.concatenate_cube()
    add_metadata(new_sice_cube)
    print 'write to ',fout_year
    iris.save(new_sice_cube, fout_year, unlimited_dimensions = ['time'], fill_value = 1.0e20)

    return fout_year

if __name__ == '__main__':
    months = range(0,12)
    # period to average HadISST sic = f(sst)
    years_past = range(2000, 2016)
    ice_mask, ice_mask_min, ice_ref = ice_maximum_extent(years_past)

    # period to create future seaice files
    years_hist = range(2016,2055)
    #years_hist = range(2016,2017)
    for year in years_hist:
        mean_freq = read_histogram(year)
        future_siconc_file = generate_future_siconc_from_sst(mean_freq, year, months, ice_mask, ice_mask_min, ice_ref)

        if year == FIRST_FUTURE_YEAR:
            merge_first_year_to_historic(future_siconc_file)


        
    
