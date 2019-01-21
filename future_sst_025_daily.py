#!/usr/local/sci/bin/python2.7
'''
NAME:
    make_sst-future

Description:
    Produce the future SST dataset for the HighResMIP experiment highresSST-future

Method:

    Overview:
        Want to use the HadISST2 variability (monthly to interannual) together with the CMIP5 SST projections (derived from the RCP8.5 - historic) to produce a continuous, 1948-2100, (daily) SST dataset
        Need to calculate a detrended HadISST dataset, then use this to add to the trend from the CMIP5 models - replacing the longer term variability

    Details:
        Read HadISST2 data for full period - sst_cube_full
        Extract years between 1965 and 2009 - sst_cube
        Convert to SST if necessary
        Calculate annual mean for each year - sst_cube_annual
        Decompose the cube into temporal components: spatial time-mean, spatial linear trend, spatial variability
            Calculate monthly mean - monthly_mean
            Remove monthly daily mean - sst_nomonth (using sst_cube) - not sure exactly what this does
            
            For each month, 
                Extract data for month between year limits
                Calculate trend for this month at each point (hadisst2.fit_trend_per_point)
                Put these monthly trends into observed_trend
            

        Start from the HadISST2 dataset (monthly 1degree)
        Decompose this into time components

    For each year (every 12 points) nt
        For each month m
            Residual SST variability (sst_variability[nt+, :]) = full SST (sst_cube[nt+m, :]) - monthly_mean[m,:] - observed_trend[m,:]*(nt/12)

    
    Copy SST cube for future: future_sst = sst_cube
    Load cmip5_trend

    For each year (every 12 points) nt of future_sst
        Calculate fraction of years into future we are

        For each month m
            future_sst[nt+m, :] = monthly_mean[m, :] 
                                + cmip5_trend[m,:]*years_into_future_p1
                                + sst_variability[nt+m, :]
            # this one is the second time interval future (so offset by nyears of the period chosen)
            future2_sst[nt+m, :] = monthly_mean[m, :] 
                                + cmip5_trend[m,:]*years_into_future_p2
                                + sst_variability[nt+m, :]

        (Needs thinking about for daily data here)
    
    Make a cube list
    Append sst_cube
    Append future_sst
    Append future_sst2

    Various assumptions about monthly data
    Not correct years
    Only able to add CMIP5 because monthly
    Need to enable scale up to daily for some accumulate routines
    Note that around ice edge the trends are huge

LAST MODIFIED:
    2017-03-10

'''

import numpy as np
import os, sys
import iris
import iris.analysis
import iris.coord_categorisation as icc
import cf_units, datetime

import sst_future as sf

DATADIR_orig = '/group_workspaces/jasmin2/primavera1/WP6/forcing/HadISST2_submit/v1.2/'
DATADIR = '/gws/nopw/j04/hrcm/cache/malcolm/HadISST2/1x1/processing_2018'  # from 1962-2009, centred on 2014
CMIP5_trend_file = '/gws/nopw/j04/hrcm/cache/malcolm/HighResMIP/sst_forcing/processing_ocean_v3/cmip5_modelmean_month_annual_1950_2120_delta_on_hadisst2_IPSL-CM5A-LR_HadGEM2-ES_GFDL-CM3_ACCESS1-0_ACCESS1-3_MPI-ESM-MR_CNRM-CM5_IPSL-CM5A-MR_running_monthlymean_025_iris.nc'

SST_MIN = -2.0

SAVEDIR = DATADIR+'/future2/'
if not os.path.exists(SAVEDIR): 
    os.makedirs(SAVEDIR)
if not os.path.exists(SAVEDIR+'sst'): 
    os.makedirs(SAVEDIR+'sst')
if not os.path.exists(SAVEDIR+'siconc'): 
    os.makedirs(SAVEDIR+'siconc')

def callback_attrib(cube, field, filename):
    # add a categorical year coordinate to the cube
    attributes = ['time_period', 'history']
    for attrib in attributes:
        try:
            del cube.attributes[attrib]
        except:
            pass


def add_time_names(cube, year=False, month_number=False, day_of_year=False):
    '''
    Add additional coordiate names to time dimension
    '''
    if year:
        icc.add_year(cube, 'time', name = 'year')

    if month_number:
        try:
            icc.add_month_number(cube, 'time', name = 'month')
        except:
            pass

    if day_of_year:
        try:
            icc.add_day_of_year(cube, 'time', name = 'day_of_year')
        except:
            pass

def extract_cube_years(cube, start_yr, end_yr):
    '''
    Extract years in cube between start and end
    '''
    years_between = iris.Constraint(coord_values = {'year' : lambda l : int(start_yr) <= l <= int(end_yr)})

    cube_subset = cube.extract(years_between)
    return cube_subset

def aggregate_mean(cube, coord = 'year'):
    '''
    Calculate annual means
    '''
    cube_agg_by_coord = cube.aggregated_by(coord, iris.analysis.MEAN)
    print 'sst_cube_year ',cube_agg_by_coord
    return cube_agg_by_coord

def calc_observed_trend(cube, sst_remove_monthmean):
    observed_trend = cube[-12:].copy()
    for m in range(0,11+1):
        print cube.coord('month').points[m],cube.coord('month').points[-12+m]
        print cube.coord('time').points[m], cube.coord('time').points[-12+m]
        time_gap = cube.coord('time').points[-12] - cube.coord('time').points[0]
        print 'time_gap ',time_gap
        month_constraint = iris.Constraint(month = m+1)
        month_extract = sst_remove_monthmean.extract(month_constraint)
        print 'month extract for obs trend ',month_extract
        month_trend = sf.fit_trend_per_point(month_extract)
        # initial calculation obtained bad values along Antarctica edge - maybe the missing data comes and goes? Check for bad values and reset to zero
        miss = np.where(np.abs(month_trend.data) >= 1.0)
        if len(miss) > 0:
            print 'dodgy trend at ',miss,month_trend.data[miss[0], miss[1]]
            month_trend.data[miss[0], miss[1]] = 0.0
        observed_trend.data[m,:] = month_trend.data[:]
    return observed_trend

def calc_variability(cube, monthly_mean, year, years_from_start, save = True, weighted = True):
    '''
    Currently assumes cube using monthly mean data
    For daily, input cube has dimensions:
        time, latitude, longitude
        where time is years, months, days
    '''
    # Calculate the observed SST variability with monthly mean and trend removed
    #sst_variability = cube.copy()
    #print 'sst_variability ',sst_variability
    
    hadisst2_mask = iris.load_cube('/gws/nopw/j04/hrcm/cache/malcolm/HadISST2.2.2.0/hadisst_0to360_alldays_sst_1963.nc','sea_surface_temperature')[0]

    sst_var_yr_l = iris.cube.CubeList()
    sst_var_yr_day_l = iris.cube.CubeList()
    print 'calc_var ',monthly_mean.coord('year')

    for m in range(0, 12):
        print 'process calc_var , m = ',m
        con_mon = iris.Constraint(coord_values = {'month' : lambda l : l == m+1})
        con_monm1 = iris.Constraint(coord_values = {'month' : lambda l : l == m})
        con_monp1 = iris.Constraint(coord_values = {'month' : lambda l : l == m+2})
        con_yr = iris.Constraint(coord_values = {'year' : lambda l : l == year})
        con_yrm1 = iris.Constraint(coord_values = {'year' : lambda l : l == year})
        con_yrp1 = iris.Constraint(coord_values = {'year' : lambda l : l == year})
        # this month
        cube_mon = cube.extract(con_yr & con_mon) # daily data
        cube_mon_mn = cube_mon.collapsed('time', iris.analysis.MEAN)

        cube_mon_m1 = cube.extract(con_yrm1 & con_monm1) # daily data
        cube_mon_p1 = cube.extract(con_yrp1 & con_monp1) # daily data

        # long term monthly mean
        monthly_mean_m1 = monthly_mean.extract(con_monm1&con_yr)
        monthly_mean_c = monthly_mean.extract(con_mon&con_yr)
        monthly_mean_p1 = monthly_mean.extract(con_monp1&con_yr)

        if m == 0:
            # previous month
            con_yrm1 = iris.Constraint(coord_values = {'year' : lambda l : l == year - 1})
            con_monm1 = iris.Constraint(coord_values = {'month' : lambda l : l == 12})
            cube_mon_m1 = cube.extract(con_yrm1 & con_monm1) # daily data

            # long term monthly mean
            monthly_mean_m1 = monthly_mean.extract(con_monm1&con_yrm1)
            if monthly_mean_m1 == None:
                print 'SHould be using the previous Dec data, but not available ',m, year
                monthly_mean_m1 = monthly_mean.extract(con_monm1&con_yr)

        elif m == 11:
            # next month
            con_yrp1 = iris.Constraint(coord_values = {'year' : lambda l : l == year + 1})
            con_monp1 = iris.Constraint(coord_values = {'month' : lambda l : l == 1})
            cube_mon_p1 = cube.extract(con_yrp1 & con_monp1) # daily data

            # long term monthly mean
            monthly_mean_p1 = monthly_mean.extract(con_monp1&con_yrp1)

        else:
            pass

        # calculate monthly means for this local month
        try:
            cube_mon_m1_mn =  cube_mon_m1.collapsed('time', iris.analysis.MEAN)
        except:
            print 'var, no cube for year, month m=0 ', year, year-1, m
            cube_mon_m1_mn = cube_mon_mn

        try:
            cube_mon_p1_mn =  cube_mon_p1.collapsed('time', iris.analysis.MEAN)
        except:
            print 'var, no cube for year, month m=11 ', year, year+1, m
            cube_mon_p1_mn = cube_mon_mn

        # extract this month's daily data
        cube_mon = cube.extract(con_yr & con_mon) # daily data
        cube_var_day = cube_mon.copy()
        cube_var_mon = cube_mon.copy()
        ndays = cube_mon.shape[0]

        # check units
        #print 'units ',monthly_mean.units, cube_mon.units
        #if cube_mon.units != monthly_mean.units:
        #    monthly_mean.units = cube_mon.units

        for iday in range(0, ndays):
            weight = float(iday)/float(ndays-1)
            w1 = (1.0 - weight) * 0.5
            w2 = 0.5
            w3 = weight * 0.5
            if weighted:
                mon_day_long = w1 * monthly_mean_m1.data[:, :] + w2 * monthly_mean_c.data[:, :] + w3 * monthly_mean_p1.data[:, :]
                mon_day = w1 * cube_mon_m1_mn.data[:, :] + w2 * cube_mon_mn.data[:, :] + w3 * cube_mon_p1_mn.data[:, :]
            else:
                mon_day_long = monthly_mean_c.data[:, :]
                mon_day = w1 * cube_mon_m1_mn.data[:, :] + w2 * cube_mon_mn.data[:, :] + w3 * cube_mon_p1_mn.data[:, :]

            #cube_var_mon.data[iday, :, :] = (cube_mon.data[iday, :, :] - cube_mon_mn.data[:,:]) \
            #                     - (cube_mon_mn.data[:,:] - mon_day)
            cube_var_day.data[iday, :, :] = (cube_mon.data[iday, :, :] - mon_day)
            cube_var_mon.data[iday, :, :] = (mon_day - mon_day_long)

                                 #- observed_trend.data[m, :, :]*float(years_from_start)
        cube_var_day.data.mask = hadisst2_mask.data.mask
        cube_var_mon.data.mask = hadisst2_mask.data.mask
        sst_var_yr_day_l.append(cube_var_day)
        sst_var_yr_l.append(cube_var_mon)

        print 'end calc_var month loop, ',m, sst_var_yr_l
    sst_var_yr_day = sst_var_yr_day_l.concatenate_cube()
    sst_var_yr = sst_var_yr_l.concatenate_cube()

    if save:
        savefile = os.path.join(SAVEDIR, 'sst', 'sst_variability_'+str(year)+'_025_daily_fixed_v1.nc')
        iris.save(sst_var_yr_day, savefile, unlimited_dimensions = ['time'], fill_value = 1.0e20)
        savefile = os.path.join(SAVEDIR, 'sst', 'sst_variability_'+str(year)+'_025_month_fixed_v1.nc')
        iris.save(sst_var_yr, savefile, unlimited_dimensions = ['time'], fill_value = 1.0e20)
    else:
        savefile = os.path.join(SAVEDIR, 'sst', 'sst_temp_'+str(year)+'_025_daily_fixed_v1.nc')
        iris.save(sst_var_yr, savefile, unlimited_dimensions = ['time'], fill_value = 1.0e20)

                                 
    return sst_var_yr, sst_var_yr_day

def add_future_time_info(start_cube):
    '''
    Fix calendar so that future_cube has dates starting from start_cube+1
    '''
    future_cube = start_cube.copy()
    
    for coord in ['year','month']:
        try:
            future_cube.remove_coord(coord)
        except:
            pass

    print 'start future year ',future_cube.coord('time').points[0:50]
    print 'start future year ',future_cube.coord('time').points[-50:]
    # transform time dimension to make this go from Jan onwards
    start_period2 =  start_cube.coord('year').points[-1]+1
    start_period1 =  start_cube.coord('year').points[0]
    print 'start period 1,2 ',start_period1, start_period2

    # find difference between date wanted and date have
    time_delta = datetime.date(start_period2,1,1) - datetime.date(start_period1,1,1)
    print 'time delta ',time_delta
    time_units = start_cube.coord('time').units
    print 'time_units ',time_units
    if 'hours' in str(time_units):
        time_delta_unit = 24.
    else:
        time_delta_unit = 1.

    # add this delta on to time coordinate
    future_cube.coord('time').points = future_cube.coord('time').points + time_delta.days*time_delta_unit
    
    icc.add_year(future_cube, 'time', name = 'year')
    icc.add_month_number(future_cube, 'time', name = 'month')

    print 'end future year ',future_cube.coord('time').points[0:50]
    print 'end future year ',future_cube.coord('year').points[0:50]

    return future_cube

def add_metadata(cube):
    '''
    Change metadata from tos to siconc
    '''
    cube.var_name = 'tos'
    cube.long_name = 'HighResMIP Future Sea Surface Temperature'
    cube.standard_name = 'sea_surface_temperature'

def get_cmip5_this_year(yr, year_hist, cmip5_trend, year_constraint, year_constraint_hist):
    '''
    check if this year is within cmip5 dataset, 
    if not use the last existing year of real data
    '''
    year_constraint_m1 = iris.Constraint(coord_values = {'year' : lambda l : l == yr-1})
    year_constraint_p1 = iris.Constraint(coord_values = {'year' : lambda l : l == yr+1})

    cmip5_yr = cmip5_trend.extract(year_constraint)
    cmip5_yrm1 = cmip5_trend.extract(year_constraint_m1)
    cmip5_yrp1 = cmip5_trend.extract(year_constraint_p1)

    print 'cmip5_yr ',cmip5_yr
    if cmip5_yr.shape[0] == 12 and cmip5_yrm1.shape[0] == 12 and cmip5_yrp1.shape[0] > 1:
        print 'cmip5_yr ',cmip5_yr
    else:
        raise Exception('Not 12 months in this year, or last year, or next year is empty '+yr)

    return cmip5_yr, cmip5_yrm1, cmip5_yrp1

def get_daily_interpolated_data(m, sst_variability, sst_variability_mon, future_c, monthly_mean_recent, cmip5_yr, cmip5_yrm1, cmip5_yrp1, year_constraint, year_constraint_hist, monthly_mean_previous, monthly_mean_next, year_last_decade, sst_variability_yrm1 = [], sst_variability_mon_yrm1 = [], adjust_jan = False):
    '''
    Calculate daily interpolated values of SST
    Have 12 months of monthly_mean_recent files
    Have year -1 , year and year + 1 of cmip5_yr
    Assuming that the months are in the correct 1-12 order (used by indexing). 
    '''
    con_mon = iris.Constraint(coord_values = {'month' : lambda l : l == m+1})
    con_yr_last = iris.Constraint(coord_values = {'year' : lambda l : l == year_last_decade})
    con_yr_lastm1 = iris.Constraint(coord_values = {'year' : lambda l : l == year_last_decade-1})
    con_yr_lastp1 = iris.Constraint(coord_values = {'year' : lambda l : l == year_last_decade+1})
    future_yrmon = future_c.extract(year_constraint & con_mon)
    sst_var = sst_variability.extract(year_constraint_hist & con_mon)
    sst_var_mon = sst_variability_mon.extract(year_constraint_hist & con_mon)
    print 'month, year ',m, year_last_decade
    print 'sst_var ',sst_var
    print 'future_yrmon ',future_yrmon
    print 'cmip5_yr ',cmip5_yr
    print 'monthly_mean ',monthly_mean_recent
    ndays = sst_var.shape[0]
    monthly_mean_c = monthly_mean_recent.extract(con_yr_last)
    if m > 0 and m < 11:
        cmip5_yr_monm1 = cmip5_yr
        cmip5_yr_monp1 = cmip5_yr
        mm1 = m - 1 
        mp1 = m + 1
        monthly_mean_mm1 = monthly_mean_recent.extract(con_yr_last)
        monthly_mean_mp1 = monthly_mean_recent.extract(con_yr_last)
    elif m == 0:
        cmip5_yr_monm1 = cmip5_yrm1
        cmip5_yr_monp1 = cmip5_yr
        mm1 = 11
        mp1 = m + 1
        if adjust_jan:
            monthly_mean_mm1 = monthly_mean_previous
        else:
            monthly_mean_mm1 = monthly_mean_recent.extract(con_yr_lastm1)
        monthly_mean_mp1 = monthly_mean_recent.extract(con_yr_last)
    elif m == 11:
        cmip5_yr_monm1 = cmip5_yr
        cmip5_yr_monp1 = cmip5_yrp1
        mm1 = m - 1
        mp1 = 0
        monthly_mean_mm1 = monthly_mean_recent.extract(con_yr_last)
        monthly_mean_mp1 = monthly_mean_next.extract(con_yr_lastp1)

    print 'month ',m
    print 'monthly_mean_c ', monthly_mean_c
    print 'monthly_mean_mm1 ', monthly_mean_mm1
    print 'monthly_mean_mp1 ', monthly_mean_mp1

    for iday in range(0, ndays):
        weight = float(iday)/float(ndays-1)
        w1 = (1.0 - weight) * 0.5
        w2 = 0.5
        w3 = weight * 0.5
        if cmip5_yr_monm1.coord('month').points[mm1] != mm1 + 1:
            raise Exception('CMIP5 month data not the right month '+str(cmip5_yr_monm1.coord('month').points[mm1])+', '+str(mm1+1))
        cmip5_day = w1 * cmip5_yr_monm1.data[mm1, :, :] + w2 * cmip5_yr.data[m, :, :] + w3 * cmip5_yr_monp1.data[mp1, :, :]
        if monthly_mean_mm1.coord('month').points[mm1] != mm1+1:
            raise Exception('monthly mean mm1 month not right month ')
        if monthly_mean_c.coord('month').points[m] != m+1:
            raise Exception('monthly mean mm1 month not right month ')
        mon_day = w1 * monthly_mean_mm1.data[mm1, :, :] + w2 * monthly_mean_c.data[m, :, :] + w3 * monthly_mean_mp1.data[mp1, :, :]

        # adjust the start of the merged year from real data
        nmonths_adjust = 1
        if adjust_jan and m <= nmonths_adjust-1:
            m_5_days = [31, 29, 31, 30, 31] # Jan-Apr
            mon_day_weight = float((sum(m_5_days[:m]) + iday)) / float(sum(m_5_days[:nmonths_adjust]))  # 0-1
            print 'month, day, mon_day_weight ', m, iday, mon_day_weight
            # merge the monthly anomaly over 5 months to gradually bring it into the first year of new data from the old obs data
            # needs to start at exactly end 2015, and end at the end-of-month 5 (m=4) monthly anomaly
            mon_merge = sst_variability_mon_yrm1.data[-1, :, :] * (1.0 - mon_day_weight)
            mon_var = (sst_var_mon.data[iday, :, :] + mon_day) * mon_day_weight
            mon_edge = monthly_mean_previous.data[11, :, :] * (1.0 - mon_day_weight)

            if m == 0:
            # use the last day of the previous year's exact variability, and average with this year's variability over the first month to merge smoothly
                #mon_merge = sst_variability_mon_yrm1.data[-1, :, :] * (1.0 - weight) + \
                #            (sst_var_mon.data[iday, :, :]) * weight
                var_merge = sst_variability_yrm1.data[-1, :, :] * (1.0 - weight) +  sst_var.data[iday, :, :] * weight

            else:
                var_merge = sst_var.data[iday, :, :]
            # then need to stop a sudden jump to the next month, so continue to merge over a few months
            future_yrmon.data[iday, :, :] = cmip5_day * mon_day_weight + mon_merge + mon_edge + mon_var + var_merge

        else:
            future_yrmon.data[iday, :, :] = cmip5_day  + mon_day + \
                                    sst_var.data[iday, :, :] + sst_var_mon.data[iday, :, :] 

    # reset points below SST_MIN
    too_cold = np.where((future_yrmon.data[:, :, :] < SST_MIN) & (future_yrmon.data.mask[:, :, :] == False))
    future_yrmon.data[too_cold[0], too_cold[1], too_cold[2]] = SST_MIN

        
    return future_yrmon

def adjust_monthly_mean_edge(cube_monthly_mean_recent, fname, adjust_dec = False, adjust_jan = False):
    '''
    Need to adjust Jan value of monthly_mean_recent to smoothly evolve from fname December
    '''
    cube = iris.load_cube(fname)
    cube.convert_units('degC')
    
    if adjust_dec:
        # need to change the Dec value to be that from the previous period
        monthly_mean = aggregate_mean(cube, coord = 'month')
        cube_previous = monthly_mean
        icc.add_year(cube_previous, 'time', name = 'year')
        return cube_previous

    if adjust_jan:
        # need to change the jan value to be that from the next period
        monthly_mean = aggregate_mean(cube, coord = 'month')
        cube_next = monthly_mean
        icc.add_year(cube_next, 'time', name = 'year')
        return cube_next

def calc_future_sst(cube, observed_trend_file, running_mean_file, datafiles):
    '''
    Currently each period (1950-2015, 2016-2081, 2082-) has a join at the end
    when the start of the repeat period + CMIP trend is used.
    Need to figure out how to smooth this join, taking into account e.g. peak
    ENSO in Dec/Jan needs to be reduced slowly

    Also note that we hope to get extra 2016-2017 data before we finalise this, 
    so probably don't want anything very specific (e.g. picking ENSO years for end and start) - end of 2015 is a big El Nino
    Could try 2015 - 1958 as a join instead
    Also could use 2050- to add as the last period, since the rate of warming varies over time
        or perhaps adjust the time in the CMIP trend to use for 2080-

       Inputs:
       cube: The cube of SST data from the start to the end of the historic period
       running_mean_file: file names containing the running mean of seasonal variability
    '''
    print 'start calc_future_sst'
    hadisst2_mask = iris.load_cube('/gws/nopw/j04/hrcm/cache/malcolm/HadISST2.2.2.0/hadisst_0to360_alldays_sst_1963.nc','sea_surface_temperature')[0]
    mask = np.where(hadisst2_mask.data.mask == True)

    cube.coord('time').bounds = None
    year_data_start = cube.coord('year').points[0] # 1952
    year_data_end = cube.coord('year').points[-1] # 2015

    future_sst = add_future_time_info(cube)
    #future2_sst = add_future_time_info(future_sst)
    #del cube
    #observed_trend = iris.load_cube(observed_trend_file)

    # load externally produced trend from CMIP5
    cmip5_trend = iris.load_cube(CMIP5_trend_file)
    try:
        icc.add_year(cmip5_trend, 'time', name = 'year')
    except:
        pass

    for future_c in [future_sst]:
        year1 = future_c.coord('year').points[0] # 2016
        year2 = future_c.coord('year').points[-1] # 2080

        for yr in range(year1, year2+1):
            adjust_jan = False  # edge of year period
            print ('processing for future year ',yr)
            future_file_output = os.path.join(SAVEDIR, 'sst', 'future_sst_'+str(yr)+'_025_daily_v1.nc')
            if os.path.exists(future_file_output):
                continue

            monthly_mean_hist_file = ''
            years_from_start = float(yr - year_data_start)
            year_hist = yr-year1+year_data_start
            year_constraint_hist = iris.Constraint(coord_values = {'year' : lambda l : l == year_hist})

            print ('extract year hist from running mean ',yr, year_hist, running_mean_file)
            running_mean_hist = iris.load_cube(running_mean_file)
            #running_mean_hist = monthly_mean_hist.extract(year_constraint_hist)
            future_yr_l = iris.cube.CubeList()

        # have months in the trend
            print 'years into future ',yr, years_from_start, year_hist
        
            year_constraint = iris.Constraint(coord_values = {'year' : lambda l : l == yr})
            cmip5_yr, cmip5_yrm1, cmip5_yrp1 = get_cmip5_this_year(yr, year_hist, cmip5_trend, year_constraint, year_constraint_hist)
                
            # calculate the variability here for hist year
    # calculate sst variability by removing monthly mean and observed trend
    # i.e. baseline variability (to add into future SST segment)
            print 'calculate variability, year, year_hist ',yr, year_hist
            sst_variability_mon, sst_variability = calc_variability(cube, running_mean_hist, year_hist, years_from_start)

            # find the most recent decadal mean for this year
            # if we are at the edge of a period, we need to adjust the value to get a smooth Dec-Jan transition
            year_constraint_end = iris.Constraint(coord_values = {'year' : lambda l : l == year_data_end-4})
            year_end_decade = year_data_end-4
            print ('extract year from running mean file for recent ',year_data_end, running_mean_file)
            running_mean_recent = iris.load_cube(running_mean_file)
            #running_mean_recent = running_mean_recent.extract(year_constraint_end)
            print ('use running mean ',running_mean_file, year_data_end)

            # if we are the first year of invented data, need to adjust the Jan value of monthly_mean_recent 
            # we also need to get the last day of variability from the previous year, to make the variability consistent across the boundary
            if yr == year1:
                print 'yr, year1, calc local means, anoms ',yr,year1
                monthly_mean_previous = adjust_monthly_mean_edge(running_mean_recent, datafiles[int(year_data_end)], adjust_dec = True)
                monthly_mean_next = running_mean_recent
                cube_2015 = iris.load_cube(datafiles[int(year_data_end)])
                cube_2015.convert_units('degC')
                icc.add_year(cube_2015, 'time')
                print 'cube 2015 ',cube_2015
                sst_variability_mon_yrm1, sst_variability_yrm1 = calc_variability(cube_2015, monthly_mean_previous, int(year_data_end), years_from_start, save = False, weighted = False)
                adjust_jan = True
            else:
                monthly_mean_previous = running_mean_recent
                monthly_mean_next = running_mean_recent
                sst_variability_yrm1 = []; sst_variability_mon_yrm1 = []

            for m in range(0,12):
                future_yrmon = get_daily_interpolated_data(m, sst_variability, sst_variability_mon, future_c, running_mean_recent, cmip5_yr, cmip5_yrm1, cmip5_yrp1, year_constraint, year_constraint_hist, monthly_mean_previous, monthly_mean_next, year_end_decade, sst_variability_yrm1 = sst_variability_yrm1, sst_variability_mon_yrm1 = sst_variability_mon_yrm1, adjust_jan = adjust_jan)
                future_yrmon.data.mask[:, mask[0], mask[1]] == True

                future_yr_l.append(future_yrmon)

            future_yr = future_yr_l.concatenate_cube()
            add_metadata(future_yr)
            savefile = os.path.join(future_file_output)
            iris.save(future_yr, savefile, unlimited_dimensions = ['time'], fill_value = 1.0e20)

    for c in [cube,future_sst]:
        try:
            c.remove_coord('year')
            c.remove_coord('month')
        except:
            pass

def calculate_time_components(datafiles, start_year, end_year):
    '''
    Decompose the time components of the historic data
    trend over each month
    variability
    observed trend
    '''
    files = []
    for year in range(start_year, end_year+1):
        files.append(datafiles[year])

    sst_cube_full_l = iris.load(files, callback = callback_attrib)
    print 'in calc time components, sst_cube_full_l ',sst_cube_full_l
    iris.util.describe_diff(sst_cube_full_l[0], sst_cube_full_l[10])
    sst_cube_full = sst_cube_full_l.concatenate_cube()

    # add various Aux coordinates to make it easier to access times
    try:
        add_time_names(sst_cube_full, year=True, month_number=True, day_of_year = True)
    except:
        pass

# extract all years between start_year and end_year
    sst_cube_t1_t2 = extract_cube_years(sst_cube_full, start_year, end_year)
    del sst_cube_full

    sst_cube_t1_t2.convert_units('degC')
    print 'sst_cube ',sst_cube_t1_t2

    # calculate the mean over each month
    monthly_mean_file = SAVEDIR+'/sst_monthly_mean_025_fixed_'+str(start_year)+'-'+str(end_year)+'.nc'
    print 'calculate monthly means into file ', monthly_mean_file
    if not os.path.exists(monthly_mean_file):
        monthly_mean = aggregate_mean(sst_cube_t1_t2, coord = 'month')
        print monthly_mean
        iris.save(monthly_mean, monthly_mean_file)
    else:
        monthly_mean = iris.load_cube(monthly_mean_file)

    # now decompose this cube into its components: 
    #    spatial time-mean, spatial linear trend, spatial variability

    # calculate linear trend for each month
    print 'calculate observed trend'
    observed_trend_file = SAVEDIR+'/observed_trend_025_fixed_'+str(start_year)+'-'+str(end_year)+'.nc'
#    if not os.path.exists(observed_trend_file):
#        sst_remove_monthmean = sf.remove_monthly_daily_mean_time_avg(sst_cube_t1_t2, monthly_mean = monthly_mean)
    
#        observed_trend = calc_observed_trend(sst_cube_t1_t2, sst_remove_monthmean)
#        print observed_trend
#        iris.save(observed_trend, observed_trend_file)
#    else:
#        observed_trend = iris.load_cube(observed_trend_file)

    return sst_cube_t1_t2, monthly_mean_file, observed_trend_file

def calculate_detrended(cube):
    data = cube.data
    order = float(len(data)/10.0)
    print order, data
    detrended = scipy.signal.detrend(data)
    cube.data = detrended
    return cube

def monthly_running_mean(datafiles, start_year, end_year):
    '''
    Try removing mean seasonal cycle from all data
    Then perhaps filtering the rest (or else extracting the trend for each month
    '''
    file_runningmean = SAVEDIR+'/hadisst2_running_monthly_mean_025_fixed_v1.nc'
    if os.path.exists(file_runningmean):
        return file_runningmean

    files = []
    for year in range(start_year-5, end_year+1):
        if year in datafiles:
            files.append(datafiles[year])

    sst_cube_full_l = iris.load(files, callback = callback_attrib)
    print 'in monthly running mean, sst_cube_full_l ',sst_cube_full_l
    c = sst_cube_full_l.concatenate_cube()

    # add various Aux coordinates to make it easier to access times
    add_time_names(c, year=True, month_number=True, day_of_year = True)
    # convert to SST if surface temperature
    c.convert_units('degC')

    cube_mean = sf.monthly_mean_running_time_avg(c)

    iris.save(cube_mean, file_runningmean, unlimited_dimensions = ['time'], fill_value = 1.0e20)

    return file_runningmean
    
def calculate_monthly_means(datafiles, start_year, end_year):
    '''
    Decompose the time components of the historic data
    Calculate monthly means over 10 years, 5 years apart
    '''
    monthly_mean_files = {}
    year_edges = []
    for y_central in range(start_year, end_year+1, 5):
        print('process monthly ',y_central)
        year_edges.append(y_central+5)
        files = []
        monthly_mean_files[y_central] = SAVEDIR+'/sst_monthly_mean_025_fixed_'+str(y_central)+'_v1.nc'
        if os.path.exists(monthly_mean_files[y_central]):
            continue

        # get years from start_year-5, start_year+5
        for year in range(y_central-4, y_central+5):
            if year in datafiles:
                files.append(datafiles[year])

        sst_cube_full_l = iris.load(files, callback = callback_attrib)
        print sst_cube_full_l
        print 'meaning central year ',y_central,', using files ',files
        print 'save to ',monthly_mean_files[y_central]
        sst_cube_full = sst_cube_full_l.concatenate_cube()

    # add various Aux coordinates to make it easier to access times
        add_time_names(sst_cube_full, year=True, month_number=True, day_of_year = True)


    # convert to SST if surface temperature
        sst_cube_full.convert_units('degC')
        #data_max = np.amax(sst_cube_full[0].data)
        #print 'data max ',data_max
        #if data_max >= 100.: sst_cube_full -= 273.15
        #print 'sst_cube ',sst_cube_full

    # need to detrend the period first, before calculating the monthly mean
        

    # calculate the mean over each month
        print 'calculate monthly means'
        monthly_mean = aggregate_mean(sst_cube_full, coord = 'month')
        print monthly_mean
        iris.save(monthly_mean, monthly_mean_files[y_central])

    return monthly_mean_files

def work():
    '''
    Main function to produce the future SST dataset
    '''
    runid = 'HadISST'
    start_year = 1976 # this is a La Nina, may cause problem with 2015 end
    start_year = 1980 # this is fairly neutral Nino
    end_year = 2015
    #end_year = 1980
    datafiles = {}
    for yr in range(start_year-5, end_year+1):
        fname = os.path.join(DATADIR, 'hadisst2_tos_daily_'+str(yr)+'_fixed_day_v1_monthly_under_ice.nc')
        if os.path.exists(fname):
            datafiles[yr] = fname

# Read SST data and calculate annual means
# the callback method only works well is there are full years of 12 months data
    print 'start analysis'
    #sst_cube_full = iris.load_cube(datafile, 'sea_surface_temperature', callback = sf.year_month_callback)
    #monthly_mean_files = calculate_monthly_means(datafiles, start_year, end_year)
    running_mean_file = monthly_running_mean(datafiles, start_year, end_year)

    sst_cube_t1_t2, monthly_mean_file, observed_trend_file = calculate_time_components(datafiles, start_year, end_year)

    # calculate future SST segment
    # by adding cmip5_trend into dataset
    print 'start calc_future_sst'
    #final_cube = calc_future_sst(sst_cube_t1_t2, observed_trend_file, monthly_mean_files, datafiles)
    final_cube = calc_future_sst(sst_cube_t1_t2, observed_trend_file, running_mean_file, datafiles)
  

if __name__ == '__main__':
    work()

