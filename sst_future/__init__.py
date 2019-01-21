#!/usr/local/sci/bin/python2.7
'''
DESCRIPTION
    Helper functions to produce the future SST forcing for HighResMIP

AUTHOR:
    Malcolm Roberts (hadom)

LAST MODIFIED:
    2017-03-10

'''

import iris.coord_categorisation as icc
import os, sys
import numpy as np
import iris
import iris.analysis
import scipy.linalg
import subprocess
from scipy.stats import norm
import datetime

verbose = True
# chunks of data to get at once
interval = 30
colourname = ['red','blue','green','orange','aqua','lime','sienna','black']

def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2 ) + 1
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


def year_month_callback(cube, field, filename):
    # add a categorical year coordinate to the cube
    try:
        cube.remove_coord('forecast_period')
        cube.remove_coord('forecast_reference_time')
    except:
        pass

def reshape_data(dat,mdi,weights=None):
    '''
    Returns matrix suitable for eof analysis.
    '''
    nt,ny,nx = dat.shape
    matr = dat.reshape((nt,ny*nx))
    matr_index = [slice(None)]*matr.ndim
    ind = np.where(matr[0] != mdi)[0]
    matr_index[1] = ind
    M = matr[matr_index]
    if weights != None:
        wts = weights.reshape(ny*nx)[ind]
        wts = wts/wts.max()
        M = M*wts
    return M 
    

def centre_columns(M):
    '''
    Removes time mean of each matrix column
    
    '''
    return M - np.mean(M,axis=0)

def mask_ice(cube):
    '''
    Updates mask to cover ice
    '''
    if not "MaskedArray" in str(type(cube.data)):
        cube.data = np.ma.MaskedArray(cube.data, shrink=False)
    
    cube.data.mask = cube.data.data < -2
    return cube

def mask_land(cube, resol_a, resol_o):
    '''
    Updates mask to land (if model surface temperature)
    Adds the mask in case something else (like ice) is already masked
    '''
    maskdir = '/project/hadgem3/hadom/hadgem3h/'
    mask = iris.load_cube(os.path.join(maskdir,'mask_frac_'+resol_a+'_'+resol_o+'_ref.pp'))
    mask.data = np.ma.MaskedArray(mask.data, mask.data > 0)
    if not "MaskedArray" in str(type(cube.data)):
        cube.data = np.ma.MaskedArray(cube.data, shrink=False)
#    cube.data.mask = mask.data > 0.
    cube.data.mask += mask.data.mask
    #print cube.data.mask
    return cube

def remove_monthly_mean_time_avg(cube):
    '''
    Removes monthly mean time mean from value
    '''
    monthly_mean = cube.aggregated_by('month', iris.analysis.MEAN)
    print 'monthly_mean ',monthly_mean
    
    sst_nomonth = cube.copy()
    for m in range(0,12):
        month_no = np.where(sst_nomonth.coord('month').points == m+1)[0]
        for index in month_no:
            sst_nomonth[index].data[:,:] = sst_nomonth[index].data[:,:] - monthly_mean[m].data[:,:]
    return sst_nomonth

def monthly_mean_running_time_avg(cube):
    '''
    Removes monthly mean time mean from value
    '''
    cube_mm = cube.aggregated_by(['month','year'], iris.analysis.MEAN)
    cube_running = cube_mm.copy()

    year_min = cube_mm.coord('year').points.min()
    year_max = cube_mm.coord('year').points.max()
    years = range(year_min, year_max+1)
    for year in years:
        # extract years 5 years either side of this year
        # calculate mean seasonal cycle
        # calculate anomaly this year
        # may be I just want to mean centred on this year
        year_from = np.max([year-4, year_min])
        year_to = np.min([year+4, year_max])
        year_con = iris.Constraint(coord_values = {'year' : lambda l : year_from <=l <= year_to})
        cube_ext = cube_mm.extract(year_con)
        monthly_mean = cube_ext.aggregated_by('month', iris.analysis.MEAN)
    
        year_no = np.where(cube_running.coord('year').points == year)[0]
        cube_running.data[year_no, :, :] = monthly_mean.data[:, :, :]

    return cube_running


def remove_monthly_daily_mean_time_avg(cube, monthly_mean = '', window = 60):
    '''
    Removes monthly mean time mean from cube with matching months
    Include monthly mean in argument list to not recalculate it
    window is the width of the rolling window 
    '''
    if monthly_mean == '':
        monthly_mean = cube.aggregated_by('month', iris.analysis.MEAN)
    print 'monthly_mean ',monthly_mean
    
    sst_nomonth = cube.copy()
    #sst_nomonth = cube.rolling_window('time', iris.analysis.MEAN, window)

    for m in range(0,12):
        month_no = np.where(sst_nomonth.coord('month').points == m+1)[0]
        for index in month_no:
            sst_nomonth[index].data[:,:] = sst_nomonth[index].data[:,:] - monthly_mean[m].data[:,:]
    return sst_nomonth

def remove_glob_time_avg(cube):
    '''
    Removes time mean from value
    '''
    #glob_areas = guess_areas(cube)
    time_avg = cube.collapsed('time', iris.analysis.MEAN)
    for nt in xrange(cube.shape[0]):
        cube.data[nt,:,:] = cube.data[nt,:,:] - time_avg.data
    return cube, time_avg

def replace_glob_time_avg(cube, time_avg):
    '''
    Removes time mean from value
    '''
    #glob_areas = guess_areas(cube)
    for nt in xrange(cube.shape[0]):
        cube.data[nt,:,:] = cube.data[nt,:,:] + time_avg.data
    return cube

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

def replace_monthly_avg(cube, sst_month):
    '''
    Removes monthly mean from each mon/ann value
    '''
    print 'sst_month ',sst_month
    print 'cube ',cube
    print 'month number ',cube.coord('month')
    print 'month number ',len(cube.coord('month').points)
    monthly_full = cube.copy()
    for m in range(0,11+1):
        inds = np.where(cube.coord('month').points == m+1)[0]
        for i in inds:
            monthly_full.data[i] += sst_month.data[m]
    
    return monthly_full

def remove_global_trend(cube, trend_fit):
    '''
    Removes global mean from each mon/ann value
    '''
    indices = np.where(cube.data.mask[0,:,:] == False)
    print indices[0]
    for nt in xrange(cube.shape[0]):
        cube.data[nt,:,:] = cube.data[nt,:,:] - trend_fit[nt]
    return cube

def add_trend(cube, trend_fit):
    '''
    Removes global mean from each mon/ann value
    '''
    for nt in xrange(cube.shape[0]):
        cube.data[nt,:,:] = cube.data[nt,:,:] + trend_fit[nt]
    return cube

def remove_coords(cube):
    '''
    Removes coordinates that are not needed
    '''
    try:
        cube.remove_coord('forecast_period')
        cube.remove_coord('forecast_reference_time')
        cube.remove_coord('year')
        cube.remove_coord('month')
    except:
        pass
    return cube

def write_cube_netcdf(cube,savefile):
    '''
    Write cubedata to netcdf
    '''
    print 'Saving: '+savefile
    iris.fileformats.netcdf.save(cube,savefile,netcdf_format='NETCDF4')

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
#    if cube.coord('latitude').bounds is None:
    try:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
    except:
        pass

#    if cube.coord('longitude').bounds is None:
#        cube.coord('longitude').guess_bounds()
        
    #slice = cube.slices(['latitude','longitude']).next()
    #grid_areas = iris.analysis.cartography.area_weights(slice)
    grid_areas = iris.analysis.cartography.area_weights(cube)

    return grid_areas

    
def derive_variability(cube, runid, lon_range_pdo, lon_range_amo, savedir, remove_global_means, neofs, period_title=''):
    sst_cube = cube.copy()
    # mask the sea-ice regions
    sst_cube = mask_ice(sst_cube)
    
    sst_cube_monthanom, tmp_avg = remove_monthly_avg(sst_cube)
    tpi_ts = tpi_timeseries(sst_cube_monthanom)
    tpi_smooth = tpi_ts['tpi'].copy()
    smooth_data = utilities_mr.smooth(tpi_ts['tpi'].data, window_len=156, window = 'hanning')
    print len(smooth_data)
    tpi_smooth.data[:len(smooth_data)] = smooth_data

    #print 'tpi', tpi_ts['tpi']
    #print 'reg2', tpi_ts['reg2']

    file_out = os.path.join(savedir, runid+'_tpi_timeseries.nc')
    tpi_ts['tpi'].long_name = 'TPI timeseries'
    iris.save(tpi_ts['tpi'], file_out, netcdf_format="NETCDF3_CLASSIC")

    #print 'sst cube for nino ',sst_cube
    try:
        icc.add_season(sst_cube, 'time', name='clim_season')
        icc.add_season_year(sst_cube, 'time', name='season_year')
    except:
        pass
    #print 'sst cube for nino ',sst_cube
    #print 'sst cube for nino ',sst_cube.coord('clim_season')
    
    sst_djf = sst_cube.extract(DJF_constraint)
    nino_timeseries = nino34_timeseries(sst_djf)
    nino_timeseries.long_name = 'NINO3.4 timeseries'
    file_out = os.path.join(savedir, runid+'_nino34_timeseries'+period_title+'.nc')
    iris.io.save(nino_timeseries,file_out, netcdf_format="NETCDF3_CLASSIC")

    # calculate annual means
    sst_ann = sst_cube.aggregated_by('year', iris.analysis.MEAN)

# mask out land when using surface temperature
    if remove_global_means: sst_ann = remove_glob_avg(sst_ann)
    
# calculate PDO EOfs, timeseries and regression
    pdo_pc = calculate_pdo_eof_timeseries(runid, sst_ann, savedir, lat_range_pdo, lon_range_pdo, neofs, remove_global_means, period_title)
    #print 'pdo_pc ',pdo_pc

# calculate PDO SST regression and write file
    pltdir = os.path.join(savedir, 'pdo_patterns/plts/gc2/')
    if not os.path.exists(pltdir):
        os.makedirs(pltdir)
    pdofile = runid+'_pdo_eof1_pc_time_series_glob_mean_sst_removed'+period_title+'.nc'
    pdo_cube = iris.load(savedir + pdofile)[0]
    calc_regression(sst_ann, pltdir, pdo_cube, runid, 'PDO', period_title)

#calculate AMO timeseries
    amo_timeseries = calc_amo_timeseries(sst_ann, runid, lon_range_amo, lat_range_amo)
    amo_timeseries.long_name='AMO timeseries'
    file_out = os.path.join(savedir, runid+'_amo_timeseries_trenberth'+period_title+'.nc')
    iris.io.save(amo_timeseries,file_out, netcdf_format="NETCDF3_CLASSIC")

    
#calculate AMO regression and write files
    amodir = os.path.join(savedir, 'amo')
    if not os.path.exists(amodir): os.makedirs(amodir)
    calculate_amo_regression(sst_ann, amo_timeseries, lon_range_amo, lat_range_amo, amodir, runid, period_title)

    resol_model = ['HadISST']; title = 'HadISST'; desc_model = ['Nino3.4','AMO','PDO','TPI']; 
    ymin = [-3, -0.3, -0.2, -1.]
    ymax = [3, 0.3, 0.2, 1.]
    fig = plt.figure(figsize=(8,11),dpi=100)#dpi=300
    for i, ts in enumerate([nino_timeseries, amo_timeseries, pdo_pc, tpi_ts['tpi']]):
        subpl=fig.add_subplot(4,1,i+1)
        period = 'year'
        if i == 3: period = 'month'
        plot_timeseries(i, runid, ts, resol_model[0], desc_model[i], title, \
                        subpl, ymin[i], ymax[i], period = period)
    plt.savefig(os.path.join(savedir,runid+'_modes_timeseries'+period_title+'.png'))
    plt.show()
    del sst_cube

def add_buffer_months(cube_past, cube_future):
    w1 = np.zeros(24); w2 = np.zeros(24)
    w1[0] = 0.9; w1[1] = 0.6; w1[2] = 0.3
    w2[21] = 0.3; w2[22] = 0.6; w2[23]= 0.9
    #print 'mask ',cube_past.data.mask
    
    # for each non-land point in cube
    for j in range(cube_past.shape[1]):
        for i in range(cube_past.shape[2]):
            if cube_past.data.mask[0 , j , i] == False:
            #if not cube_past.data[0 , j , i] == -1.0e30:
                mu, std = norm.fit(cube_past.data[:, j, i])
                ts = norm.rvs(loc=mu, scale=std, size=24)
                tsnew = np.zeros(24)
                val0 = cube_past.data[-1, j, i]
                val1 = cube_future.data[24, j, i]
                for i in range(24): tsnew[i] = ts[i]*(1. - np.amax([w1[i], w2[i]])) + w1[i]*val0 + w2[i]*val1
                cube_future.data[0:24, j, i] = tsnew
                cube_future.data.mask[0:24, j, i] = False
            else:
                pass
                
    return cube_future

def calculate_sst_trend(cube):
    cube_yr = cube.aggregated_by('year', iris.analysis.MEAN)
    grid_areas = guess_areas(cube_yr)
    cube_ts = cube_yr.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
    coeffs = np.polyfit(cube_ts.coord('time').points, cube_ts.data, 1)
    print 'coeffs ',coeffs
    trend_fit = np.polyval(coeffs, cube_ts.coord('time').points)
    trend_per_decade = (trend_fit[-1] - trend_fit[-2]) * 10.
    print 'trend_fit', trend_fit, trend_per_decade
    print 'trend_per_decade', trend_per_decade
    del cube_yr
                
    # get timeseries at that point
    # calculate mean and std
    # generate timeseries with ts = norm.rvs(loc=mu, scale=std, size=24)
    
    # val0= = last point in first timeseries
    # val1 = first point in new timeseries
    # for i in range(24): tsnew[i] = ts[i]*(1. - np.amax([w1[i], w2[i]])) + w1[i]*val0 + w2[i]*val1
    # use these points to join together two timeseries
    # return cube

def fit_trend(cube, order):
    # fit a trend (least squares fit)
    grid_areas = guess_areas(cube)
    cube_ts = cube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
    coeffs = np.polyfit(cube_ts.coord('time').points, cube_ts.data, order)
    print 'coeffs for global trend',coeffs
    trend_fit = np.polyval(coeffs, cube_ts.coord('time').points)
    return trend_fit

def fit_trend_per_point(cube_past):
    if not "MaskedArray" in str(type(cube_past.data)):
        cube_past.data = np.ma.MaskedArray(cube_past.data, shrink=False)
        cube_past.data.mask = False
    cube_trend = cube_past[0].copy()
    #print cube_past.coord('year').points
    for j in range(cube_past.shape[1]):
        for i in range(cube_past.shape[2]):
            if cube_past.data.mask[0 , j , i] == False:
            #if not cube_past.data[0 , j , i] == -1.0e30:
                coeffs = np.polyfit(range(len(cube_past.coord('year').points)), cube_past.data[:, j, i], 1)
                #print 'coeffs ',coeffs
                cube_trend.data[j,i] = coeffs[0]
                #trend_fit = np.polyval(coeffs, range(len(cube_past.coord('year').points)))
                #print j,i,10.*(trend_fit[-1] - trend_fit[0]) / float(len(cube_past.coord('year').points))
                #cube_trend.data[j,i] = 10.*(trend_fit[-1] - trend_fit[0]) / float(len(cube_past.coord('year').points))
                #cube_trend.data[j,i] = (trend_fit[-1] - trend_fit[0])

            else:
                pass
                
    return cube_trend

def remove_point_trend(cube):
    '''
    Removes trend at each spatial point
    '''
    cube_trend_fit = cube.copy()
    cube_trend_coeffs = cube[0:3].copy()
    for j in range(cube.shape[1]):
        for i in range(cube.shape[2]):
            if cube.data.mask[0 , j , i] == False:
                coeffs = np.polyfit(range(len(cube.coord('year').points)), cube.data[:, j, i], 2)
                trend_fit = np.polyval(coeffs, range(len(cube.coord('year').points)))
                #print j,i,10.*(trend_fit[-1] - trend_fit[0]) / float(len(cube_past.coord('year').points))
                cube_trend_fit.data[:,j,i] = trend_fit
                cube_trend_coeffs.data[:,j,i] = coeffs

    cube -= cube_trend_fit
    del cube_trend_fit
    return cube, cube_trend_coeffs

def add_point_trend(cube, trend_coeffs):
    '''
    Removes global mean from each mon/ann value
    '''
    for j in range(cube.shape[1]):
        for i in range(cube.shape[2]):
            if cube.data.mask[0 , j , i] == False:
                coeffs = trend_coeffs.data[:,j,i]
                trend_fit = np.polyval(coeffs, range(len(cube.coord('year').points)))
                cube.data[:, j, i] += trend_fit
    
    return cube


if __name__ == '__main__':
    pass
