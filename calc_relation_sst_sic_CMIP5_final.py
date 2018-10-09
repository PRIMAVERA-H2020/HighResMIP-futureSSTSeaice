'''
Purpose:
    Fix HadISST2 SSTs (removing stripes in Arctic)
        by calculating relationship between SST and sea-ice concentration by month
    Produce monthly mean SST using this fix
    The monthly mean data will be used to produce the future SST values
    The SST-sea-ice-month relation will also be used to produce the future sea-ice concentrations
'''

import os, sys
import iris
import glob
import iris.coord_categorisation as icc
import numpy as np
import cPickle as pickle #cPickle is an optimised version of Pickle and is O(100) times faster
import calendar
from scipy.interpolate import interp1d

DATADIR = '/home/users/mjrobert/hrcm/cache/malcolm/HadISST2/1x1/processing_2018/'
BIN = {}; BIN['SST'] = {}; BIN['SICE'] = {}
BIN['SST']['MAX_MIN'] = [-2.8, 5.0]
BIN['SST']['SIZE'] = 0.02
BIN['SST']['XBINS'] = np.arange(BIN['SST']['MAX_MIN'][0], BIN['SST']['MAX_MIN'][1], BIN['SST']['SIZE'])
BIN['SST']['NBINS'] = len(BIN['SST']['XBINS'])-1

BIN['SICE']['MAX_MIN'] = [0.0, 101.5]
BIN['SICE']['SIZE'] = 0.5
BIN['SICE']['XBINS'] = np.arange(BIN['SICE']['MAX_MIN'][0], BIN['SICE']['MAX_MIN'][1], BIN['SICE']['SIZE'])
BIN['SICE']['NBINS'] = len(BIN['SICE']['XBINS'])-1

pole_lat_crit = {}; pole_lon_crit = {}
pole_lat = {}; pole_lon = {}
pole_lat['npole'] = [50., 90.]; pole_lat['spole'] = [-90.0, -50.0]
pole_lat_crit['npole'] = iris.Constraint(coord_values = {'latitude' : lambda l : 50 < (l.point) <= 90})
pole_lon_crit['npole'] = iris.Constraint(coord_values = {'longitude' : lambda l : 0 < (l.point) <= 360})
pole_lat_crit['spole'] = iris.Constraint(coord_values = {'latitude' : lambda l : -90 < (l.point) <= -50})
pole_lon_crit['spole'] = iris.Constraint(coord_values = {'longitude' : lambda l : 0 < (l.point) <= 360})
POLES = [(pole_lon_crit['npole'], pole_lat_crit['npole']), (pole_lon_crit['spole'], pole_lat_crit['spole'])]
pole_title = ['Arctic, daily HadISST2','Antarctic, daily HadISST2']
ICE_TEMP_MIN = -2.00
ICE_CONC_MAX = 100.0

CMIP5_ref = {'ACCESS1-0': 'CSIRO-BOM', 'ACCESS1-3': 'CSIRO-BOM', 'GFDL-CM3': 'NOAA-GFDL', 'IPSL-CM5A-MR': 'IPSL', 'CNRM-CM5': 'CNRM-CERFACS', 'HadGEM2-ES':'MOHC'}
#CMIP5_ref = {'CNRM-CM5': 'CNRM-CERFACS'}

def history_period_callback(cube, field, filename):
    # remove attributes preventing cube concatenation
    attrib = ['history','time_period']
    for att in attrib:
        try:
            del cube.attributes[att]
        except:
            pass

def interpolate_histogram(year, mean_freq):
    '''
    Interpolate the monthly histogram to generate daily values, for the number
    of days in this year
    '''
    iyear = int(year)
    if calendar.isleap(iyear):
        ndays = 366
    else:
        ndays = 365

    midmonth = np.zeros((14))
    midmonth[0] = 1
    daynumber = 0
    for im in range(0,12):
        monthdays = calendar.monthrange(iyear, im+1)[1]
        if im > 0:
            monthadd = calendar.monthrange(iyear, im)[1]
        else:
            monthadd = 0
        daynumber += monthadd
        midmonth[im+1] = monthdays // 2 + daynumber
    midmonth[-1] = ndays
    alldays = np.arange(1, ndays+1)

    mean_freq_tmp = np.ma.zeros((14))
    mean_freq_days = np.ma.zeros((BIN['SICE']['NBINS'], 2, ndays))
    
    print 'shapes ',mean_freq.shape, mean_freq_tmp.shape, midmonth.shape
    print midmonth

    for ir in range(0,2):
        for ib, bin in enumerate(BIN['SICE']['XBINS'][:-1]):
            mean_freq_tmp[1:13] = mean_freq[ib, ir, :]
            # set wrap points for year
            mean_freq_tmp[0] = (mean_freq[ib, ir, 0] + mean_freq[ib, ir, 11]) / 2.0
            mean_freq_tmp[13] = mean_freq_tmp[0]
            f = interp1d(midmonth, mean_freq_tmp, kind = 'linear')
            
            mean_freq_days[ib, ir, :] = f(alldays)

    return mean_freq_days
 

def calculate_histogram(x_values, x_histogram, y_values, xbins):
    # sea-ice relationship to SST
    sy, _ = np.histogram(x_values, bins=xbins, weights=y_values)
    sy2, _ = np.histogram(x_values, bins=xbins, weights=y_values*y_values)
    
    mean = sy / x_histogram          
    #print 'in histogram'
    #for ix, bin in enumerate(xbins[:-1]):
    #    print 'x_histogram ',bin, x_histogram[ix], sy[ix], mean[ix]

    std = np.sqrt(sy2/x_histogram - mean*mean)
    var = (sy2/x_histogram - mean*mean)
    return mean, var

def cube_histogram(cube_ind, cube_dep, nbins, xbins, xrange=[-1000,1000], yrange=[-1000,1000], sst_dep = False):
    '''
    Return a histogram of the frequency, mean and stdev of the independent cube
    Independent cube cube_ind
    Dependent cube cube_dep
    What I want to do:
        exclude masked values
        exclude values which are masked
        make sure to do this the same for both cubes
        would like to remove these values from the data, then can do histogram
    with the 2 cubes, I need to make sure that the mask is the same in both
    then I can use 
    data = cube.data
    data.compressed() - removes masked values
    data[~data.mask] - should do same
    '''
    histogram = np.ma.zeros((nbins))
    mean = np.ma.zeros((nbins))
    var = np.ma.zeros((nbins))
    mean[...] = np.ma.masked
    var[...] = np.ma.masked
    histogram[...] = np.ma.masked

    print 'cube shapes ',cube_ind.shape, cube_dep.shape
    if not "MaskedArray" in str(type(cube_dep)):
        cube_dep.data = np.ma.MaskedArray(cube_dep.data, shrink=False)
    cube_dep.data.mask = cube_ind.data.mask

    #mask = np.where((cube_ind.data.mask == True) | (cube_dep.data.mask == True))
    #mask = np.where((cube_ind.data.mask == True))
    #cube_ind.data.mask[mask[0], mask[1], mask[2]] = True
    #cube_dep.data.mask[mask[0], mask[1], mask[2]] = True

#        edge = np.where((cube_dep.data.mask == False) & (cube_dep.data > 0.0) & (cube_dep.data < 15.0) & (cube_ind.data > 0.5))
#        cube_dep.data[edge[0], edge[1], edge[2]] = 0.01
        
    cube_ind_data = cube_ind.data
    cube_dep_data = cube_dep.data
    
    #------------Reshape to 1D arrays--------------------------------------------
    x1 = cube_ind_data.compressed()
    y1 = cube_dep_data.compressed()
    #x1 = np.reshape(cube_ind_data, -1)
    #y1 = np.reshape(cube_dep_data, -1)

    
    print 'shape ',x1.shape, y1.shape
    print 'max x1 ',np.amax(x1), np.amin(x1)
    print 'max y1 ',np.amax(y1), np.amin(y1)

    data_ind = x1
    data_dep = y1
    
    if data_ind.shape != data_dep.shape:
        print 'xshape, yshape ',data_ind.shape, data_dep.shape
        sys.exit('sst_shape shape not equal to wind_shape shape')
    
    max_value = data_dep.max()
    min_value = data_dep.min()
    
    # make a histogram using SST
    hist, _ = np.histogram(data_ind, bins=xbins)
    histogram[...] = hist

    # make minimum 1 since divide by number of elements
    ma = np.where(histogram < 1)[0]
    
    mean, var = calculate_histogram(data_ind, histogram, data_dep, xbins)

    if not sst_dep:
        mean[ma] = max_value
    else:
        mean[ma] = min_value
            
    return histogram, mean, var, max_value, data_ind, data_dep

def extract_sst_sice_data(sst, pole, pole_lat, model):
    '''
    Extract data when not regular grid (and hence cannot use iris.extract)
    Simply mask the non-relevant points
    '''
    sst_copy = sst.copy()
    print 'pole_lat ',pole_lat
    if pole == 'spole':
        # SH
        miss = np.where(sst.coord('latitude').points[:] > pole_lat[1])
        if not 'HadGEM' in model:
            sst_copy.data.mask[:, miss[0], miss[1]] = True
        else:
            sst_copy.data.mask[:, miss[0], :] = True
    else:
        # NH
        miss = np.where(sst.coord('latitude').points < pole_lat[0])
        if not 'HadGEM' in model:
            sst_copy.data.mask[:, miss[0], miss[1]] = True
        else:
            sst_copy.data.mask[:, miss[0], :] = True
    return sst_copy
                   

def sst_ice_relationship(sst, ice):
    '''
    Generate the relationship between sst, sea-ice conc and month
    '''
    histogram_freq_sst = np.ma.zeros((BIN['SST']['NBINS'], 2))
    histogram_freq_sst[...] = np.ma.masked
    histogram_freq_sice = np.ma.zeros((BIN['SICE']['NBINS'], 2))
    histogram_freq_sice[...] = np.ma.masked

    mean_freq_sst = np.ma.zeros((BIN['SST']['NBINS'], 2))
    mean_freq_sice = np.ma.zeros((BIN['SICE']['NBINS'], 2))

    #for ir, (pole_lon, pole_lat) in enumerate(POLES):
    for ir, pole in enumerate(['npole', 'spole']):
        cube_sst = extract_sst_sice_data(sst, pole, pole_lat[pole], model)
        cube_sice = ice.copy()

        print 'cube_sst max/min',cube_sst.data.max(), cube_sst.data.min()
        print 'cube_sice max/min',cube_sice.data.max(), cube_sice.data.min()

        nbins = BIN['SST']['NBINS']
        xbins = BIN['SST']['XBINS']
        print 'calc sst = f(sice), ir =  ',ir
        freq, mean, var, max_value, x_data, y_data = \
            cube_histogram(cube_sst, cube_sice, nbins, xbins)
        mean_freq_sst[:, ir] = mean

#    for ir, (pole_lon, pole_lat) in enumerate(POLES):
#        cube_sst = sst.extract(pole_lat)
#        cube_sst = cube_sst.extract(pole_lon)
#        cube_sice = ice.extract(pole_lat)
#        cube_sice = cube_sice.extract(pole_lon)

#        print 'calc sice = f(sst), ir =  ',ir
#        nbins = BIN['SICE']['NBINS']
#        xbins = BIN['SICE']['XBINS']
#        freq, mean, var, max_value, x_data, y_data = \
#            cube_histogram(cube_sice, cube_sst, nbins, xbins, sst_dep = True)
#        mean_freq_sice[:, ir] = mean
        

    return mean_freq_sst

def calc_sst_func_sice_relationship(dir_in, model, filename_sst_func_siconc_month_pickle, years):

    ftype = 'tos'
    sst_fixed = os.path.join(DATADIR, model+'_monthly_1950_2100_tos.nc')
    fnames_sst = os.path.join(dir_in, model+'_monthly_1950_2100_tos.nc')
    fnames_ice = os.path.join(dir_in, model+'_monthly_1950_2100_sic.nc')

    print 'choose ',fnames_sst

    cube_mon = iris.cube.CubeList()

    mean_freq_sst = np.ma.zeros((BIN['SST']['NBINS'], 2, 12))
    
    print 'process ',fnames_sst, fnames_ice
    sst1 = iris.load_cube(fnames_sst)
    sst1.convert_units('degC')
    ice1 = iris.load_cube(fnames_ice)
    sst_coords = [co.name() for co in sst1.aux_coords]
    ice_coords = [co.name() for co in ice1.aux_coords]
    if not 'year' in sst_coords:
        icc.add_year(sst1, 'time', name = 'year')
    if not 'year' in ice_coords:
        icc.add_year(ice1, 'time', name = 'year')

    years_con = iris.Constraint(coord_values = {'year' : lambda l : years[0] <= l <= years[-1]})
    sst = sst1.extract(years_con)
    ice = ice1.extract(years_con)

    if not 'month' in sst_coords:
        icc.add_month_number(sst, 'time', name = 'month')
    if not 'month' in ice_coords:
        icc.add_month_number(ice, 'time', name = 'month')

    for im in range(0,12):
        print 'processing month ',im
        month = im+1
        month_no = iris.Constraint(coord_values = {'month' : lambda l : l == month})
        sst_mon = sst.extract(month_no)
        ice_mon = ice.extract(month_no)

        mean_freq_sst_month_year = sic_functions.sst_ice_relationship(sst_mon, ice_mon, do_freq_sst = True)
        mean_freq_sst[:, :, im] = mean_freq_sst_month_year

        mean_freq_sst[0, :, im] = 100.0

        ir = 0
        for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
                print 'nh, im, ir, bin, siconc ',im, ir, bin, mean_freq_sst[ib, ir, im]

    do_mono = True

    if do_mono:
    # want to make the relationship monotonic
        for im in range(0,12):
            reset_nh = False; reset_sh = False
            for ir in range(0,2):
                for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
                    if ib > 0 and ib < BIN['SST']['NBINS']-2:
                        current_val = mean_freq_sst[ib, ir, im]
                        max_rest = np.amax(mean_freq_sst[ib+1:, ir, im])
                        if max_rest > current_val:
                            mean_freq_sst[ib, ir, im] = max_rest
                # once the siconc for a given sst hits zero, make sure all the rest are zero
                if mean_freq_sst[ib, 0, im] == 0:
                    reset_nh = True
                if reset_nh:
                    mean_freq_sst[ib:, 0, im] = 0.0
                if mean_freq_sst[ib, 1, im] == 0:
                    reset_sh = True
                if reset_sh:
                    mean_freq_sst[ib:, 1, im] = 0.0

    # try making the NH concentration less than 15% equal to zero (to remove wide edge in summer)
#    min_ice_conc = 20.0
#    for im in range(0,12):
#        for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
#            if ib > 0 and ib < BIN['SST']['NBINS']-2:
#                for ir in range(0,1):
#                    if (mean_freq_sst[ib, ir, im] < min_ice_conc):
#                            mean_freq_sst[ib+4:, ir, im] = 0.0
    

    im = 4
    ir = 0
    for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
        print 'after mono, nh, im, bin, siconc ',im, bin, mean_freq_sst[ib, ir, im]

    fh = open(filename_sst_func_siconc_month_pickle, 'wb')
    pickle.dump(mean_freq_sst, fh)
    fh.close()
    

def process_sst_seaice(dir_in, years, months, year_last_real_data, sst_fixed_year, model, sst_ice_only = False):

    cube_mon = iris.cube.CubeList()
    # 271.32 is what HadISST2 uses under full seaice
    #ice_temp = 271.32 # 273.15 - 1.8

    #filename_sice_func_sst_month_pickle = os.path.join(DATADIR, 'siconc_func_sst_month_relationship_'+model+'.pkl')
    filename_sst_func_siconc_month_pickle = os.path.join(DATADIR, 'sst_func_siconc_month_relationship_'+model+'.pkl')

    # then calculate sst func sice using fixed SST under sea-ice

    fnames_sst = os.path.join(dir_in, model+'_monthly_1950_2100_tos.nc')
    fnames_ice = os.path.join(dir_in, model+'_monthly_1950_2100_sic.nc')

    sic_functions.calc_sst_func_sice_relationship(dir_in, model, filename_sst_func_siconc_month_pickle, years, fnames_sst, fnames_ice)

if __name__ == '__main__':
    dir_in = '/home/users/mjrobert/hrcm/cache/malcolm/HighResMIP/sst_forcing/processing_ocean_v2/'
    year_last_real_data = '2015'
    #years = range(2015, 2016)
    years = range(2010, 2020)
    months = range(0,12)
    #years = range(2000, 2014)
    #months = range(0,1)

    model = CMIP5_ref
    for model in CMIP5_ref:
        sst_fixed = os.path.join(DATADIR, model+'_monthly_1950_2100_tos.nc')
        sst_fixed_year = sst_fixed
        process_sst_seaice(dir_in, years, months, year_last_real_data, sst_fixed_year, model, sst_ice_only = False)
    

                       
