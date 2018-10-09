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
import sic_functions

DATADIR = '/home/users/mjrobert/hrcm/cache/malcolm/HadISST2/1x1/processing_2018/'
BIN = define_histogram_bins()

npole_lat = iris.Constraint(coord_values = {'latitude' : lambda l : 50 < (l.point) <= 90})
npole_lon = iris.Constraint(coord_values = {'longitude' : lambda l : 0 < (l.point) <= 360})
spole_lat = iris.Constraint(coord_values = {'latitude' : lambda l : -90 < (l.point) <= -50})
spole_lon = iris.Constraint(coord_values = {'longitude' : lambda l : 0 < (l.point) <= 360})
POLES = [(npole_lon, npole_lat), (spole_lon, spole_lat)]
pole_title = ['Arctic, daily HadISST2','Antarctic, daily HadISST2']
ICE_TEMP_MIN = -2.00
ICE_CONC_MAX = 100.0


def history_period_callback(cube, field, filename):
    # remove attributes preventing cube concatenation
    attrib = ['history','time_period']
    for att in attrib:
        try:
            del cube.attributes[att]
        except:
            pass

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

    for ir, (pole_lon, pole_lat) in enumerate(POLES):
        cube_sst = sst.extract(pole_lat)
        cube_sst = cube_sst.extract(pole_lon)
        cube_sice = ice.extract(pole_lat)
        cube_sice = cube_sice.extract(pole_lon)

        #print 'subcube shape ',cube_sst

        nbins = BIN['SST']['NBINS']
        xbins = BIN['SST']['XBINS']
        print 'calc sst = f(sice), ir =  ',ir
        freq, mean, var, max_value, x_data, y_data = \
            sic_functions.cube_histogram(cube_sst, cube_sice, nbins, xbins)
        mean_freq_sst[:, ir] = mean

    for ir, (pole_lon, pole_lat) in enumerate(POLES):
        cube_sst = sst.extract(pole_lat)
        cube_sst = cube_sst.extract(pole_lon)
        cube_sice = ice.extract(pole_lat)
        cube_sice = cube_sice.extract(pole_lon)

        print 'calc sice = f(sst), ir =  ',ir
        nbins = BIN['SICE']['NBINS']
        xbins = BIN['SICE']['XBINS']
        freq, mean, var, max_value, x_data, y_data = \
            sic_functions.cube_histogram(cube_sice, cube_sst, nbins, xbins, sst_dep = True)
        mean_freq_sice[:, ir] = mean
        
    return mean_freq_sst, mean_freq_sice

def fix_sst_based_on_siconc(sst, ice, mean_freq_sice, im):
    '''
    Loop through ice concentration bins
    Find points in which ice concentration is near to this bin
    Only do this when the sea-ice conc is > 60% (we don't whant to effect the edge, just the Arctic to remove the stripes
    '''
    lat_size = ice.shape[1]
    ice_min_change_nh = 75.0
    ice_min_change_sh = 75.0

    start_bin = int(ice_min_change_nh / BIN['SICE']['SIZE'])-1

    coord_names = [coord.name() for coord in sst.coords()]
    print coord_names
    if 'day_of_year' not in coord_names:
        icc.add_day_of_year(sst, 'time', name = 'day_of_year')

    ndays = sst.shape[0]
    day_of_year = sst.coord('day_of_year').points
    print 'in fix_sst ', day_of_year

    for day in range(ndays):
        day_number = day_of_year[day] - 1
        print ('process day, month ',day, im)
        for ib, bin in enumerate(BIN['SICE']['XBINS'][:-1]):
            if bin < ice_min_change_nh-1:
                continue
            ice_range = np.where((np.abs(ice.data[day, :, :] - bin) < BIN['SICE']['SIZE']) & (ice.data[day, :, :] > ice_min_change_nh) & (ice.data.mask[day, :, :] == False))
            nhemi = np.where(ice_range[0] > lat_size / 2)
            print 'nhemi ',nhemi, ib, bin
            sst.data[day, ice_range[0][nhemi], ice_range[1][nhemi]] = mean_freq_sice[ib, 0, day_number]

            ice_range_sh = np.where((np.abs(ice.data[day, :, :] - bin) < BIN['SICE']['SIZE']) & (ice.data[day, :, :] > ice_min_change_sh) & (ice.data.mask[day, :, :] == False))
            shemi = np.where(ice_range_sh[0] < lat_size / 2)
            print 'shemi ',shemi, ib, bin
            sst.data[day, ice_range_sh[0][shemi], ice_range_sh[1][shemi]] = mean_freq_sice[ib, 1, day_number]
    return sst

def calc_sice_func_sst_relationship(dir_in, filename_sice_func_sst_month_pickle, year):

    fnames_sst = sorted(glob.glob(os.path.join(dir_in, 'tos*'+year+'0101-*.nc')))
    fnames_ice = sorted(glob.glob(os.path.join(dir_in, 'siconc*'+year+'0101-*.nc')))

    print 'choose ',fnames_sst

    NYEARS = len(fnames_sst)
    print 'NYEARS ',NYEARS

    cube_mon = iris.cube.CubeList()

    mean_freq_sice = np.ma.zeros((BIN['SICE']['NBINS'], 2, 12))
    
    for f_sst, f_ice in zip(fnames_sst, fnames_ice):
        year = f_sst.split('_')[-1][0:4]
        year1 = f_ice.split('_')[-1][0:4]
        print 'process ',f_sst, year
        if year != year1:
            raise Exception('Paired SST and sea-ice years no same '+year+' '+year1)
        sst = iris.load_cube(f_sst)
        ice = iris.load_cube(f_ice)
        icc.add_month_number(sst, 'time', name = 'month')
        icc.add_month_number(ice, 'time', name = 'month')

        for im in range(0,12):
            print 'processing month ',im
            month = im+1
            month_no = iris.Constraint(coord_values = {'month' : lambda l : l == month})
            sst_mon = sst.extract(month_no)
            ice_mon = ice.extract(month_no)

            mean_freq_sst_month_year, mean_freq_sice_month_year = sic_functions.sst_ice_relationship(sst_mon, ice_mon, do_freq_sst = True, do_freq_sic = True)
            mean_freq_sice[:, :, im] = mean_freq_sice_month_year

    # want to make the relationship monotonic
    # temperature decreases as siconc bin increases
    for im in range(0,12):
        reset_nh = False; reset_sh = False
        for ib, bin in enumerate(BIN['SICE']['XBINS'][:-1]):
            if ib > 0 and ib < BIN['SICE']['NBINS']-2:
                for ir in range(0,2):
                    #if (mean_freq_sice[ib, ir, im] > np.amax(mean_freq_sice[ib:ib+5, ir, im])):
                    if (mean_freq_sice[ib, ir, im] > (mean_freq_sice[ib-1, ir, im])):
                        mean_freq_sice[ib, ir, im] = mean_freq_sice[ib-1, ir, im]

    im = 0
    ir = 0
    for ib, bin in enumerate(BIN['SICE']['XBINS'][:-1]):
        print 'nh, im, bin, sst ',im, bin, mean_freq_sice[ib, ir, im]

    fh = open(filename_sice_func_sst_month_pickle, 'wb')
    pickle.dump(mean_freq_sice, fh)
    fh.close()

def calc_sst_func_sice_relationship(dir_in, filename_sst_func_siconc_month_pickle, year):

    print 'choose ',fnames_sst

    NYEARS = len(fnames_sst)
    print 'NYEARS ',NYEARS

    cube_mon = iris.cube.CubeList()

    mean_freq_sst = np.ma.zeros((BIN['SST']['NBINS'], 2, 12))
    
    for f_sst, f_ice in zip(fnames_sst, fnames_ice):
        print 'process ',f_sst, year
        sst = iris.load_cube(f_sst)
        ice = iris.load_cube(f_ice)
        sst_coords = [co.name() for co in sst.aux_coords]
        ice_coords = [co.name() for co in ice.aux_coords]
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

            mean_freq_sst_month_year, mean_freq_sice_month_year = sic_functions.sst_ice_relationship(sst_mon, ice_mon, do_freq_sst = True, do_freq_sic = True)
            mean_freq_sst[:, :, im] = mean_freq_sst_month_year

            ir = 1
            for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
                print 'nh, im, ir, bin, siconc ',im, ir, bin, mean_freq_sst[ib, ir, im]

    # want to make the relationship monotonic
    for im in range(0,12):
        reset_nh = False; reset_sh = False
        for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
            if ib > 0 and ib < BIN['SST']['NBINS']-2:
                for ir in range(0,2):
                    if (mean_freq_sst[ib, ir, im] < np.amax(mean_freq_sst[ib:ib+5, ir, im])):
                        mean_freq_sst[ib, ir, im] = mean_freq_sst[ib-1, ir, im]
                    elif (mean_freq_sst[ib, ir, im] > mean_freq_sst[ib-1, ir, im]):
                        mean_freq_sst[ib, ir, im] = mean_freq_sst[ib-1, ir, im]
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
    min_ice_conc = 20.0
    for im in range(0,12):
        for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
            if ib > 0 and ib < BIN['SST']['NBINS']-2:
                for ir in range(0,1):
                    if (mean_freq_sst[ib, ir, im] < min_ice_conc):
                            mean_freq_sst[ib+4:, ir, im] = 0.0
    

    im = 0
    ir = 1
    for ib, bin in enumerate(BIN['SST']['XBINS'][:-1]):
        print 'nh, im, bin, siconc ',im, bin, mean_freq_sst[ib, ir, im]

    fh = open(filename_sst_func_siconc_month_pickle, 'wb')
    pickle.dump(mean_freq_sst, fh)
    fh.close()

def smooth_sst_under_ice(sst_daily, ice_daily, ice_mask, fsst_yrm1, fsst_yrp1):
    '''
    Under ice, the SST can be very discontinuous, this coming from the original data
    Use monthly means and interpolate to create an alternative daily dataset of tos under ice (at some point in the year)

    '''
    sst_monmn = sst_daily.aggregated_by('month', iris.analysis.MEAN)
    cube_var_l = iris.cube.CubeList()

    for m in range(0, 12):
        if m > 0 and m < 11:
            mm1 = m - 1 
            mp1 = m + 1
        elif m == 0:
            con_mon = iris.Constraint(coord_values = {'month' : lambda l : l == 12})
            sstm1_yr = iris.load_cube(fsst_yrm1)
            try:
                icc.add_month_number(sstm1_yr, 'time', name = 'month')
            except:
                pass
            sstm1_mon = sstm1_yr.extract(con_mon)
            sstm1_monmn = sstm1_mon.collapsed('time', iris.analysis.MEAN)
            mm1 = 11
            mp1 = m + 1
        elif m == 11:
            mm1 = m - 1
            con_mon = iris.Constraint(coord_values = {'month' : lambda l : l == 1})
            sstp1_yr = iris.load_cube(fsst_yrp1)
            try:
                icc.add_month_number(sstp1_yr, 'time', name = 'month')
            except:
                pass
            sstp1_mon = sstp1_yr.extract(con_mon)
            sstp1_monmn = sstp1_mon.collapsed('time', iris.analysis.MEAN)
            mp1 = 0

        con_mon = iris.Constraint(coord_values = {'month' : lambda l : l == m+1})
        cube_mon = sst_daily.extract(con_mon) # daily data
        #cube_ice_mon = cube_ice.extract(con_yr & con_mon) # daily data
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
            ice_points = np.where(ice_mask.data == 1.0)
            
            if m > 0 and m < 11:
                mon_day = w1 * sst_monmn.data[mm1, :, :] + w2 * sst_monmn.data[m, :, :] + w3 * sst_monmn.data[mp1, :, :]
            elif m == 0:
                mon_day = w1 * sstm1_monmn.data[:, :] + w2 * sst_monmn.data[m, :, :] + w3 * sst_monmn.data[mp1, :, :]
            elif m == 11:
                mon_day = w1 * sst_monmn.data[mm1, :, :] + w2 * sst_monmn.data[m, :, :] + w3 * sstp1_monmn.data[:, :]

            cube_var_mon.data[iday, ice_points[0], ice_points[1]] = mon_day[ice_points[0], ice_points[1]]

        cube_var_l.append(cube_var_mon)
    cube_var = cube_var_l.concatenate_cube()
    
    return cube_var
    

def fix_sst_under_ice(dir_in, filename_sice_func_sst_month_pickle, year, months, year_last_real_data, sst_fixed_year):

    fnames_sst = sorted(glob.glob(os.path.join(dir_in, 'tos*'+year+'0101-*.nc')))
    fnames_ice = sorted(glob.glob(os.path.join(dir_in, 'siconc*'+year+'0101-*.nc')))
    fout_year = sst_fixed_year

    if os.path.exists(fout_year):
        return

    year_min = int(os.path.basename(fnames_sst[0]).split('_')[-1][0:4])
    NYEARS = len(fnames_sst)
    print 'NYEARS ',NYEARS

    cube_mon = iris.cube.CubeList()

    # read calculated relationship
    fh = open(filename_sice_func_sst_month_pickle, 'r')
    mean_freq_sice = pickle.load(fh)
    fh.close()

    mean_freq_sice_daily = sic_functions.interpolate_histogram(year, mean_freq_sice)

    for f_sst, f_ice in zip(fnames_sst, fnames_ice):
        year = f_sst.split('_')[-1][0:4]
        year1 = f_ice.split('_')[-1][0:4]
        iy = int(year) - year_min
        print 'process ',f_sst, year, iy

        if year != year1:
            raise Exception('Paired SST and sea-ice years not same '+year+' '+year1)
        sst = iris.load_cube(f_sst)
        ice = iris.load_cube(f_ice)
        icc.add_month_number(sst, 'time', name = 'month')
        icc.add_month_number(ice, 'time', name = 'month')

        if os.path.exists(fout_year):
            continue

        files_fixed = []
        for im in months:
            print 'processing month ',im
            month = im+1
            month_no = iris.Constraint(coord_values = {'month' : lambda l : l == month})
            sst_mon = sst.extract(month_no)
            ice_mon = ice.extract(month_no)

            print ('calc the fixed SST under siconc ')
            sst_fixed = fix_sst_based_on_siconc(sst_mon, ice_mon, mean_freq_sice_daily, im)
            fout = fout_year[:-3]+'_'+str(month).zfill(2)+'.nc'
            iris.save(sst_fixed, fout, unlimited_dimensions = ['time'])
            print ('saved file ', fout)
            files_fixed.append(fout)
                                
        cmd = 'ncrcat -O '+' '.join(files_fixed)+' '+fout_year
        print 'cmd ',cmd
        os.system(cmd)
	for f in files_fixed:
	    os.remove(f)

def smooth_sst(dir_in, sst_fixed_year, sst_fixed_year_m1, sst_fixed_year_p1, year):
    # now smooth the fixed SST under sea-ice to remove discontinuities
    fice = os.path.join(dir_in, 'siconc_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_'+str(int(year))+'0101-'+str(int(year))+'1231.nc')
    ice = iris.load_cube(fice)

    sst_daily = iris.load_cube(sst_fixed_year)
    ice_daily = ice
    ice_mask = sic_functions.ice_maximum_extent(ice)

    print 'fsst ',sst_fixed_year_m1, sst_fixed_year_p1

    cube_smooth = smooth_sst_under_ice(sst_daily, ice_daily, ice_mask, sst_fixed_year_m1, sst_fixed_year_p1)
    iris.save(cube_smooth, sst_fixed_year[:-3]+'_monthly_under_ice.nc')
                
def process_sst_seaice(dir_in, year, months, year_last_real_data, sst_fixed_year, sst_ice_only = False):

    cube_mon = iris.cube.CubeList()
    # 271.32 is what HadISST2 uses under full seaice
    #ice_temp = 271.32 # 273.15 - 1.8

    filename_sice_func_sst_month_pickle = os.path.join(DATADIR, 'siconc_func_sst_month_relationship_'+year+'.pkl')
    filename_sst_func_siconc_month_pickle = os.path.join(DATADIR, 'sst_func_siconc_month_relationship_'+year+'.pkl')

    # first calculate sice func sst, to use to fix SST under sea-ice
    if not sst_ice_only:
        calc_sice_func_sst_relationship(dir_in, filename_sice_func_sst_month_pickle, year)
        
        fix_sst_under_ice(dir_in, filename_sice_func_sst_month_pickle, year, months, year_last_real_data, sst_fixed_year)

    # then calculate sst func sice using fixed SST under sea-ice
    fnames_sst = sorted(glob.glob(os.path.join(dir_in, 'tos*'+year+'0101-*.nc')))
    fnames_ice = sorted(glob.glob(os.path.join(dir_in, 'siconc*'+year+'0101-*.nc')))

    sic_functions.calc_sst_func_sice_relationship(dir_in, filename_sst_func_siconc_month_pickle, year, fnames_sst, fnames_ice, restrict_ice_edge = True)

if __name__ == '__main__':
    dir_in = '/group_workspaces/jasmin2/primavera1/WP6/forcing/HadISST2_submit/v1.2/'
    sst_fixed = os.path.join(DATADIR, 'hadisst2_tos_daily_{}_fixed_day_v1.nc')
    year_last_real_data = '2015'
    #years = range(2015, 2016)
    #years = range(1967, 2016)
    months = range(0,12)
    years = range(1967, 1972)
    #months = range(0,1)

    for year in years:
        year_str = str(year)
        sst_fixed_year = sst_fixed.format(year_str)
        process_sst_seaice(dir_in, year_str, months, year_last_real_data, sst_fixed_year, sst_ice_only = False)
    
    for year in years[1:]:
        year_str = str(year)
        sst_fixed_year = sst_fixed.format(year_str)
        sst_fixed_year_m1 = sst_fixed.format(str(year-1))
        sst_fixed_year_p1 = sst_fixed.format(str(year+1))
        if not os.path.exists(sst_fixed_year_p1):
            if str(year) == year_last_real_data:
                sst_fixed_year_p1 = sst_fixed.format(str(year))
            else:
                raise Exception('Not last year but no file found '+year+' '+sst_fixed_year_p1)
        
        smooth_sst(dir_in, sst_fixed_year, sst_fixed_year_m1, sst_fixed_year_p1, year_str)



                       
