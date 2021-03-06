import matplotlib
matplotlib.use('Agg')
import iris, os, sys, glob
import iris.coord_categorisation as icc
import iris.analysis
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
import iris.quickplot as qplt
import matplotlib.dates as mdates

dir_in = '/gws/nopw/j04/hrcm/cache/malcolm/HadISST2/1x1/processing_2018/future2/'

dir_in_hist = '/group_workspaces/jasmin2/primavera1/WP6/forcing/HadISST2_submit/v1.2/'

fice = 'future_ice_conc.nc'

npole_lat = iris.Constraint(coord_values = {'latitude' : lambda l : 60 < (l.point) <= 90})
npole_lon = iris.Constraint(coord_values = {'longitude' : lambda l : 0 < (l.point) <= 360})
spole_lat = iris.Constraint(coord_values = {'latitude' : lambda l : -90 < (l.point) <= -50})
spole_lon = iris.Constraint(coord_values = {'longitude' : lambda l : 0 < (l.point) <= 360})

def history_period_callback(cube, field, filename):
    # remove attributes preventing cube concatenation
    attrib = ['history','time_period']
    for att in attrib:
        try:
            del cube.attributes[att]
        except:
            pass

    coords = ['year', 'month', 'day_of_year']

    for coord in coords:
        try:
            cube.remove_coord(coord)
        except:
            pass

    cube.coord('time').bounds = None
    cube.long_name = cube.var_name

def guess_areas(cube):
    coords = ['latitude','longitude']
    for coord in coords:
        cube.coord(coord).guess_bounds()

    grid_areas = iris.analysis.cartography.area_weights(cube)
    return grid_areas

def adjust_sst_under_ice(cube, fname_hist_sst, fname_hist_sic):
    '''
    Change the SST under ice to be the same as in the historic, to check that we do indeed get the same global SST timeseries
    '''
    hist_sic = iris.load_cube(fname_hist_sic)
    ice = np.where(hist_sic.data > 0.0)
    del hist_sic

    hist_sst = iris.load_cube(fname_hist_sst)
    cube.data[ice[0], ice[1], ice[2]] = hist_sst.data[ice[0], ice[1], ice[2]]

    return cube

def future_sst(mon, years, pole = 'global'):
    if pole == 'npole':
        pole_lat = npole_lat
        pole_lon = npole_lat
    else:
        pole_lat = spole_lat
        pole_lon = spole_lat

    files = []
    sst_l = iris.cube.CubeList()
    for year in years:
        if year < 2016:
           dirin = os.path.join(dir_in, '../')
           fname = os.path.join(dirin, 'hadisst2_tos_daily_'+str(year)+'_fixed_day_v1_monthly_under_ice.nc')
        else:
           dirin = os.path.join(dir_in, 'sst')
           fname = os.path.join(dirin, 'future_sst_'+str(year)+'_025_daily_v1.nc')

        if os.path.exists(fname):
            cube = iris.load_cube(fname, callback = history_period_callback)
            #if year < 2016:
            #    cube = adjust_sst_under_ice(cube, fname_hist_sst, fname_hist_sic)
            sst_l.append(cube)
    print 'files ',files

    #sst_l = iris.load(files, callback = history_period_callback)
    print 'sst_1 ',sst_l
    print sst_l[0]
    print sst_l[-1]
    print sst_l[0].coord('time')
    print sst_l[-1].coord('time')
    for ic, c in enumerate(sst_l):
        print 'ic ',ic, iris.util.describe_diff(sst_l[0], c)
    sst = sst_l.concatenate_cube()
    print 'sst ',sst
        
    try:
        icc.add_year(sst, 'time', name = 'year')
        icc.add_month_number(sst, 'time', name = 'month')
    except:
        pass

    con_mon = iris.Constraint(coord_values = {'month' : lambda l : l == mon})
    #sst_mon = sst.extract(con_mon)
    sst_mon = sst
    del sst, sst_l
    print 'sst_mon ',sst_mon
    print 'sst_mon years ',sst_mon.coord('year').points

    sst_yr = sst_mon.aggregated_by(['month','year'], iris.analysis.MEAN)
    grid_areas = guess_areas(sst_yr)

    sst_ts_global = sst_yr.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights = grid_areas)

    sst_ts = sst_yr.collapsed(['longitude'], iris.analysis.MEAN)
    print sst_ts

    sst_ts_anom = sst_ts.copy()

    sst_ts_anom.data[:, :] -= sst_ts.data[0,:]

    #iplt.plot(ice_ts)
    return sst_ts, sst_ts_anom, sst_ts_global

def historic_sst(mon, years, pole = 'global'):
    if pole == 'npole':
        pole_lat = npole_lat
        pole_lon = npole_lat
    else:
        pole_lat = spole_lat
        pole_lon = spole_lat

    files = []
    for year in years:
        fname = os.path.join(dir_in_hist, 'tos_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_'+str(year)+'0101-'+str(year)+'1231.nc')
        print fname
        if os.path.exists(fname):
            files.append(fname)
    print files

    sst_l = iris.load(files, callback = history_period_callback)
    sst = sst_l.concatenate_cube()
    print 'sst ',sst
        
    try:
        icc.add_year(sst, 'time', name = 'year')
        icc.add_month_number(sst, 'time', name = 'month')
    except:
        pass

    print 'sst with coords ',sst

    con_mon = iris.Constraint(coord_values = {'month' : lambda l : l == mon})
    #sst_mon = sst.extract(con_mon)
    sst_mon = sst
    print 'sst_mon ',sst_mon

    sst_yr = sst_mon.aggregated_by(['month','year'], iris.analysis.MEAN)
    grid_areas = guess_areas(sst_yr)

    sst_ts_global = sst_yr.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights = grid_areas)

    sst_ts = sst_yr.collapsed(['longitude'], iris.analysis.MEAN)
    print sst_ts
    sst_ts_anom = sst_ts.copy()

    sst_ts_anom.data[:, :] -= sst_ts.data[0,:]

    del sst, sst_l

    #iplt.plot(ice_ts)
    return sst_ts, sst_ts_anom, sst_ts_global

def plot_point_series():
    lat = iris.Constraint(coord_values = {'latitude' : lambda l : -5 < (l.point) <= 5})
    lon = iris.Constraint(coord_values = {'longitude' : lambda l : 190 < (l.point) <= 240})

    years_hist = range(1972, 2016)
    files_hist = []
    dirin = '/gws/nopw/j04/hrcm/cache/malcolm/HadISST2/1x1/processing_2018/'
    for year in years_hist:
        f = os.path.join(dirin, 'hadisst2_tos_daily_'+str(year)+'_fixed_day_v1_monthly_under_ice.nc')
        #f = os.path.join(dir_in_hist, 'tos_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_'+str(year)+'0101-'+str(year)+'1231.nc')
        files_hist.append(f)

    chl = iris.load(files_hist, lon&lat, callback = history_period_callback)
    ch = chl.concatenate_cube()
    chm = ch.collapsed(['latitude','longitude'], iris.analysis.MEAN)
    #print 'chm', chm
    #iris.save(chm, './hist_sst_ts.nc')

    years_fut = range(2016, 2050)
    sst_future_dir = os.path.join(dir_in, 'sst/')

    sst_file = 'future_sst_{}_025_daily_v1.nc'

    files = []
    for year in years_fut:
        f = os.path.join(sst_future_dir, sst_file.format(str(year)))
        files.append(f)

    cl = iris.load(files, lon&lat)
    c = cl.concatenate_cube()
    cm = c.collapsed(['latitude','longitude'], iris.analysis.MEAN)
    print 'cm', cm
    iris.save(cm, './future_sst_ts.nc')

    iplt.plot(cm)
    iplt.plot(chm)
    print c.coord('year').points
    plt.show()
        
    iplt.plot(chm)
    plt.savefig('./hist_future2_sst_region_ts_fixed.png')
    plt.show()

if __name__ == '__main__':

    yearstart = 1952
    years_fut = range(1972, 2059)
    years_hist = range(yearstart, 2016)

    #yearstart = 2010
    years_fut = range(2000, 2035)
    #years_hist = range(yearstart, 2016)

    #plot_point_series()
    #stop

    for pole in ['global']:
        for mon in range(1,13):
            sst_hist, sst_hist_anom, sst_hist_ts = historic_sst(mon, years_hist, pole = pole)
            sst_fut, sst_fut_anom, sst_fut_ts = future_sst(mon, years_fut, pole = pole)
            #plt.plot(sst_hist.coord('year').points, sst_hist.data)
            #plt.plot(sst_fut.coord('year').points, sst_fut.data)
            levels = np.arange(-1, 2.5, 0.05) 
            fig = plt.figure()
            ax = fig.add_subplot(2,1,1)
            qplt.contourf(sst_hist_anom, levels)

    # Stop matplotlib providing clever axes range padding
            plt.axis('tight')
    # As we are plotting annual variability, put years as the x ticks
            plt.gca().xaxis.set_major_locator(mdates.YearLocator(10))

    # And format the ticks to just show the year
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

            ax = fig.add_subplot(2,1,2)
            qplt.contourf(sst_fut_anom, levels)

    # Stop matplotlib providing clever axes range padding
            plt.axis('tight')
    # As we are plotting annual variability, put years as the y ticks
            plt.gca().xaxis.set_major_locator(mdates.YearLocator(10))

    # And format the ticks to just show the year
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

            plt.title('Month '+str(mon)+', '+pole)
            plt.savefig('./hist_future3_sst_'+str(mon)+'_'+pole+'_fixed.png')

            fig = plt.figure()
            plt.plot(sst_hist_ts.coord('year').points, sst_hist_ts.data, label = 'Historic')
            plt.plot(sst_fut_ts.coord('year').points, sst_fut_ts.data, label = 'Future')
            plt.title('Month '+str(mon))
            plt.legend(loc = 'best')
            plt.savefig('./hist_future3_sst_global_ts_'+str(mon)+'_'+pole+'_fixed.png')

        #plt.show()

