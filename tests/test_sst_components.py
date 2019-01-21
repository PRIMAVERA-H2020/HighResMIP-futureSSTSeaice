import iris
import numpy as np
import matplotlib.pyplot as plt
import iris.coord_categorisation as icc

dir_in_hadisst = '/home/users/mjrobert/hrcm/cache/malcolm/HadISST2/1x1/processing_2018/future/sst/'
fname_had = 'hadisst2_running_monthly_mean_025_fixed_v1.nc'
fname_format = 'sst_variability_{}_025_{}_fixed_v1.nc'
fname_format1 = 'future_sst_{}_025_{}_v1.nc'
fname_format = '/../../hadisst2_tos_daily_{}_fixed_day_v1_monthly_under_ice.nc'

lat = iris.Constraint(coord_values = {'latitude' : lambda l : 80 < (l.point) <= 85})
lat = iris.Constraint(coord_values = {'latitude' : lambda l : -5 < (l.point) <= 5})
#lat = iris.Constraint(coord_values = {'latitude' : lambda l : 25 < (l.point) <= 35})
#lat = iris.Constraint(coord_values = {'latitude' : lambda l : -25 < (l.point) <= 25})
lon = iris.Constraint(coord_values = {'longitude' : lambda l : 190 < (l.point) <= 240})
#lon = iris.Constraint(coord_values = {'longitude' : lambda l : 0 < (l.point) <= 360})
con_mon = iris.Constraint(coord_values = {'month' : lambda l : l == 7})
con_year = iris.Constraint(coord_values = {'year' : lambda l : 2000 <= l <= 2015})

mask_ref = iris.load_cube(dir_in_hadisst + fname_format1.format('2016', 'daily'))[0]
years = range(1976, 1980)
years = range(2014, 2018)
fn = []; yr = []
offset = 0
for year in years:
    if year > 2015:
        fname = fname_format1.format(str(year), 'daily')
    else:
        fname = fname_format.format(str(year))
    c = iris.load_cube(dir_in_hadisst + fname)
    for coord in ['latitude', 'longitude']:
        c.coord(coord).guess_bounds()
                          
#icc.add_month_number(cf, 'time', name = 'month')
    ce = c
    #ce = ce.extract(lat&lon)
    #ce = ce.extract(con_mon)
                          
    grid = iris.analysis.cartography.area_weights(ce)
    ce.data.mask[:, :, :] = mask_ref.data.mask[:,:]
    cm = ce.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights = grid)
    print year
    days = range(cm.shape[0])
    d = [x + offset for x in days]
    plt.plot(d, cm.data)
    offset += len(days)
plt.savefig('./sst_'+str(year)+'.png')
plt.show()

