import os
import iris
import glob
import iris.coord_categorisation as icc

def make_hadisst_1x1degree(dir_in, dir_out):
    '''
   Filling in mask over land
   Horizontal remapping to grid specified in interp_setup
   (within metoffice_conversion) - Calendar adjustment (360 day for Met Office) and pp format
    '''

    for yr in YEARS:
        for variable in interp_setup['variables']:
            file = LATLONG_YR.format(dir_in, interp_setup['variables'][variable], str(yr))
            #cube = iris.load_cube(file)
            # is it better to fill in land before remapping, or remap first and then fill?
            #cube_filled = calculate_filled(cube, file[:-3]+'_filled.nc')
        
            # remapping to model grid
            file_regrid = REGRID_1DEGREE_YR.format(dir_out, interp_setup['variables'][variable], str(yr))
            #cube_remap = remap_cube_with_iris(cube_filled, interp_setup, variable)
            remap_with_cdo(file, file_regrid, '1x1')

def make_hadisst_1x1degree(dir_in, dir_out):
    '''
   Filling in mask over land
   Horizontal remapping to grid specified in interp_setup
   (within metoffice_conversion) - Calendar adjustment (360 day for Met Office) and pp format
    '''

    for yr in YEARS:
        for variable in interp_setup['variables']:
            file = LATLONG_YR.format(dir_in, interp_setup['variables'][variable], str(yr))
            #cube = iris.load_cube(file)
            # is it better to fill in land before remapping, or remap first and then fill?
            #cube_filled = calculate_filled(cube, file[:-3]+'_filled.nc')
        
            # remapping to model grid
            file_regrid = REGRID_1DEGREE_YR.format(dir_out, interp_setup['variables'][variable], str(yr))
            #cube_remap = remap_cube_with_iris(cube_filled, interp_setup, variable)
            remap_with_cdo(file, file_regrid, '1x1')


def remap_with_cdo(file_in, file_out, grid_out):
    namelist_file = make_cdo_namelist(grid_out)
    cmd = 'cdo remapbil,'+namelist_file+' '+file_in+' '+file_out
    os.system(cmd)

def make_cdo_namelist(grid_out):
    namelist_file = os.path.join(os.getcwd(), 'namelist_grid1x1')
    with open(namelist_file,'w') as outp:
        outp.write('gridtype = lonlat \n')
        if grid_out == '1x1':
            outp.write('xsize = 360 \n')
            outp.write('ysize = 180 \n')
            outp.write('xfirst = 0.0 \n')
            outp.write('yfirst = -89.5 \n')
            outp.write('yinc = 1 \n')
        elif grid_out == '1/4x1/4':
            outp.write('xsize = 1440 \n')
            outp.write('ysize = 720 \n')
            outp.write('xfirst = 0.125 \n')
            outp.write('yfirst = -89.875 \n')
            outp.write('yinc = 0.25 \n')
            outp.write('\n')
        else:
            raise ValueError('Unknown grid for output '+grid_out)
    return namelist_file

if __name__ == '__main__':
    dir_in = '/home/users/mjrobert/hrcm/cache/malcolm/HighResMIP/sst_forcing/processing_ocean/'
    cmip_file = os.path.join(dir_in, 'cmip5_mean_trend_filtered_modelmean_masked_time.nc')
    datadir = '/home/users/mjrobert/hrcm/cache/malcolm/HadISST2/1x1/processing_2018'    
    datafile = os.path.join(datadir, 'hadisst2_monthly_1948-2015_tos_01-12_1x1.nc')

    hadisst2_dir = '/home/users/mjrobert/hrcm/cache/malcolm/HadISST2/1x1/'
    cmip_1x1 = os.path.join(hadisst2_dir, 'cmip5_trend_1x1.nc')

    if not os.path.exists(cmip_1x1):
        trend_025 = iris.load_cube(cmip_file)
        trend_025.coord('longitude').circular = True
        c_ref = iris.load(datafile)[0]
        for coord in ['latitude','longitude']:
            c_ref.coord(coord).guess_bounds()
            trend_025.coord(coord).guess_bounds()
        trend_1x1 = trend_025.regrid(c_ref, iris.analysis.AreaWeighted())
        iris.save(trend_1x1, cmip_1x1)

    trend = iris.load_cube(cmip_1x1)
    icc.add_day_of_year(trend, 'time')
    icc.add_year(trend, 'time')

    hadisst2_files = glob.glob(os.path.join(hadisst2_dir, 'HadISST2_1x1_regrid_sst*'))
    

    '''
    for each year, read in the HadISST2 daily data
    read in the monthly data, dec year-1 to jan year +1
    do time interpolation from one time to the other
    '''
    hadisst_data = iris.load_cube(os.path.join(hadisst2_dir, 'HadISST2_1x1_regrid_sst_2002.nc'))
    icc.add_day_of_year(hadisst_data, 'time')
    YEARS = iris.Constraint(coord_values = {'year' : lambda l : l == 2002})
    sub_trend = trend.extract(YEARS)
    print  sub_trend
    try:
        icc.add_day_of_year(sub_trend, 'time')
    except:
        pass
    points_to_calc = [('day_of_year', hadisst_data.coord('day_of_year').points)]
    print points_to_calc
    print sub_trend.coord('day_of_year').points
    sub_trend_daily = sub_trend.interpolate(points_to_calc, iris.analysis.Linear())
    iris.save(sub_trend_daily, os.path.join(hadisst2_dir, 'test.nc'))

