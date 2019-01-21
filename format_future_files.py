import os, sys
import iris
import numpy as np

DATADIR = '/home/users/mjrobert/hrcm/cache/malcolm/HadISST2/1x1/processing_2018'
HIST_DIR = '/group_workspaces/jasmin2/primavera1/WP6/forcing/HadISST2_submit/v1.2'
FUTURE_DIR = '/group_workspaces/jasmin2/primavera1/WP6/forcing/FutureSSTandSeaice/c.1/'
if not os.path.exists(FUTURE_DIR):
    os.makedirs(FUTURE_DIR)

hist_cmip_file = os.path.join(HIST_DIR, '{}_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_{}0101-{}1231.nc')
future_cmip_file = os.path.join(FUTURE_DIR, '{}_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-highresSST-future_gn_{}0101-{}1231.nc')

def metadata(cf_name, variable_name, variable_longname, time_period, realm):
    '''
    Metadata from PCMDI monthly file
    '''

    metadata = {
    'version': '20190120',
    'activity_id': 'input4MIPs',
    'comment': '',
    'contact': 'Met Office (malcolm.roberts@metoffice.gov.uk)',
    'creation_date': '2019-01-20T12:00:00Z',
    'data_structure': 'grid',
    'dataset_category': 'SSTsAndSeaIce',
    'dataset_version_number': 'r0',
    'frequency': 'day',
    'further_info_url': 'http://collab.knmi.nl/highresmip',
    'grid': "0.25x0.25 degree latitude x longitude",
    'grid_label': 'gn',
    'history': 'File processed on 2019-01-20, converting original netCDF3 files to deflated CF-compliant netCDF4 with additional metadata for CMIP6',
    'institution': 'Met Office Hadley Centre, Fitzroy Road, Exeter, Devon, EX1 3PB, UK',
    'institution_id': 'MOHC',
    'mip_era': 'CMIP6',
    'nominal_resolution': '25 km',
    'product': 'forcing',
    'realm': realm,
    'references': 'Haarsma et al, 2017: High Resolution Model Intercomparison Project (HighResMIP).',
    'source': 'Processed from combining past HadISST2.2 observations with future SST from CMIP5 coupled historic and rcp85 simulations ',
    'source_id': 'highresSST-future-r0',
    'supplementary_information': '',
    'table_id': 'input4MIPs',
    'target_mip': 'HighResMIP',
    'time_period': time_period, 
    'title': 'highresSST-future-0.25x0.25 dataset prepared by Met Office for CMIP6 HighResMIP',
    'variable_id': variable_name, 
    'license': 'highresSST-future dataset is created only for CMIP6 HighResMIP, it is not a projection of future climate and should not be used as such. The material may be downloaded for the purposes of private study and scientific research. Any other proposed use of the material is subject to a copyright licence available from the Met Office. Licences and further information can be obtained from the Met Office IPR Officer, Met Office, FitzRoy Road, Exeter, Devon, EX1 3PB. E-mail: ipr@metoffice.gov.uk.'
    }
    return metadata

def add_metadata(cube, var, period, variable_names, realm):
    '''
    Add netcdf metadata to file
    Change from HadISST2.2 to highresSST-future in metadata
    contact, creation_date, dataset_version_number, license, product, references, source, source_id, title, version
    '''
    if var == 'tos':
        cube.var_name = 'tos'
        cube.long_name = 'HighResMIP Future Sea Surface Temperature'
        cube.standard_name = 'sea_surface_temperature'
    else:
        cube.var_name = 'siconc'
        cube.long_name = 'HighResMIP Future Sea Ice Concentration'
        cube.standard_name = 'sea_ice_area_fraction'
        cube.units = '%'


    pcmdi_future_metadata = metadata(cube.standard_name, variable_names['variable_name'], variable_names['long_name'], period, realm)

    for meta in pcmdi_future_metadata:
        cube.attributes[str(meta)] = pcmdi_future_metadata[meta]
    #cube.coord('time').guess_bounds()

    extra_coords = ['year', 'month', 'day_of_year']
    for coord in extra_coords:
        try:
            cube.remove_coord(coord)
        except:
            pass
    cube.coord('time').bounds = None

    extra_attribs = ['nco_openmp_thread_number']
    for attrib in extra_attribs:
        try:
            del cube.attributes[attrib]
        except:
            pass

    return cube

def add_common_mask(cube):
    '''
    Make the mask the same between tos and sic
    '''
    sst_hist = hist_cmip_file.format('tos', '1950', '1950')
    sic_hist = hist_cmip_file.format('siconc', '1950', '1950')
    sst_mask = iris.load_cube(sst_hist)[0]
    sic_mask = iris.load_cube(sic_hist)[0]

    mask1 = cube[0].data.mask
    miss = np.where(sst_mask.data.mask == True)
    miss1 = np.where(sic_mask.data.mask == True)

    mask1[miss[0], miss[1]] = True
    mask1[miss1[0], miss1[1]] = True

    for t in range(cube.shape[0]):
        cube.data.mask[t, :, :] = mask1[:, :]

    return cube

def work():
    '''
    Process tos and siconc future files to a more reasonable name and format
    '''
    sea_ice = {'cf_standard_name': 'sea_ice_area_fraction', 'variable_name':'siconc', 'long_name': 'HighResMIP Future Sea Ice Concentration', 'realm': 'seaIce'}
    sst = {'cf_standard_name': 'sea_surface_temperature', 'variable_name':'tos', 'long_name': 'HighResMIP Future Sea Surface Temperature', 'realm': 'ocean'}
    variables = [sst, sea_ice]

    year = 2015
    for var in ['tos' ,'siconc']:
        if var == 'tos':
            future_file_orig = os.path.join(DATADIR, 'hadisst2_tos_daily_'+str(year)+'_fixed.nc')
            variable_names = variables[0]
        else:
            future_file_orig = os.path.join(HIST_DIR, 'siconc_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-HadISST-2-2-0-0-0_gn_'+str(year)+'0101-'+str(year)+'1231.nc')
            variable_names = variables[1]
        future_file_final = future_cmip_file.format(var, str(year), str(year))
        if os.path.exists(future_file_orig):
            print future_file_orig
            cube = iris.load_cube(future_file_orig)
            period = str(year)+'0101-'+str(year)+'1231'
            add_metadata(cube, var, period, variable_names, variable_names['realm'])
            add_common_mask(cube)
            iris.save(cube, future_file_final, unlimited_dimensions = ['time'], fill_value = 1.0e20)

    years = range(2016, 2051)
    #years = range(2016, 2017)
    for var in ['tos' ,'siconc']:
        for year in years:
            if var == 'tos':
                future_file_orig = os.path.join(DATADIR, 'future2', 'sst', 'future_sst_'+str(year)+'_025_daily_v1.nc')
                variable_names = variables[0]
            else:
                future_file_orig = os.path.join(DATADIR, 'future2', 'siconc', 'future_siconc_'+str(year)+'_025_daily_v1.nc')
                variable_names = variables[1]
            future_file_final = future_cmip_file.format(var, str(year), str(year))
            if os.path.exists(future_file_orig):
                print future_file_orig
                cube = iris.load_cube(future_file_orig)
                period = str(year)+'0101-'+str(year)+'1231'
                add_metadata(cube, var, period, variable_names, variable_names['realm'])
                add_common_mask(cube)
                iris.save(cube, future_file_final, unlimited_dimensions = ['time'], fill_value = 1.0e20)


if __name__ == '__main__':
    work()
