'''
Created on Feb 15, 2021

@author: Erik Johansson
'''
import numpy as np
import netCDF4 as nc  # @UnresolvedImport
import pdb
import os
import h5py  # @UnresolvedImport
if __name__ == '__main__':
#     mainDir = '/scratch/b.legras/NWCGEO-APPLI-2018.1-9/export/CMA'
    mainDir = '/scratch/b.legras/OLD-V2016/NWCGEO-APPLI-4/export/CMA'
#     fil = 'S_NWC_CMA_MSG1_FULLAMA-VISIR_20170515T053000Z.nc'
    fil = 'S_NWC_CMA_HIMAWARI08_FULLAMA-NR_20170430T190000Z.nc'
    filename = mainDir + '/' + fil
    ncf = nc.Dataset(filename)
    
    def createLatLon(ncf, fn):
        basename = os.path.basename(fn)
        saveDir = '/scratch/erikj/Data/Nwcgeo/LatLon'
        lon = ncf['lon'][:].data
        lat = ncf['lat'][:].data
        satname = basename.split('_')[3].lower()
        sfn = '%s/%s_LonLat.h5' %(saveDir, satname)
        if not os.path.isfile(sfn):
            f = h5py.File(sfn, 'w')
            f.create_dataset('lon', data = lon)
            f.create_dataset('lat', data = lat)
            f['lon'].attrs['long_name'] = str(ncf['lon'].long_name).encode('utf-8')
            f['lon'].attrs.__setitem__('FillValue',  ncf['lon']._FillValue)
#             f['lon'].attrs.__setitem__('long_name',  ncf['lon'].long_name)
            f['lon'].attrs.__setitem__('units',  ncf['lon'].units)
            f['lat'].attrs.__setitem__('FillValue',  ncf['lat']._FillValue)
            f['lat'].attrs.__setitem__('long_name',  ncf['lat'].long_name)
            f['lat'].attrs.__setitem__('units',  ncf['lat'].units)
            f.close()
#             root_grp = nc.Dataset(sfn, 'w', format='NETCDF4')
#             root_grp.description = 'Latitude and Longitude for %s' %satname.title()
#             
            
        f5 = h5py.File(sfn, 'r')   
        pdb.set_trace()
    createLatLon(ncf, filename)
    
    
    cma = ncf['cma'][:].data.astype(float)
    cma[ncf['cma'][:].mask] = np.nan
    
    
    pdb.set_trace()
    