'''
Created on Feb 15, 2021

@author: Erik Johansson
'''
import numpy as np
import netCDF4 as nc  # @UnresolvedImport
import pdb
import os
import h5py  # @UnresolvedImport
import glob


try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen
from io import BytesIO

import cartopy.crs as ccrs
import matplotlib.pyplot as plt



#----------------------------------------------------

def saveLatLon(ncf, fn):
    basename = os.path.basename(fn)
    LatLonDir = '/scratch/erikj/Data/Nwcgeo/LatLon'
    satname = basename.split('_')[3].lower()
    sfn = '%s/%s_LonLat.h5' %(LatLonDir, satname)
    if not os.path.isfile(sfn):
        lon = ncf['lon'][:].data
        lat = ncf['lat'][:].data
        f = h5py.File(sfn, 'w')
        f.create_dataset('lon', data = lon)
        f.create_dataset('lat', data = lat)
#             f['lon'].attrs['long_name'] = str(ncf['lon'].long_name).encode('utf-8')
        f['lon'].attrs.__setitem__('FillValue',  ncf['lon']._FillValue)
        f['lon'].attrs.__setitem__('long_name',  ncf['lon'].long_name)
        f['lon'].attrs.__setitem__('units',  ncf['lon'].units)
        f['lat'].attrs.__setitem__('FillValue',  ncf['lat']._FillValue)
        f['lat'].attrs.__setitem__('long_name',  ncf['lat'].long_name)
        f['lat'].attrs.__setitem__('units',  ncf['lat'].units)
        f.close()
    return sfn


def getLatLon(gfn):
#     lon = ncf['lon'][:].data
#     lat = ncf['lat'][:].data
    f = h5py.File(gfn, 'r')
    lat = np.asarray(f['lat'])
    lon = np.asarray(f['lon'])
    lat = np.where(lat == f['lat'].attrs['FillValue'], np.nan, lat)
    lon = np.where(lon == f['lon'].attrs['FillValue'], np.nan, lon)
    f.close()
    return lat, lon


#----------------------------------------------------

if __name__ == '__main__':
    mainDir_new = '/scratch/b.legras/NWCGEO-APPLI-2018.1-9/export'
    mainDir_old = '/scratch/b.legras/OLD-V2016/NWCGEO-APPLI-3/export'
    plotDir = '/scratch/erikj/Data/Nwcgeo/Compare/Plots'
#     fil = 'S_NWC_CMA_MSG1_FULLAMA-VISIR_20170515T053000Z.nc'
    
    sat='MSG1'
    datum = '20170429T1400'
    product = 'CTTH'
#     HIMAWARI08
    files_new = glob.glob('%s/%s/S_NWC_%s_%s_FULLAMA-*_%s*.nc' %(mainDir_new, product, product, sat, datum))
    files_old = glob.glob('%s/%s/S_NWC_%s_%s_FULLAMA-*_%s*.nc' %(mainDir_old, product, product, sat, datum))
    files_new.sort()
    files_new.sort()
    fil_new = files_new[0]
    fil_old = files_old[0]
#     filename = mainDir + '/' + fil
    ncf_new = nc.Dataset(fil_new)
    ncf_old = nc.Dataset(fil_old)
    #: Read the required data
    nc_prod = product.lower()
    if nc_prod == 'ctth':
        nc_prod = 'ctth_pres'
    prod_new = ncf_new[nc_prod][:].data.astype(float)
    prod_old = ncf_old[nc_prod][:].data.astype(float)
    
    prod_new[ncf_new[nc_prod][:].mask] = np.nan
    prod_old[ncf_old[nc_prod][:].mask] = np.nan
    
    #: Take care of latitude and longitude
    sfn = saveLatLon(ncf_old, fil_old)
    ncf_old.close()
    ncf_new.close()
    lat, lon = getLatLon(sfn)
    pdb.set_trace()
    #: Plotting    
    fig = plt.figure()
    ax = fig.add_subplot(3,1,1)#, projection=ccrs.Miller(central_longitude=180.0))#projection=ccrs.PlateCarree())
#     ax.coastlines()
    im = ax.imshow(prod_new)#, extent=(np.nanmin(lon), np.nanmin(lon), np.nanmin(lat), np.nanmin(lat)))#transform=ccrs.Geostationary(satellite_height=35786000))
    cbar = fig.colorbar(im)
    ax.set_title('2018')
    
    ax = fig.add_subplot(3,1,2)#, projection=ccrs.Miller(central_longitude=180.0))#projection=ccrs.PlateCarree())
    im = ax.imshow(prod_old)
    cbar = fig.colorbar(im)
    ax.set_title('2016')
    
    ax = fig.add_subplot(3,1,3)#, projection=ccrs.Miller(central_longitude=180.0))#projection=ccrs.PlateCarree())
    im = ax.imshow((prod_new - prod_old), cmap='RdBu_r')
    cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)  # @UnusedVariable
    ax.set_title('2018 - 2016')
    
    plt.tight_layout()
    
    figname = '%s/%s_%s' %(plotDir, nc_prod, datum)
    fig.savefig(figname + '.png')
    fig.show()
    
    
    pdb.set_trace()
    
    
    
    
    
    
    
    
#     im = ax.imshow(rhdLat, origin='lower', cmap='RdBu_r', aspect=aspect, vmin=valminmax * -1, vmax=valminmax)
#     ax.set_yticks(yticks['lre'])
#     ax.set_yticklabels(['-30', '-15', '0', '15', '30'])
# 
#     ax.set_xticks(xticks['lre'])
#     ax.set_title('High - Low')
#     ax.set_ylabel('Latitude [deg]')
# 
#     if use_datline_center:
#         ax.set_xticklabels(['0', '45', '90', '135', '180', '225', '270', '315', '360'])
#     else:
#         ax.set_xticklabels(['-180', '-135', '-90', '-45', '0', '45', '90', '135', '180'])
#     ax.set_xlabel('Longitude [deg]')
#     barticks = [valminmax*-1, valminmax*-0.75, valminmax*-0.5, valminmax*-0.25, 0, valminmax*0.25, valminmax*0.5, valminmax*0.75, valminmax]
#     cbar_ax = fig.add_axes([0.2, 0.25, 0.6, 0.01])
#     cbar = fig.colorbar(im, orientation='horizontal', cax=cbar_ax, ticks=barticks)  # @UnusedVariable
#         else:
#             ax.set_xticklabels(['', '', '', '', '', '', '', '', ''])
#     figname = figname_st + '_top'
#     if not use_datline_center:
#         figname = os.path.dirname(figname) + '/map_' + os.path.basename(figname)
#     if useClim:
#         figname = figname + '_anom'
#     fig.savefig(figname + '.png')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(1,1,1, projection=ccrs.Miller(central_longitude=180.0))#projection=ccrs.PlateCarree())
    ax1.coastlines()
    
    ax1.imshow(cma, extent=(np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)))#transform=ccrs.Geostationary(satellite_height=35786000))
    fig1.show()
    pdb.set_trace()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.Miller())
    ax.coastlines()
    ax.set_global()
#     img, crs, extent, origin = geos_image()
    img_proj = ccrs.Geostationary(satellite_height=35786000)
    pdb.set_trace()
    plt.imshow(cma, transform=img_proj, origin='lower', cmap='gray')#, extent=extent, origin=origin, cmap='gray')
    pdb.set_trace()
    plt.show()
    
    pdb.set_trace()
    