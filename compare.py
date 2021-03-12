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
import sys


# try:
#     from urllib2 import urlopen
# except ImportError:
#     from urllib.request import urlopen
# from io import BytesIO
import cartopy.crs as ccrs  # @UnresolvedImport
from cartopy import feature  # @UnresolvedImport
import matplotlib.pyplot as plt  # @UnresolvedImport

sys.path.append('/home/erikj/Projects/STC/pylib')
import geosat
from SAFNWCnc import SAFNWC_CTTH
from datetime import datetime

satConversion = {'himawari': 'hima08', 'hima08': 'hima08', 'msg1': 'msg1', 'msg3': 'msg3'}

#----------------------------------------------------

def chart(obt,field,cmap='jet',clim=[190.,300.],txt='',subgrid=None, block=True, xlocs=None, figsize= None, show=True):
        # test existence of key field
        if field not in obt.var.keys():
            print ('undefined field')
            return
        if subgrid == None:
            geogrid = obt.geogrid
        else:
            geogrid = subgrid
        if 'FullAMA' in geogrid.gridtype:
            fig = plt.figure(figsize=[10, 6])
        elif figsize is not None:
            fig = plt.figure(figsize=figsize)
        else:
            fig = plt.figure(figsize=[11,4])
#         fig.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
        fs = 15
        # it is unclear how the trick with cm_lon works in imshow but it does
        # the web says that it is tricky to plot data accross dateline with cartopy
        # check https://stackoverflow.com/questions/47335851/issue-w-image-crossing-dateline-in-imshow-cartopy
        cm_lon =0
        # guess that we want to plot accross dateline
        if geogrid.box_range[0,1]> 181:
            cm_lon = 180
        proj = ccrs.PlateCarree(central_longitude=cm_lon)
#         ax = plt.axes(projection = proj)
        ax = fig.add_subplot(111, projection = proj)
        if subgrid == None:
            plotted_field = obt.var[field]
        else:
            # extraction in subgrid
            plotted_field = obt.var[field][geogrid.corner[1]:geogrid.corner[1]+geogrid.box_biny,
                                            geogrid.corner[0]:geogrid.corner[0]+geogrid.box_binx]
        iax = ax.imshow(plotted_field, transform=proj, interpolation='nearest',
                    extent=geogrid.box_range.flatten()-np.array([cm_lon,cm_lon,0,0]),
                    origin='lower', aspect=1.,cmap=cmap,clim=clim)
        ax.add_feature(feature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none'))
        ax.coastlines('50m')
        #ax.add_feature(feature.BORDERS)
        # The grid adjusts automatically with the following lines
        # If crossing the dateline, superimposition of labels there
        # can be suppressed by specifying xlocs

        if (cm_lon == 180) & (xlocs == None): xlocs = [0,30,60,90,120,150,180,-150,-120,-90,-60,-30]
        gl = ax.gridlines(draw_labels=True, xlocs=xlocs,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        #gl.xformatter = LONGITUDE_FORMATTER
        #gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': fs}
        gl.ylabel_style = {'size': fs}
        #gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
        plt.title(txt,fontsize=fs)
        # plot adjusted colorbar and show
        axpos = ax.get_position()
        pos_x = axpos.x0 + axpos.x0 + axpos.width + 0.01
        pos_cax = fig.add_axes([pos_x,axpos.y0,0.04,axpos.height])
        cbar=fig.colorbar(iax,cax=pos_cax)
        cbar.ax.tick_params(labelsize=fs)
        if show:
            fig.show()
        return ax




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
    mainDir_old = '/scratch/b.legras/sats/msg1/safnwc-SAFBox/netcdf'
#     mainDir_old = '/scratch/b.legras/OLD-V2016/NWCGEO-APPLI-3/export'
    plotDir = '/scratch/erikj/Data/Nwcgeo/Compare/Plots'
#     fil = 'S_NWC_CMA_MSG1_FULLAMA-VISIR_20170515T053000Z.nc'
    
    sat='MSG1'
    year = 2017
    mon = 5
    day = 5
    hour = 14
    min = 0
    datum = '20170429T1400'
    product = 'CTTH'
#     HIMAWARI08
    files_new = glob.glob('%s/%s/S_NWC_%s_%s_FULLAMA-*_%d%02d%02dT%02d%02d*.nc' %(mainDir_new, product, product, sat, year, mon, day, hour, min))
#     files_old = glob.glob('%s/%s/S_NWC_%s_%s_FULLAMA-*_%s*.nc' %(mainDir_old, product, product, sat, datum))
    files_old = glob.glob('%s/%d_%02d_%02d/S_NWC_%s_%s_FULLAMA-*_%d%02d%02dT%02d%02d*.nc' %(mainDir_old, year, mon, day, product, sat, year, mon, day, hour, min))
    files_new.sort()
    files_old.sort()
    fil_new = files_new[0]
    fil_old = files_old[0]
#     filename = mainDir + '/' + fil
    ncf_new = nc.Dataset(fil_new)
    ncf_old = nc.Dataset(fil_old)
    #: Read the required data
    nc_prod = product.lower()
    if nc_prod == 'ctth':
        nc_prod = 'ctth_pres'
    mask_prod_new = ncf_new[nc_prod][:]
    mask_prod_old = ncf_old[nc_prod][:]
    prod_new = mask_prod_new.data.astype(float)
    prod_old = mask_prod_old.data.astype(float)
    
    prod_new[ncf_new[nc_prod][:].mask] = np.nan
    prod_old[ncf_old[nc_prod][:].mask] = np.nan
    
    #: Take care of latitude and longitude
    sfn = saveLatLon(ncf_old, fil_old)
    ncf_old.close()
    ncf_new.close()
    lat, lon = getLatLon(sfn)
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
    prod_diff = prod_new - prod_old
    ax = fig.add_subplot(3,1,3)#, projection=ccrs.Miller(central_longitude=180.0))#projection=ccrs.PlateCarree())
    im = ax.imshow((prod_diff), cmap='RdBu_r')
    cbar = fig.colorbar(im)#, orientation='horizontal', cax=cbar_ax, ticks=barticks)  # @UnusedVariable
    ax.set_title('2018 - 2016')
    
    plt.tight_layout()
    
    figname = '%s/%s_%s' %(plotDir, nc_prod, datum)
    fig.savefig(figname + '.png')
    fig.show()
    
    
    
    #: Initiate the data
    date = datetime(year,day,mon,hour,min)
    dat_new = SAFNWC_CTTH(date,'msg1',BBname='SAFBox',fullname=fil_new)
    dat_old = SAFNWC_CTTH(date,'msg1',BBname='SAFBox',fullname=fil_old)
    #: Change to new
    
    gg_new = geosat.GeoGrid('FullAMA_SAFBox')
    gg_old = geosat.GeoGrid('FullAMA_SAFBox')
    dat_new._CTTH_PRESS()
    dat_old._CTTH_PRESS()
#     dat.show('CTTH_PRESS')
    p1_new = geosat.SatGrid(dat_new,gg_new)
    p1_old = geosat.SatGrid(dat_old,gg_old)
#     p1.var.update({'CTTH_PRESS':mask_prod_new})
    p1_new._sat_togrid('CTTH_PRESS')
    p1_old._sat_togrid('CTTH_PRESS')
    import copy
    chart(p1_new, 'CTTH_PRESS')
    chart(p1_old, 'CTTH_PRESS')
    p1_diff = copy.copy(p1_new)
    pdb.set_trace()
    p1_diff.var['CTTH_PRESS'].data[:] = np.where((p1_new.var['CTTH_PRESS'].data==mask_prod_new.fill_value) | (p1_old.var['CTTH_PRESS'].data==mask_prod_old.fill_value), mask_prod_new.fill_value, p1_new.var['CTTH_PRESS'].data-p1_old.var['CTTH_PRESS'].data)
#     np.where(np.isnan(prod_diff), mask_prod_new.fill_value, prod_diff)
    chart(p1_diff, 'CTTH_PRESS')
    print('hmm')
    pdb.set_trace()
    p1.chart('CTTH_PRESS')
    pdb.set_trace()
    p1_diff.var['CTTH_PRESS'].data
    
    
    
    
    
    
    
    
    
#     (Pdb) mask_prod_new
# masked_array(
#   data=[[--, --, --, ..., --, --, --],
#         [--, --, --, ..., --, --, --],
#         [--, --, --, ..., --, --, --],
#         ...,
#         [--, --, --, ..., --, --, --],
#         [--, --, --, ..., --, --, --],
#         [--, --, --, ..., --, --, --]],
#   mask=[[ True,  True,  True, ...,  True,  True,  True],
#         [ True,  True,  True, ...,  True,  True,  True],
#         [ True,  True,  True, ...,  True,  True,  True],
#         ...,
#         [ True,  True,  True, ...,  True,  True,  True],
#         [ True,  True,  True, ...,  True,  True,  True],
#         [ True,  True,  True, ...,  True,  True,  True]],
#   fill_value=65535,
#   dtype=float32)
# (Pdb) dat.var['CTTH_PRESS']
# masked_array(
#   data=[[--, --, --, ..., --, --, --],
#         [--, --, --, ..., --, --, --],
#         [--, --, --, ..., --, --, --],
#         ...,
#         [65535.0, 65535.0, 65535.0, ..., 65535.0, 65535.0, 65535.0],
#         [65535.0, 65535.0, 65535.0, ..., 65535.0, 65535.0, 65535.0],
#         [65535.0, 65535.0, 65535.0, ..., 65535.0, 65535.0, 65535.0]],
#   mask=[[ True,  True,  True, ...,  True,  True,  True],
#         [ True,  True,  True, ...,  True,  True,  True],
#         [ True,  True,  True, ...,  True,  True,  True],
#         ...,
#         [False, False, False, ..., False, False, False],
#         [False, False, False, ..., False, False, False],
#         [False, False, False, ..., False, False, False]],
#   fill_value=65535.0,
#   dtype=float32)
# (Pdb) dat.copy()
# *** AttributeError: 'SAFNWC_CTTH' object has no attribute 'copy'
# (Pdb) import copy
# (Pdb) copy(dat)
# *** TypeError: 'module' object is not callable
# (Pdb) copy.copy(dat)
# <SAFNWCnc.SAFNWC_CTTH object at 0x7f98af2356d0>
# (Pdb) copy.deepcopy(dat)
# *** NotImplementedError: Dataset is not picklable
# (Pdb) copy.deep_copy(dat)
# *** AttributeError: module 'copy' has no attribute 'deep_copy'

    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
    