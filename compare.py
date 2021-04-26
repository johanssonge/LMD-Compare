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
from datetime import datetime, timedelta

sys.path.append('/home/erikj/Projects/STC/pylib')
import geosat  # @UnresolvedImport
from SAFNWCnc import SAFNWC_CTTH  # @UnresolvedImport

satConversion = {'hima08': 'himawari', 'msg1': 'msg1', 'msg3': 'msg3'}
#----------------------------------------------------

def chart(obts,field,cmap='jet',clim=[180.,300.],txt=[''],subgrid=None, block=True, xlocs=None, figsize= None, show=True, typ1='', typ2='', datum='', sat=''):
        if not isinstance(obts, list):
            obts = [obts]
        obt_t = obts[0]
        # test existence of key field
        if field not in obt_t.var.keys():
            print ('undefined field')
            return
        if subgrid == None:
            geogrid = obt_t.geogrid
        else:
            geogrid = subgrid
        if 'FullAMA' in geogrid.gridtype:
            fig = plt.figure(figsize=[10, 6])
        elif figsize is not None:
            fig = plt.figure(figsize=figsize)
        else:
            fig = plt.figure(figsize=[11,4])
#         fig.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
        fig.suptitle('%s, %s, %s' %(sat, field, datum))
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
        for f, obt in enumerate(obts):
            ax = fig.add_subplot(len(obts), 1, f+1, projection = proj)
            if subgrid == None:
                if isinstance(obt, dict):
                    plotted_field = obt[field]
                    plotted_field = np.where(plotted_field == 65535.0, np.nan, plotted_field)
                    clim = [-20, 20]
                    cmap='RdBu_r'
                else:
                    plotted_field = obt.var[field]
                    cmap=None
            else:
                # extraction in subgrid
                plotted_field = obt.var[field][geogrid.corner[1]:geogrid.corner[1]+geogrid.box_biny,
                                                geogrid.corner[0]:geogrid.corner[0]+geogrid.box_binx]
            im = ax.imshow(plotted_field, transform=proj, interpolation='nearest',
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
    
            if (cm_lon == 180) & (xlocs == None): 
                xlocs = [0,30,60,90,120,150,180,-150,-120,-90,-60,-30]
            gl = ax.gridlines(draw_labels=True, xlocs=xlocs,
                          linewidth=2, color='gray', alpha=0.5, linestyle='--')
            if f in [0,1]:
                gl.bottom_labels = False
            gl.top_labels = False
            gl.right_labels = False
            #gl.xformatter = LONGITUDE_FORMATTER
            #gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': fs}
            gl.ylabel_style = {'size': fs}
            #gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
            ax.set_title(txt[f],fontsize=fs)
#             plt.title(txt,fontsize=fs)
            # plot adjusted colorbar and show
            
#             barticks = [valminmax*-1, valminmax*-0.75, valminmax*-0.5, valminmax*-0.25, 0, valminmax*0.25, valminmax*0.5, valminmax*0.75, valminmax]
#             cbar_ax = fig.add_axes([0.2, 0.25, 0.6, 0.01])
#             cbar = fig.colorbar(im, orientation='horizontal', cax=cbar_ax, ticks=barticks)  # @UnusedVariable
            
            axpos = ax.get_position()
#             pos_x = axpos.x0 + axpos.x0 + axpos.width + 0.01
            pos_x = axpos.x0 + axpos.width + 0.01
#             pos_cax = fig.add_axes([pos_x,axpos.y0,0.04,axpos.height])
            if f in [0,1]:
                barticks = [clim[0], clim[0] + 40, clim[0] + 80, clim[1]]
            else:
                barticks = [clim[0], clim[0] + 10, 0, clim[0] + 30, clim[1]]
            cbar_ax = fig.add_axes([pos_x, axpos.y0, 0.01, axpos.height])
            cbar = fig.colorbar(im, orientation='vertical', cax=cbar_ax, ticks=barticks)  # @UnusedVariable
#             cbar=fig.colorbar(im)#, cax=pos_cax)
            cbar.ax.tick_params(labelsize=fs)
        figname = '/scratch/erikj/Proj1/Plots/map_diff_%s_%s_%s_%s_%s' %(sat, field, txt[0], txt[2], datum)
        fig.savefig(figname + '.png')
        if show:
            fig.show()
            pdb.set_trace()
        return ax


def bars(obts,field,clim=[180.,300.],txt=[''], show=True, datum='', sat='', v2018s=None, days=1, hours=1):
        if isinstance(obts, dict) and (v2018s is not None):
            obts_use = []
            txt = []
            for v in v2018s:
                obts_use.append(np.asarray(obts['v2018-%d_diff_agg' %v]))
                txt.append('v2018-%d' %v)
            obts_use.append(np.asarray(obts['v2016_diff_agg']))
            txt.append('v2016')
            figname_txt = 'v2018s_%s' %datum
        elif not isinstance(obts, list):
            obts_use = [obts]
            figname_txt = '%s_%s_%s' %(txt[0], txt[1], datum)
            
        else:
            obts_use = obts
            figname_txt = '%s_%s_%s' %(txt[0], txt[1], datum)
            
        # test existence of key field
#         pdb.set_trace()
#         if field not in obt_t.var.keys():
#             print ('undefined field')
#             return
        
        suptitle_txt = '%s, %s, %s' %(sat, field, datum)
        if days > 1:
            suptitle_txt = suptitle_txt + ', days=%d' %days
            figname_txt= figname_txt + '_d%d' %days
        if hours > 1:
            suptitle_txt = suptitle_txt + ', hours=%d' %hours
            figname_txt = figname_txt + '_h%d' %hours
        fz = (11, len(v2018s)*4)#, 4)
        fig = plt.figure(figsize=fz)
        fig.suptitle(suptitle_txt)
        fs = 15
        for f, obt in enumerate(obts_use):
            ax = fig.add_subplot(len(obts_use), 1, f+1)
            if isinstance(obt, list) or isinstance(obt, np.ndarray):
                plotted_field = obt
            elif isinstance(obt, dict):
                plotted_field = obt[field]
            else:
                plotted_field = obt.var[field].data[:]
            if (txt[f] == 'Diff') or ((v2018s is not None) and (f > 0)):
                rang = (-30, 30)
                bins = 100
            else:
                rang = (clim[0], clim[1])
                bins = 150
            mapv = ~(plotted_field == 65535.0)
            #: hPa
            try:
                plotted_field = plotted_field[mapv]# / 100.
            except:
                pdb.set_trace()
#             plotted_field = np.where(plotted_field == 65535.0, np.nan, plotted_field)
#             ax.hist(plotted_field, range=rang)#, aspect=1.)
            ax.hist(plotted_field, bins=bins, range=rang)#, aspect=1.)
            ax.set_title(txt[f],fontsize=fs)
#             if f in [0,1]:
#             ax.set_ylim((0, 150))
#             ax.set_yticks([0, 50, 100, 150])
#             gl.xlabel_style = {'size': fs}
#             gl.ylabel_style = {'size': fs}
#             #gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
# #             plt.title(txt,fontsize=fs)
#             # plot adjusted colorbar and show
#             
# #             barticks = [valminmax*-1, valminmax*-0.75, valminmax*-0.5, valminmax*-0.25, 0, valminmax*0.25, valminmax*0.5, valminmax*0.75, valminmax]
# #             cbar_ax = fig.add_axes([0.2, 0.25, 0.6, 0.01])
# #             cbar = fig.colorbar(im, orientation='horizontal', cax=cbar_ax, ticks=barticks)  # @UnusedVariable
#             
#             axpos = ax.get_position()
# #             pos_x = axpos.x0 + axpos.x0 + axpos.width + 0.01
#             pos_x = axpos.x0 + axpos.width + 0.01
# #             pos_cax = fig.add_axes([pos_x,axpos.y0,0.04,axpos.height])
#             if f in [0,1]:
#                 barticks = [clim[0], clim[0] + 40, clim[0] + 80, clim[1]]
#             else:
#                 barticks = [clim[0], clim[0] + 10, 0, clim[0] + 30, clim[1]]
#             cbar_ax = fig.add_axes([pos_x, axpos.y0, 0.01, axpos.height])
#             cbar = fig.colorbar(im, orientation='vertical', cax=cbar_ax, ticks=barticks)  # @UnusedVariable
# #             cbar=fig.colorbar(im)#, cax=pos_cax)
#             cbar.ax.tick_params(labelsize=fs)
        plt.tight_layout()
        figname = '/scratch/erikj/Proj1/Plots/hist_diff_%s_%s_%s' %(sat, field, figname_txt)
        fig.savefig(figname + '.png')
        if show:
            fig.show()
            pdb.set_trace()
        return ax



def saveLatLon(ncf, fn):
    basename = os.path.basename(fn)
    LatLonDir = '/scratch/erikj/Data/Nwcgeo/LatLon'
    satname = basename.split('_')[3].lower()
    sfn = '%s/%s_LonLat.h5' %(LatLonDir, satname)
    if not os.path.isfile(sfn):
        pdb.set_trace()
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
    mainDir_new = '/scratch/b.legras/NWCGEO-APPLI-2018.1-x/export'
#     mainDir_old = '/scratch/b.legras/sats/himawari/safnwc/netcdf/2017/'
    plotDir = '/scratch/erikj/Data/Nwcgeo/Compare/Plots'
#     fil = 'S_NWC_CMA_MSG1_FULLAMA-VISIR_20170515T053000Z.nc'
    sat='MSG1'# 'HIMA08' #'MSG1'
    product = 'CTTH'
    clim=[0.,500.]
    if sat == 'HIMA08':
        sat_dir = 'himawari'
        sat_fil = 'HIMAWARI08'
    else:
        sat_dir = sat.lower()
        sat_fil = sat
    mainDir_old = '/scratch/b.legras/sats/%s/safnwc-SAFBox/netcdf' %sat_dir
    if mainDir_new.split('/')[3] == 'sats':
        typ1 = '2016'
    else:
        typ1 = mainDir_new.split('/')[3].replace('NWCGEO-APPLI-', '')
    if mainDir_old.split('/')[3] == 'sats':
        typ2 = '2016'
    else:
        typ2 = mainDir_old.split('/')[3].replace('NWCGEO-APPLI-', '')
    
    nc_prod = product.lower()
    if nc_prod == 'ctth':
        nc_prod = 'ctth_pres'
        grid_type = 'CTTH_PRESS'
    p1_new_agg = []
    p1_old_agg = []
    p1_diff_agg = []
#     all_dates = [datetime(2017,7,22,20,0), datetime(2017,7,22,20,20), datetime(2017,7,22,21,0), \
#                 datetime(2017,7,24,8,40), datetime(2017,7,24,9,0), datetime(2017,7,24,9,20), \
#                 datetime(2017,7,26,15,40), datetime(2017,7,26,16,0), datetime(2017,7,26,16,20)]
    #: Start date
    st_date = datetime(2017,8,3,8,0)
    st_datum = st_date.isoformat().replace('-', '').replace(':', '')[0:-2]
    
    all_dates = []
    #: Add days
    nrdays = 10
    nrhours = 8
    for d in range(nrdays):
        #: Add hours
        for h in range(nrhours):
            all_dates.append(st_date + timedelta(days=d) + timedelta(hours=h))
    #: List for wrong files
    wrong_files = []
    #: 2018 versions
    if sat == 'MSG1':
        v_2018 = list(range(1,13+1,2))
    elif sat == 'HIMA08':
        v_2018 = list(range(2,10+1,2))
    
#     v_2018 = [1,2,3,4,5,6,7,8,9,10,11,12,13]
    v_0 = v_2018[0]
    #: Dict for results.
    results = {}
    for v in v_2018:
        results.update({'v2018-%d_p1' %v: ''})
        results.update({'v2018-%d_diff' %v: ''})
        results.update({'v2018-%d_diff_agg' %v: []})
    results.update({'v2016_p1': ''})
    results.update({'v2016_diff': ''})
    results.update({'v2016_diff_agg': []})

    for date in all_dates:
        year = date.year
        mon = date.month
        day = date.day
#         hour = 20
#         minut = 0
        datum = date.isoformat().replace('-', '').replace(':', '')[0:-2]
        print(datum)
#         datum = '%d%02d%02dT%02d%02d' %(year, mon, day, hour, minut)
        #: V2018
        for v in v_2018:
            print('2018.1-%d' %v)
            v2018_filestr = ('%s/%s/S_NWC_%s_%s_FULLAMA-*_%s*.nc' %(mainDir_new, product, product, sat, datum)).replace('2018.1-x', '2018.1-%d' %v)
            v2018_filedic = glob.glob(v2018_filestr)
            if len(v2018_filedic) == 0:
                wrong_files.append(v2018_filestr)
                continue
            v2018_filedic.sort()
            v2018_file = v2018_filedic[0]
            v2018_ncf = nc.Dataset(v2018_file)
            #: Read the required data
            v2018_maskprod = v2018_ncf[nc_prod][:]
            v2018_prod = v2018_maskprod.data.astype(float)
            v2018_prod[v2018_ncf[nc_prod][:].mask] = np.nan
            v2018_ncf.close()
            #: --- New Grid ---
            #: Initiate the data
            v2018_data = SAFNWC_CTTH(date,sat_dir,BBname='SAFBox',fullname=v2018_file)
            #: Change to new grid
            v2018_gg = geosat.GeoGrid('FullAMA_SAFBox')
            v2018_data._CTTH_PRESS()
            v2018_p1 = geosat.SatGrid(v2018_data, v2018_gg)
            v2018_p1._sat_togrid(grid_type)
            
            results['v2018-%d_p1' %v] = v2018_p1
            v2018_valid_mask = ((v2018_p1.var[grid_type].data != v2018_maskprod.fill_value) & \
                                (v2018_p1.var[grid_type].data >= clim[0]) & \
                                (v2018_p1.var[grid_type].data <= clim[1]))
            if v == v_0:
                results.update({'v2018-%d_fill_value' %v_0: v2018_maskprod.fill_value})
                results.update({'v2018-%d_valid_mask' %v_0: v2018_valid_mask})
                #: Not diff for v_0
                p1_diff_var = v2018_p1.var[grid_type].data[v2018_valid_mask]
                results['v2018-%d_diff_agg' %v].extend(p1_diff_var)
                continue
            valid_mask = ((results['v2018-%d_valid_mask' %v_0]) & (v2018_valid_mask))
#             p1_diff_var = np.where(valid_mask, 
#                                    results['v2018-%d_p1' %v_0].var[grid_type].data - v2018_p1.var[grid_type].data, 
#                                    v2018_maskprod.fill_value)
            p1_diff_var = results['v2018-%d_p1' %v_0].var[grid_type].data[valid_mask] - v2018_p1.var[grid_type].data[valid_mask]
            #: Aggregate
            results['v2018-%d_diff_agg' %v].extend(p1_diff_var)
            
            plot_single = False
            if plot_single:    
                p1_diff = {grid_type: p1_diff_var}
                bars([results['v2018-%d_p1' %v_0], v2018_p1, p1_diff], grid_type,clim=clim, txt=['v2018.1-%d' %v_0, 'v2018.1-%d' %v, 'Diff'], show=True, datum=datum, sat=sat)
                chart([results['v2018-%d_p1' %v_0], v2018_p1, p1_diff], grid_type,clim=clim, txt=['v2018.1-%d' %v_0, 'v2018.1-%d' %v, 'Diff'], show=True, datum=datum, sat=sat)
            
#     files_old = glob.glob('%s/%s/S_NWC_%s_%s_FULLAMA-*_%s*.nc' %(mainDir_old, product, product, sat, datum))
#     files_old = glob.glob('%s/%d_%02d_%02d/S_NWC_%s_%s_FULLAMA-*_%d%02d%02dT%02d%02d*.nc' %(mainDir_old, year, mon, day, product, sat, year, mon, day, hour, minut))
        v2016_filestr = '%s/%d/%d_%02d_%02d/S_NWC_%s_%s_FULLAMA-*_%s*.nc' %(mainDir_old, year, year, mon, day, product, sat_fil, datum)
        v2016_filedic = glob.glob(v2016_filestr)
        if len(v2016_filedic) == 0:
            wrong_files.append(v2016_filestr)
            continue
        v2016_filedic.sort()
        v2016_file = v2016_filedic[0]
        v2016_ncf = nc.Dataset(v2016_file)
        #: Read the required data
        v2016_maskprod = v2016_ncf[nc_prod][:]
        v2016_prod = v2016_maskprod.data.astype(float)
        v2016_prod[v2016_ncf[nc_prod][:].mask] = np.nan
        #: Take care of latitude and longitude
        sfn = saveLatLon(v2016_ncf, v2016_file)
        lat, lon = getLatLon(sfn)
        pdb.set_trace()
        v2016_ncf.close()
        #: --- New Grid ---
        #: Initiate the data
        v2016_data = SAFNWC_CTTH(date,sat_dir,BBname='SAFBox',fullname=v2016_file)
        #: Change to new grid
        v2016_gg = geosat.GeoGrid('FullAMA_SAFBox')
        v2016_data._CTTH_PRESS()
        v2016_p1 = geosat.SatGrid(v2016_data, v2016_gg)
        v2016_p1._sat_togrid(grid_type)
        
        results['v2016_p1'] = v2016_p1
        v2016_valid_mask = ((v2016_p1.var[grid_type].data != v2016_maskprod.fill_value) & \
                                (v2016_p1.var[grid_type].data >= clim[0]) & \
                                (v2016_p1.var[grid_type].data <= clim[1])) 
        valid_mask = ((results['v2018-%d_valid_mask' %v_0]) & (v2016_valid_mask))
#         p1_diff_var = np.where(valid_mask, 
#                                results['v2018-%d_p1' %v_0].var[grid_type].data - v2016_p1.var[grid_type].data, 
#                                v2016_maskprod.fill_value)
        p1_diff_var = results['v2018-%d_p1' %v_0].var[grid_type].data[valid_mask] - v2016_p1.var[grid_type].data[valid_mask]
        #: Aggregate
        results['v2016_diff_agg'].extend(p1_diff_var)
        plot_single = False
        if plot_single:
            p1_diff = {grid_type: p1_diff_var}
            bars([results['v2018-%d_p1' %v_0], v2016_p1, p1_diff], grid_type,clim=clim, txt=['v2018.1-%d' %v_0, 'v2016', 'Diff'], show=True, datum=datum, sat=sat)
            chart([results['v2018-%d_p1' %v_0], v2016_p1, p1_diff], grid_type,clim=clim, txt=['v2018.1-%d' %v_0, 'v2016', 'Diff'], show=True, datum=datum, sat=sat)

        
    #     pdb.set_trace()
        if False:
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
    
#     figname = '%s/%s_%s' %(plotDir, nc_prod, datum)
#     fig.savefig(figname + '.png')
#     fig.show()
    
    
    #     chart(p1_diff, 'CTTH_PRESS')
#     p1_new_agg_dict = {grid_type: p1_new_agg}
#     p1_old_agg_dict = {grid_type: p1_old_agg}
#     p1_diff_agg_dict = {grid_type: p1_diff_agg}
    bars(results, grid_type,clim=clim, txt=[typ1, typ2, 'Diff'], show=True, datum=st, sat=sat, v2018s=v_2018, days=nrdays, hours=nrhours)
    pdb.set_trace()
    
    sys.exit()
    
    
    
    
    
    
    
    
    
    gg = geosat.GeoGrid('FullAMA_SAFBox')
    fig = plt.figure(figsize=[10, 6])
    proj = ccrs.PlateCarree(central_longitude=0)
    ax = fig.add_subplot(111, projection=proj)
    ax.coastlines()
    ax.imshow(prod_new)
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
    