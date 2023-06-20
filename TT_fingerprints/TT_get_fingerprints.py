import numpy as np
import xarray as xr

def get_area_values(area_coords, data,lat_str = 'latitude',lon_str = 'longitude', inv=False,lon360=False):
    if inv:
        lats = data[lat_str].values[::-1]
    else:
        lats = data[lat_str].values
    if lon360:
        lons = data[lon_str].values-180
    else:
        lons = data[lon_str].values
    la = len(lats)
    lo = len(lons)
    n = la * lo # number of grid points
    indices = np.arange(n).reshape((la, lo)) # array of shape of grid points
    # get array of all coordinate pairs by repeating and tiling:
    coords = np.transpose([np.repeat(lats, lo), np.tile(lons, la)])
    lat_max = lats[np.argmin(np.abs(lats - np.max(area_coords[:, 0])))]
    lat_min = lats[np.argmin(np.abs(lats - np.min(area_coords[:, 0])))]
    lon_max = lons[np.argmin(np.abs(lons - np.max(area_coords[:, 1])))]
    lon_min = lons[np.argmin(np.abs(lons - np.min(area_coords[:, 1])))]
    coord_indices = indices[np.where(lats == lat_min)[0][0] : np.where(lats == lat_max)[0][0] + 1 , np.where(lons == lon_min)[0][0] : np.where(lons == lon_max)[0][0] + 1 ]
    coord_indices = coord_indices.flatten()
    ## lat-lon coorinates in same area as area_coords:
    sq_ds_coords = coords[coord_indices]
    d = np.zeros((len(area_coords), len(coord_indices)))
    
    for i in range(len(area_coords)):
        for j in range(len(coord_indices)):
            d[i, j] = np.sum(np.abs(area_coords[i] - sq_ds_coords[j]))
    indices_area = np.unique(coord_indices[np.argmin(d, axis = 1)])
    
    ds_coords = coords[indices_area]
    area_values = []
    for i in range(len(ds_coords)):
        if lat_str == 'latitude':
            area_values.append(data.sel(longitude=ds_coords[i,1]).sel(latitude=ds_coords[i,0]).values)
        elif lat_str == 'lat':
            area_values.append(data.sel(lon=(ds_coords[i,1]+180)).sel(lat=ds_coords[i,0]).values)
        else:
            print('Different lon/lat keys')
    return ds_coords, np.array(area_values)

##########

# First run the bash cdo_script.sh to get the global means, gridarea files and the seasonal dipole fingerprint, then run this script
# NB: as the dataset files are regularly updated with the latest months, some of the numbers in this file are dependant on when the data was downloaded. When opened the data should be indexed to only take years with full data - in this script some files end at 2023, some at 2022
# Start years are: HadISST1 - 1870; ERSTTv5 - 1854; HadCRUT and HadSST - 1850 

#############
# SPG fingerprints
###########

######## HadISST1

sst = xr.open_dataset('HadISST_sst_masked.nc').sst
sst_masked = sst.where(sst!=-1000)[:-11]
weights = xr.open_dataset('HadISST_sst_masked_gridarea.nc').cell_area
Hgl = xr.open_dataset('HadISST_sst_masked_gl.nc').sst.squeeze()[:-11]
area_caesar_coords = np.loadtxt('area_ceaser.txt')
ds_coords, area_values = get_area_values(area_caesar_coords,sst_masked,inv=True)
ds_coords, area_weights = get_area_values(area_caesar_coords,weights,inv=True)

Hspg = np.nanmean((area_values*(np.repeat(area_weights,1824).reshape(744,1824))),axis=0)/np.nanmean(area_weights)
Hspg_seas = np.tile(np.nanmean(Hspg.reshape(152,12),axis=0),152)
Hgl_seas = np.tile(np.nanmean(Hgl.values.reshape(152,12),axis=0),152)

Hexact_amoc = Hspg-Hspg_seas -  Hgl+Hgl_seas
Hexact_amoc2 = Hspg-Hspg_seas -  2*(Hgl-Hgl_seas)

ds = Hgl.copy(data=Hexact_amoc)
ds.to_netcdf('final/HadISST_C18.nc')
ds = Hgl.copy(data=Hexact_amoc2)
ds.to_netcdf('final/HadISST_C18_2GMT.nc')

########## ERSSTv5

Esst = xr.open_dataset('sst.mnmean.nc').sst[:-4]
Eweights = xr.open_dataset('sst.mnmean_gridarea.nc').cell_area
area_caesar_coords = np.loadtxt('area_ceaser.txt')
ds_coords, area_values = get_area_values(area_caesar_coords,Esst,inv=True,lat_str='lat',lon_str='lon',lon360=True)
ds_coords, area_weights = get_area_values(area_caesar_coords,Eweights,inv=True,lat_str='lat',lon_str='lon',lon360=True)
Espg = np.nanmean((area_values*(np.repeat(area_weights,2028).reshape(208,2028))),axis=0)/np.nanmean(area_weights)
Espg_seas = np.tile(np.nanmean(Espg.reshape(169,12),axis=0),169)
Egl = xr.open_dataset('sst.mnmean_gl.nc').sst[:-4].squeeze().values
Egl_seas = np.tile(np.nanmean(Egl.reshape(169,12),axis=0),169)
Eexact_amoc = Espg-Espg_seas -  Egl+Egl_seas
Eexact_amoc2 = Espg-Espg_seas -  2*(Egl-Egl_seas)
ds = xr.open_dataset('sst.mnmean_gl.nc').sst[:-4].squeeze().copy(data=Eexact_amoc)
ds.to_netcdf('final/ERSSTv5_C18.nc')
ds = xr.open_dataset('sst.mnmean_gl.nc').sst[:-4].squeeze().copy(data=Eexact_amoc2)
ds.to_netcdf('final/ERSSTv5_C18_2GMT.nc')

######## HadCRUT5

Csst = xr.open_dataset('HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean_SST.nc').tas_mean[:-2]
Cweights = xr.open_dataset('HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean_SST_gridarea.nc').cell_area
area_caesar_coords = np.loadtxt('area_ceaser.txt')
ds_coords, area_values = get_area_values(area_caesar_coords,Csst,inv=False,lon360=False)
ds_coords, area_weights = get_area_values(area_caesar_coords,Cweights,inv=False,lon360=False)
Cspg = np.nanmean((area_values*(np.repeat(area_weights,2064).reshape(43,2064))),axis=0)/np.nanmean(area_weights)

Cspg_seas = np.tile(np.nanmean(Cspg.reshape(172,12),axis=0),172)
Cgl = xr.open_dataset('HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean_SST_gl.nc').tas_mean[:-2].squeeze().values
Cgl_seas = np.tile(np.nanmean(Cgl.reshape(172,12),axis=0),172)
Cexact_amoc = Cspg-Cspg_seas -  Cgl+Cgl_seas
Cexact_amoc2 = Cspg-Cspg_seas -  2*(Cgl-Cgl_seas)

Cgl = xr.open_dataset('HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean_SST_gl.nc').tas_mean[:-2].squeeze()
ds = Cgl.copy(data=Cexact_amoc)
ds.to_netcdf('final/HadCRUT5_C18.nc')
ds = Cgl.copy(data=Cexact_amoc2)
ds.to_netcdf('final/HadCRUT5_C18_2GMT.nc')

########## HadSST4

Dsst = xr.open_dataset('HadSST.4.0.1.0_median.nc').tos[:-2]
Dweights = xr.open_dataset('HadSST.4.0.1.0_median_gridarea.nc').cell_area
area_caesar_coords = np.loadtxt('area_ceaser.txt')
ds_coords, area_values = get_area_values(area_caesar_coords,Dsst,inv=False,lon360=False)
ds_coords, area_weights = get_area_values(area_caesar_coords,Dweights,inv=False,lon360=False)
Dspg = np.nanmean((area_values*(np.repeat(area_weights,2064).reshape(43,2064))),axis=0)/np.nanmean(area_weights)

Dspg_seas = np.tile(np.nanmean(Dspg.reshape(172,12),axis=0),172)
Dgl = xr.open_dataset('HadSST.4.0.1.0_median_gl.nc').tos[:-2].squeeze().values
Dgl_seas = np.tile(np.nanmean(Dgl.reshape(172,12),axis=0),172)
Dexact_amoc = Dspg-Dspg_seas -  Dgl+Dgl_seas
Dexact_amoc2 = Dspg-Dspg_seas -  2*(Dgl-Dgl_seas)

Dgl = xr.open_dataset('HadSST.4.0.1.0_median_gl.nc').tos[:-2].squeeze()
ds = Dgl.copy(data=Dexact_amoc)
ds.to_netcdf('final/HadSST4_C18.nc')
ds = Dgl.copy(data=Dexact_amoc2)
ds.to_netcdf('final/HadSST4_C18_2GMT.nc')

#############
# rectangle area fingerprints
###########

datasets= ['HadISST','ERSSTv5','HadCRUT5','HadSST4']
keys = ['sst','sst','tas_mean','tos']
paths = ['HadISST_sst_masked', 'sst.mnmean', 'HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean_SST', 'HadSST.4.0.1.0_median']
cuts = [-11,-4,-2,-2]

for i, dataset in enumerate(datasets):
    ts = xr.open_dataset(paths[i]+'_M20.nc')[keys[i]].squeeze()[:cuts[i]]
    ts_dseas = ts.groupby('time.month')-ts.groupby('time.month').mean(dim="time")
    ts_dseas.to_netcdf('final/{}_M20.nc'.format(dataset))

    ts = xr.open_dataset(paths[i]+'_dipole.nc')[keys[i]].squeeze()[:cuts[i]]
    ts_dseas = ts.groupby('time.month')-ts.groupby('time.month').mean(dim="time")
    ts_dseas.to_netcdf('final/{}_dipole.nc'.format(dataset))