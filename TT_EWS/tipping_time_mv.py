import numpy as np
import statsmodels.api as sm
import scipy.stats as st
from scipy.optimize import curve_fit
import xarray as xr
from scipy import optimize

def runMLE(x, w):
    n = x.shape[0]
    mu = np.full(len(x),np.nan)
    ar1 = np.full(len(x),np.nan)
    var = np.full(len(x),np.nan)
    for i in range(w // 2, n - w // 2):
        xw = x[i - w // 2 : i + w // 2 + 1]
        if len(np.where(np.isnan(xw))[0]) < w//2:
            xw2 = xw[np.where(~np.isnan(xw))]
            nw = len(xw2)
            mu[i] = np.sum(xw2)/(nw+1)
            ar1[i] = sm.tsa.acf(xw2)[1]
            var[i]= np.sum(np.square(xw2[1:]-xw2[:-1]*ar1[i]-mu[i]*(np.full(nw-1,1)-ar1[i])))/(nw*(1-ar1[i]**2))

        else:
            mu[i] = np.nan
            ar1[i] = np.nan
            var[i]= np.nan

    return mu,ar1,var

def ffit(tt, t0, y0, tau_r):
    return y0*(1-(np.heaviside(tt-t0,1)*(tt-t0))/tau_r)

def get_EWS_fit(Tw, ts, initial_guess, bounds):
# for a given Tw get the EWS timeseries and fit ffit to Alambda
    mu, ar1, var = runMLE(ts,w=Tw)
    alpha = -np.log(ar1)
    sigma2 = 2*alpha*var
    A_lamb=(sigma2/(4*var))**2
    y = A_lamb[Tw//2:-Tw//2]
    xx = np.arange(Tw//2,len(y)+Tw//2)
    idc = np.where(~np.isnan(y))
    p, e = optimize.curve_fit(ffit, xx[idc], y[idc],p0=initial_guess,bounds=bounds)
    return mu, ar1, var, alpha, sigma2, A_lamb, p[0], p[1], p[2]

def get_Tw_lsqs(Twss, ts, initial_guess, bounds):
# for each window Tw get the best fit t0, y0 and tau_r, and calculate least square residual
    Tw_lsqs = np.full(len(Twss),np.nan)
    for i, Tw in enumerate(Twss):
        mu, ar1, var, alpha, sigma2, A_lamb, t0, y0, tau_r = get_EWS_fit(Tw, ts, initial_guess, bounds)
        y = A_lamb[Tw//2:-Tw//2]
        xx = np.arange(Tw//2,len(y)+Tw//2)
        idc = np.where(~np.isnan(y))
        Tw_lsqs[i] = np.sum(np.square(y[idc]-ffit(xx[idc],t0, y0, tau_r)))
    return Tw_lsqs

datasets= ['HadISST','ERSSTv5','HadCRUT5','HadSST4']
keys = ['sst','sst','tas_mean','tos']
lags = [20*12,4*12,0,0]
indices = ['C18_2GMT','C18','dipole','M20']
##############
dataset='HadSST4'
index='dipole'

key_dict=dict(zip(datasets, keys))
lag_dict = dict(zip(datasets, lags))


ds = xr.open_dataset('../TT_fingerprints/final/{}_{}.nc'.format(dataset,index))
amoc = ds[key_dict[dataset]].values
amoc = amoc - np.nanmean(amoc)
ts = amoc

i1910 = 60*12 - lag_dict[dataset]  # each timeseries starts at a different year, get lag to 1910

initial_guess = [i1910+20*12,0.02, 100*12]
bounds = ([i1910, 0, 0],[i1910+40*12, 1, np.inf])

Twss = np.arange(45*12,65*12,12) # range of windows from 45 to 65 years

Tw_lsqs = get_Tw_lsqs(Twss, ts, initial_guess, bounds)

idx = np.argmin(Tw_lsqs)
Tw_fin = Twss[idx]

mu, ar1, var, alpha, sigma2, A_lamb, t0, y0, tau_r = get_EWS_fit(Tw_fin,ts, initial_guess, bounds)

array = xr.Dataset(
            data_vars = dict(
                mu=(['time'],mu),
                var=(['time'],var),
                ar1=(['time'],ar1),
                alpha=(['time'],alpha),
                sigma2=(['time'],sigma2),
                A_lamb=(['time'],A_lamb),
                params=(['pdim'],[Tw_fin,t0, y0, tau_r])
                ),
            coords = dict(
                time = ds.time,
                pdim=xr.DataArray(['Tw','t0','y0','tau_r'], dims="pdim", coords=dict(pdim=("pdim", ['Tw','t0','y0','tau_r'])))
                ),
        )

array.to_netcdf('{}_{}_EWSfit.nc'.format(dataset,index))