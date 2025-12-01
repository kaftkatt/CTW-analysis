#!/usr/bin/env python
# coding: utf-8

# # Animation

# In[1]:

import SVBfunc 

import xarray as xr
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
from matplotlib.animation import FuncAnimation
from matplotlib import animation
import ffmpeg
import pylab as pl
from math import radians, cos
from scipy.integrate import trapezoid


coast='original'
tstart=2

dirn = '/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
dirw = '/home/athelandersson/NETCDFs/' + str(coast) + '/'

dsw, dsn = SVBfunc.loadNetCDFs(dirw, dirn, 'dynVars',tstart)

dswPhi, dsnPhi = SVBfunc.loadNetCDFs(dirw, dirn, 'phiHyd',tstart)

pathVav='/home/athelandersson/CTW-analysis/Files/' + str(coast)+ '/MeanVVEL.nc'
dsVav=xr.open_dataset(pathVav)

pathPav='/home/athelandersson/CTW-analysis/Files/' + str(coast)+ '/MeanPHIHYD.nc'
dsPav=xr.open_dataset(pathPav)

Z=dswPhi[0].Z.values

hFacCw = dsnPhi[0].hFacC
hFacCusew=hFacCw.values

hfa = np.ma.masked_values(hFacCusew, 0)
maskw = np.ma.getmask(hfa)

Zmat=np.ma.masked_array(np.tile(Z,(len(maskw[0,0,:]),len(maskw[0,:,0]),1)).T,mask=maskw)
dz=Zmat[:-1]-Zmat[1:]

timeout=np.arange(2880,14400,20)

fluxout=np.zeros((len(dsw),len(maskw[0,:,0]),len(maskw[0,0,:])))
for i in np.arange(0,len(dsw),1):
    fluxoutPart=np.zeros((len(dsw[0].time.values),len(maskw[0,:,0]),len(maskw[0,0,:])))
    for t in np.arange(0,len(dsw[0].time.values),1):
        for xin in range(len(maskw[0,0,:])):
             velanoma=np.zeros((len(Z),len(maskw[0,:,0])))
             phianoma=np.zeros((len(Z),len(maskw[0,:,0])))
             for depin in range(len(Z)):
                  velb=dsw[i].VVEL[t,depin,:,xin].values
                  veln=dsn[i].VVEL[t,depin,:,xin].values
                  phib=dswPhi[i].PHIHYD[t,depin,:,xin].values 
                  phin=dsnPhi[i].PHIHYD[t,depin,:,xin].values  
                  phianoma[depin,:]=(phib-phin)-dsPav.avePHIHYD[depin,:,xin].values
                  velanoma[depin,:]=(velb-veln)-dsVav.aveVVEL[depin,:,xin].values
             fluxoutPart[t,:,xin]=SVBfunc.fluxInternalWave((np.ma.masked_array(np.array(phianoma),mask=maskw[:,:,xin])),np.ma.masked_array(np.array(velanoma),mask=maskw[:,:,xin]),dz[:,:,xin])
    print(t)
    fluxout[i,:,:]=np.nanmean(fluxoutPart,axis=0)

pathE='/home/athelandersson/CTW-analysis/Files/' + str(coast)+ '/EfluxVNov.nc'
dsE = xr.Dataset({"Energyflux": (("time","lat","lon"), np.array(fluxout)),
					    },
				   coords ={
					 "lon" : LON.values,
					 "lat": LAT.values,
					 "time": np.array(timeout)
				      },
				       )

dsE.to_netcdf(pathE)



