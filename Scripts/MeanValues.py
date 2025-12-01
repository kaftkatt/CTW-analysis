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

var='PHIHYD'
coast='smooth'
tstart=2

dirn = '/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
dirw = '/home/athelandersson/NETCDFs/' + str(coast) + '/'

dsw, dsn = SVBfunc.loadNetCDFs(dirw, dirn, 'phiHyd',tstart)

ind=0 #0 is day 2-3, 1 is day 3-4 until index 7 (day 9-10)
dep=0  #483.2 meter depth is the 55th element
t=0
tt=(((72*ind+t)*20)+2880)/60 # Gives amount of hours from start of the model, starts at hour 48 if ind=0 and t=0


Z=dsw[ind].Z.values

LON=dsw[ind].XC-360
LAT=dsw[ind].YC

aveVal=np.zeros((len(Z),len(LAT),len(LON)))
for dep in np.arange(0,len(Z),1):
	print('Depth nr' + str(dep))
	aveValpart=np.zeros((len(dsw),len(LAT),len(LON)))
	for tt in np.arange(0,len(dsw),1):
		Vals=np.zeros((72,len(LAT),len(LON)))
		exec(f'global Valw; Valw=dsw[tt].{var}')
		exec(f'global Valn; Valn=dsn[tt].{var}')
		for i in np.arange(0,72,1):
			Vals[i] = Valw[i,dep,:,:].values-Valn[i,dep,:,:].values
			print(i)
		aveValpart[tt]=np.nanmean(Vals,axis=0)
	aveVal[dep]=np.nanmean(aveValpart,axis=0)

pathE='/home/athelandersson/CTW-analysis/Files/' + str(coast)+ '/Mean' + str(var) + '.nc'
dsE = xr.Dataset({"ave"+str(var):(("Z","YC","XC"),aveVal)},
				   coords ={
					 "XC" : LON.values,
					 "YC": LAT.values,
					 "Z": Z,
					 "times": np.arange(0,7,1)
				      },
				       )

dsE.to_netcdf(pathE)





