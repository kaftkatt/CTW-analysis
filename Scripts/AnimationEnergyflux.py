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

dsw, dsn = SVBfunc.loadNetCDFs(dirw, dirn, 'velrot',tstart)

dswPhi, dsnPhi = SVBfunc.loadNetCDFs(dirw, dirn, 'phiHyd',tstart)


def fluxInternalWave(p,vel,dz,mask):
    Flux=np.ma.masked_array((p*vel)*dz,mask=mask)
    Flux[mask==True]=np.nan
    dz[mask==True]=np.nan
    return np.nansum(Flux,axis=0)/-np.nansum(dz,axis=0)

params = {'font.size': 22,
          'figure.figsize': (12, 8),
         'font.family':'sans'}
pl.rcParams.update(params)

#Index to call from the list of netcdfs
ind=0 #0 is day 2-3, 1 is day 3-4 until index 7 (day 9-10)
dep=0  #483.2 meter depth is the 55th element
t=0
tt=(((72*ind+t)*20)+2880)/60 # Gives amount of hours from start of the model, starts at hour 48 if ind=0 and t=0


depth=dswPhi[0].Depth[1:-1,1:-1].values
Z=dswPhi[ind].Zp1.values

LON=dswPhi[ind].XC[1:-1]-360
LAT=dswPhi[ind].YC[1:-1]

hFacCw = dsnPhi[ind].hFacC[:,1:-1,1:-1]
hFacCusew=hFacCw.values

hfa = np.ma.masked_values(hFacCusew, 0)
maskw = np.ma.getmask(hfa)

dz=np.ma.masked_array(np.tile(Z[:-1]-Z[1:],(len(maskw[0,0,:]),len(maskw[0,:,0]),1)).T,mask=maskw)

timeout=np.arange(2880,14400,20)
print(len(timeout))

aveP=[]
aveVel=[]
for i in np.arange(0,576,1):
    t=i
    vel= SVBfunc.get_snapshot_at_level( t,dep,dsw,dsn,'CSHORE')
    phi=SVBfunc.get_snapshot_at_level( t,dep,dswPhi,dsnPhi,'PHIHYDiw')
    print(t)
    aveP.append(np.nanmean(np.ma.masked_array(phi,mask=maskw)))
    aveVel.append(np.nanmean(np.ma.masked_array(vel,mask=maskw)))

fluxout=np.zeros((576,len(LAT),len(LON)))
for i in np.arange(0,576,1):
    t=i
    vel= SVBfunc.get_snapshot_at_level( t,dep,dsw,dsn,'CSHORE') 
    phi=SVBfunc.get_snapshot_at_level( t,dep,dswPhi,dsnPhi,'PHIHYDiw')
    print(t)
    fluxPlot=fluxInternalWave(phi-np.mean(aveP),vel-np.mean(aveVel),dz,maskw)
    fluxout[i,:,:]=fluxPlot

fluxIN=np.mean(fluxout,axis=0)

pathE='/home/athelandersson/CTW-analysis/Files/' + str(coast)+ '/EfluxCshore.nc'
dsE = xr.Dataset({"Energyflux": (("time","lat","lon"), np.array(fluxout)),
					 "avePhi":(("time"),np.array( aveP)),
					 "aveVelocity":(("time"), np.array(aveVel))
					    },
				   coords ={
					 "lon" : LON.values,
					 "lat": LAT.values,
					 "time": np.array(timeout)
				      },
				       )

dsE.to_netcdf(pathE)



# In[29]:


fig, ax = plt.subplots()

vmin=-1e-8 #np.max(np.ma.masked_array(fluxIN,mask=maskw[0,:,:]))
vmax=1e-8 #np.max(np.ma.masked_array(fluxIN,mask=maskw[0,:,:]))

xlab='Longitude [°]'
ylab='Latitude [°]'

ax.set_facecolor('wheat')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(fluxIN,mask=maskw[0,:,:]),
                    cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth,  colors=['0.2','0.4'], 
                levels=[0,500])
ax.set(xlabel=xlab, ylabel=ylab)

ind_lon = [-115.11813068276555, -115.939167, -116.605833, -117.1625, -118.24368, -119.714167, -120.471439,
           -120.7586085906775]
ind_lat = [27.850440699318973, 30.556389, 31.857778, 32.715, 34.05223, 34.425833, 34.448113, 35.17364705813524]

for kk, ll, lab in zip(ind_lon, ind_lat,
                       ['Punta \n Eugenia', 'San Quintín', 'Ensenada', 'San Diego', 'Los Angeles', 'Santa Barbara',
                        'Point Conception', 'Port San Luis']):
	ax.plot(kk, ll, 'o', color='r', markeredgecolor='k',zorder=5)
	if lab == 'Point Conception':
		ax.text(kk - 0.06, ll + 0.25, lab)
	elif lab == 'Santa Barbara':
		ax.text(kk + 0.2, ll - 0.05, lab)
	elif lab == 'Port San Luis':
		ax.text(kk + 0.2, ll - 0.1, lab)
	elif lab == 'Punta \n Eugenia':
		ax.text(kk + 0.6, ll - 0.1, lab, horizontalalignment='center')
	else:
		ax.text(kk + 0.16, ll - 0.05, lab)

ax.set_title(f'After {tt:.1f} hours')
ax.set_xlim(238 - 360, 246 - 360)
ax.set_ylim(27, 35.3)    
cbar = plt.colorbar(cax)
cbar.set_label('Energy Flux [Wm$^{-1}$]')

plt.show()
plt.savefig('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + 'EnergyFluxCshore.png')

#anim = FuncAnimation(fig, animateInternalWave,frames=56, repeat=False)

#anim.save('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + 'EnergyFluxIW.mp4', writer=writer, dpi=600)



