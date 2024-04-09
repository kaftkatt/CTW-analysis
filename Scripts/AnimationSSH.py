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

coast='both'

#SSH
if coast == 'both':
        dirn = '/home/athelandersson/NETCDFs/original_NO/'
        dirw = '/home/athelandersson/NETCDFs/original/'
        dsw, dsn = SVBfunc.loadNetCDFs(dirw, dirn, 'dynVars',2)
        pathETA='/home/athelandersson/NETCDFs/original/ETA.nc'
        dsETA= xr.open_dataset(pathETA)
        pathETASm='/home/athelandersson/NETCDFs/smooth/ETA.nc'
        dsETASm= xr.open_dataset(pathETASm)
else:
        dirn = '/home/athelandersson/NETCDFs/' +  str(coast) + '_NO/'
        dirw = '/home/athelandersson/NETCDFs/' +  str(coast) + '/'
        dsw, dsn = SVBfunc.loadNetCDFs(dirw, dirn, 'dynVars',2)
        pathETA='/home/athelandersson/NETCDFs/' +  str(coast) + '/ETA.nc'
        dsETA= xr.open_dataset(pathETA)


if coast == 'original':
	etafiltall=dsETA.VALfilt.values
	TIME=dsETA.time.values
elif coast == 'smooth':
	etafiltall=dsETA.VALfilt.values[72:]
	TIME=dsETA.time.values[72:]
else:
	etafiltallOrig=dsETA.VALfilt.values
	TIME=dsETA.time.values
	etafiltallSm=dsETASm.VALfilt.values[72:]
	etafiltall=etafiltallSm-etafiltallOrig

LON=dsETA.x.values
LAT=dsETA.y.values
depth=dsw[0].Depth.values
hFacC = dsn[0].hFacC.values
hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)


# In[30]:


def animateETA(t):
    tin=t
    tt=(tin*20+2880)/60
    dep=0
    vmin=-0.2
    vmax=0.2
    print(tin)
    eta = etafiltall[tin,:,:]
    

    cax.set_array(np.ma.masked_array(eta*1000,mask=mask[dep,:,:]))
    ax.set_title(f'After {tt:.1f} hours')


# In[31]:


params = {'font.size': 22,
          'figure.figsize': (12, 8),
         'font.family':'sans'}
pl.rcParams.update(params)


Writer = animation.writers['ffmpeg']
writer = Writer(fps=4, metadata=dict(artist='AM'), bitrate=2000)


# In[34]:


fig, ax = plt.subplots()
    
vmin=-0.2
vmax=0.2
t=0
tt=((t*20)+2880)/60
xlab='Longitude [°]'
ylab='Latitude [°]'

ax.set_facecolor('wheat')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(etafiltall[0,:,:]*1000, mask=mask[0,:,:]),cmap=cmocean.cm.curl,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth,  colors=['0.2'], 
                levels=[0])
ax.set(xlabel=xlab, ylabel=ylab)

ax.set_title(f'After {tt:.1f} hours')
ax.set_xlim(-122,-114)     
cbar = plt.colorbar(cax)
cbar.set_label('SSH [mm]')
ax.set_ylim(27,35.3)

anim = FuncAnimation(fig, animateETA,frames=575, repeat=False)

    
anim.save('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + '/SSH.mp4', writer=writer, dpi=600)






