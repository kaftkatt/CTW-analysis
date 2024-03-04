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


# In[3]:
dirn = '/media/amelia/Trillian/SVB/exp06_512x612x100_ORL/01_noSVB_febTS/'
dirw = '/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01_SVB_febTS/'

dsw, dsn = SVBfunc.loadNetCDFs(dirw, dirn, 'DYNVARS')

Writer = animation.writers['ffmpeg']
writer = Writer(fps=4, metadata=dict(artist='AM'), bitrate=2000)

params = {'font.size': 22,
          'figure.figsize': (12, 8),
         'font.family':'sans'}
pl.rcParams.update(params)


def animate(t):
    t=t+100
    tt=(t*20+2880)/60
    dep=55
    vmin=-0.000002
    vmax=0.000002
    print(t)
    W = get_snapshot_at_level( t,dep,dsw,dsn)
    

    cax.set_array(np.ma.masked_array(W,mask=maskw[dep,:,:]))
    ax.set_title(f'At depth {Z[dep].values:.2f} m. After {tt:.1f} hours')



#Index to call from the list of netcdfs
ind=1 #0 is day 2-3, 1 is day 3-4 until index 7 (day 9-10)
dep=55  #483.2 meter depth is the 55th element
t=100-72
tt=(((72*ind+t)*20)+2880)/60 # Gives amount of hours from start of the model, starts at hour 48 if ind=0 and t=0


Ww=dsw[ind].WVEL
Wn=dsn[ind].WVEL
Win=Ww[t,dep,:,:]-Wn[t,dep,:,:]
LON=dsw[ind].XC-360
LAT=dsw[ind].YC
Z=dsw[ind].Zl

hFacCw = dsn[ind].hFacC
hFacCusew=hFacCw.values

hfa = np.ma.masked_values(hFacCusew, 0)
maskw = np.ma.getmask(hfa)


# In[29]:


fig, ax = plt.subplots()

vmin=-0.000002
vmax=0.000002

xlab='Longitude [°]'
ylab='Latitude [°]'

ax.set_facecolor('wheat')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(Win,mask=maskw[dep,:,:]),
                    cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth,  colors=['0.2','0.4'], 
                levels=[0,500])
ax.set(xlabel=xlab, ylabel=ylab)

ax.set_title(f'At depth {Z[dep].values:.2f} m. After {tt:.1f} hours')
ax.set_xlim(-122,-114)     
cbar = plt.colorbar(cax)
cbar.set_label('Vertical velocity [m/s]')
ax.set_ylim(27,35.3)

anim = FuncAnimation(fig, animate,frames=475, repeat=False)

    
anim.save('WVEL_attempt3.mp4', writer=writer, dpi=600)


# ## SSH

# In[9]:


pathETA='/media/amelia/Trillian/SVB/ETA.nc'
dsETA= xr.open_dataset(pathETA)
depth=dsw[0].Depth.values


# In[10]:


etafiltall=dsETA.ETAfiltall.values
LON=dsETA.x.values
LAT=dsETA.y.values
TIME=dsETA.time.values
hFacC = dsn[0].hFacC.values
hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)


# In[30]:


def animateETA(t):
    tin=t+100
    tt=(tin*20+2880)/60
    dep=0
    vmin=-0.01*10
    vmax=0.01*10
    print(tin)
    eta = etafiltall[tin,:,:]
    

    cax.set_array(np.ma.masked_array(eta*1000,mask=mask[dep,:,:]))
    ax.set_title(f'After {tt:.1f} hours')



# In[31]:


params = {'font.size': 22,
          'figure.figsize': (12, 8),
         'font.family':'sans'}
pl.rcParams.update(params)


# In[8]:


np.shape(etafiltall)


# In[32]:


Writer = animation.writers['ffmpeg']
writer = Writer(fps=4, metadata=dict(artist='AM'), bitrate=2000)


# In[34]:


fig, ax = plt.subplots()
    
vmin=-0.01*10
vmax=0.01*10
t=100
tt=((t*20)+2880)/60
xlab='Longitude [°]'
ylab='Latitude [°]'

ax.set_facecolor('wheat')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(etafiltall[100,:,:]*1000, mask=mask[0,:,:]),cmap=cmocean.cm.curl,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth,  colors=['0.2'], 
                levels=[0])
ax.set(xlabel=xlab, ylabel=ylab)

ax.set_title(f'After {tt:.1f} hours')
ax.set_xlim(-122,-114)     
cbar = plt.colorbar(cax)
cbar.set_label('SSH [mm]')
ax.set_ylim(27,35.3)

anim = FuncAnimation(fig, animateETA,frames=475, repeat=False)

    
anim.save('SSHfeb.mp4', writer=writer, dpi=600)





