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


coast='smooth'
tstart=2

dirn = '/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
dirw = '/home/athelandersson/NETCDFs/' + str(coast) + '/'

dsw, dsn = SVBfunc.loadNetCDFs(dirw, dirn, 'dynVars',tstart)



def animate(t):
    t=t
    tt=(t*20+2880)/60
    dep=55
    vmin=-0.000002
    vmax=0.000002
    print(t)
    W = SVBfunc.get_snapshot_at_level( t,dep,dsw,dsn)
    

    cax.set_array(np.ma.masked_array(W,mask=maskw[dep,:,:]))
    ax.set_title(f'At depth {Z[dep].values:.2f} m. After {tt:.1f} hours')


# In[5]:


Writer = animation.writers['ffmpeg']
writer = Writer(fps=4, metadata=dict(artist='AM'), bitrate=2000)


# In[27]:


params = {'font.size': 22,
          'figure.figsize': (12, 8),
         'font.family':'sans'}
pl.rcParams.update(params)


# In[7]:


#Index to call from the list of netcdfs
ind=0 #0 is day 2-3, 1 is day 3-4 until index 7 (day 9-10)
dep=55  #483.2 meter depth is the 55th element
t=0
tt=(((72*ind+t)*20)+2880)/60 # Gives amount of hours from start of the model, starts at hour 48 if ind=0 and t=0


Ww=dsw[ind].WVEL
Wn=dsn[ind].WVEL
Win=Ww[t,dep,:,:]-Wn[t,dep,:,:]
LON=dsw[ind].XC-360
LAT=dsw[ind].YC
Z=dsw[ind].Zl

depth=dsw[0].Depth.values

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

anim = FuncAnimation(fig, animate,frames=575, repeat=False)

    
anim.save('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + 'WVEL.mp4', writer=writer, dpi=600)



