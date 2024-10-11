import SVBfunc 

import xarray as xr
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean
from matplotlib.animation import FuncAnimation
from matplotlib import animation
from matplotlib.ticker import FormatStrFormatter
import ffmpeg
import pylab as pl
from math import radians, cos
from functools import partial
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

from xmitgcm import open_mdsdataset
import scipy.interpolate as sciint

def animate(tt,VAL,ik,ax,caxin,caxin1,caxin2):
    print(tt)
    for coll in caxin.collections: 
        caxin.collections.remove(coll)
        
    for colll in caxin1.collections: 
        colll.remove

    if caxin2 != 0:
        for collll in caxin2.collections: 
            collll.remove
        caxin2=ax[2].contourf(X[ik+2],Z[ik+2],np.ma.masked_array(VAL[ik+2][tt]*1e4,mask=mask[ik+2]) ,cmap=colormap,levels=levels)
        ax[2].set_title(f'At {corrinds[ik+2]}°N \n Hour {hour[ik+2][tt]:.1f}', fontdict={'fontsize': 20})

    caxin=ax[0].contourf(X[ik],Z[ik],np.ma.masked_array(VAL[ik][tt]*1e4,mask=mask[ik]) ,cmap=colormap,levels=levels)
    caxin1=ax[1].contourf(X[ik+1],Z[ik+1],np.ma.masked_array(VAL[ik+1][tt]*1e4,mask=mask[ik+1]) ,cmap=colormap,levels=levels)
    
    ax[0].set_title(f'At {corrinds[ik]}°N \n Hour {hour[ik][tt]:.1f}', fontdict={'fontsize': 20})
    ax[1].set_title(f'At {corrinds[ik+1]}°N \n Hour {hour[ik+1][tt]:.1f}', fontdict={'fontsize': 20})

Writer = animation.writers['ffmpeg']
writer = Writer(fps=4, metadata=dict(artist='AM'), bitrate=2000)

corrinds=[30.49,30.77,31.13,31.69,32.11,32.65,33.02]

hej=[35,54,79,120,154,194,219]


coast='smooth'
var='ashore'
clabelB='Alongshore Velocity [cms$^{-1}$]' 
colormap=cmocean.cm.balance

file = '/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/Locations/' + str(var) 

indicing=np.arange(0,648,3)

mask=[]
lat=[]
lon=[]
Z=[]
hour=[]
X=[]
VAL=[]

for latind in corrinds:
    ds = xr.open_dataset(file + str(latind) + '.nc')
    VALin=ds.ashore[indicing]
    Xin=ds.X
    Zin=ds.Y


    latin=ds.LAT
    lonin=ds.LON

    hourin=ds.TIME.values[indicing]/(60*60)
    
    maski=np.ma.masked_values(VALin[0].values, 0)
    maskin = np.ma.getmask(maski)
    
    mask.append(maskin)
    hour.append(hourin)
    lat.append(latin)
    lon.append(lonin)
    Z.append(Zin)
    X.append(Xin)
    VAL.append(VALin)


params = {'font.size': 20,
          'figure.figsize': (20, 8),
         'font.family':'sans'}
pl.rcParams.update(params)

t=0

xlab='Cross-shore distance [km]'
ylab='Depth [m]'


for ik in [0,3,5]:
    plt.close()
    fig = plt.figure()
    if ik==0:
        gs = GridSpec(nrows=1, ncols=3, width_ratios=[0.9,0.9,1.1],hspace=0.35)
    else: 
        gs = GridSpec(nrows=1, ncols=2, width_ratios=[0.9,1.1],hspace=0.35)
        
    vmin=-np.nanmax(abs(VAL[ik][24]))*1e4
    vmax=np.nanmax(abs(VAL[ik][24]))*1e4

    levels=np.linspace(vmin,vmax,15)

    ax = fig.add_subplot(gs[0,0])

    ax.set_facecolor('tan')
    cax=ax.contourf(X[ik],Z[ik],np.ma.masked_array(VAL[ik][t]*1e4,mask=mask[ik]) ,cmap=colormap,levels=levels)

    ax.set_title(f'At {corrinds[ik]}°N \n Hour {hour[ik][t]:.1f}', fontdict={'fontsize': 20})

    ax.set_ylim([-1000,0])

    ax.set(ylabel=ylab,xlabel=xlab)

    ax1 = fig.add_subplot(gs[0,1])

    ax1.set_facecolor('tan')
    cax1=ax1.contourf(X[ik+1],Z[ik+1],np.ma.masked_array(VAL[ik+1][t]*1e4,mask=mask[ik+1]) ,cmap=colormap,levels=levels)

    ax1.set_title(f'At {corrinds[ik+1]}°N \n Hour {hour[ik+1][t]:.1f}', fontdict={'fontsize': 20})

    ax1.set_ylim([-1000,0])

    ax1.set(xlabel=xlab)
    ax1.yaxis.set_tick_params(labelleft=False)

    if ik !=0:
        cbar_ax = plt.colorbar(cax1)
        cbar_ax.ax.set_ylabel(clabelB)
        cbar_ax.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        cax2=0
        arrax=[ax,ax1]
    else:
        ax2 = fig.add_subplot(gs[0,2])

        ax2.set_facecolor('tan')
        cax2=ax2.contourf(X[ik+2],Z[ik+2],np.ma.masked_array(VAL[ik+2][t]*1e4,mask=mask[ik+2]) ,cmap=colormap,levels=levels)

        ax2.set_title(f'At {corrinds[ik+2]}°N \n Hour {hour[ik+2][t]:.1f}', fontdict={'fontsize': 20})

        ax2.set_ylim([-1000,0])

        ax2.set(xlabel=xlab)
        ax2.yaxis.set_tick_params(labelleft=False)


        cbar_ax = plt.colorbar(cax2)
        cbar_ax.ax.set_ylabel(clabelB)
        cbar_ax.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        
        arrax=[ax,ax1,ax2]

    anim = FuncAnimation(fig, partial(animate, VAL=VAL, ik=ik,ax=arrax,caxin=cax,caxin1=cax1,caxin2=cax2),frames=len(hour[0]), repeat=False)


    anim.save('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + '/' + str(var) + str(corrinds[ik]) +'.mp4', writer=writer, dpi=600)
