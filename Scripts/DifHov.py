import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import pylab as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec

import SVBfunc
hej2=[35,54,79,120,154,194,219]

pathETA='/home/athelandersson/NETCDFs/original/ETANAC.nc'
ds= xr.open_dataset(pathETA)

pathETASm='/home/athelandersson/NETCDFs/smooth/ETANAC.nc'
dsSm= xr.open_dataset(pathETASm)

pathVEL='/home/athelandersson/NETCDFs/original/WVELAC.nc'
dsVEL= xr.open_dataset(pathVEL)

pathVELSm='/home/athelandersson/NETCDFs/smooth/WVELAC.nc'
dsVELSm= xr.open_dataset(pathVELSm)

WVELOrg=dsVEL.ValfiltAll.values
WVELSm=dsVELSm.ValfiltAll.values[72:,:len(WVELOrg[0,:])]
WVEL=WVELOrg-WVELSm

ETASm=dsSm.ValfiltAll.values[72:]
ETAOrg=ds.ValfiltAll.values[:,:len(ETASm[0,:])]

ETA=ETAOrg-ETASm

TIMEVEL=dsVEL.time2.values
distAC=dsSm.dist.values

distVEL=dsVEL.dist.values

lat_acVEL=dsVEL.latAC.values
lon_acVEL=dsVEL.lonAC.values

lon_ac=ds.lonAC.values
lat_ac=ds.latAC.values

params = {'font.size': 16,
          'figure.figsize': (10, 5),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999','tab:blue']

fig = plt.figure()
gs = GridSpec(nrows=1, ncols=2)


ax = fig.add_subplot(gs[0, 0])

vmin=-10
vmax=10
cbarall=0
SVBfunc.plot_HOVMOLLER(ax,distVEL,TIMEVEL,WVEL*1e6,'','Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,1,cbarall,'(a)')

for i in range(len(hej2)):
	ax.axhline(y=distAC[hej2[i]],color=colors[i],linewidth=2,alpha=0.7)

ax = fig.add_subplot(gs[0, 1])

vmin=-0.2
vmax=0.2
cbarall=0
SVBfunc.plot_HOVMOLLER(ax,distAC,TIMEVEL,ETA*1e3,'','SSH  [mm]',vmin,vmax,fig,lat_ac,lon_ac,1,cbarall,'(b)')

for i in range(len(hej2)):
	ax.axhline(y=distAC[hej2[i]],color=colors[i],linewidth=2,alpha=0.7)

fig.tight_layout()

plt.savefig('/home/athelandersson/CTW-analysis/Figures/CompOp.png')	
