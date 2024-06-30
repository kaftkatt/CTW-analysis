import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import xarray as xr
import pylab as pl
from SVBfunc_Figures import plot_HOVMOLLER, plot_Map


eta=dsETA.ValfiltAll.values
dist=dsETA.dist.values
TIME=dsETA.time2.astype(int).values*1e-9
lon_ac=dsETA.lonAC.values
lat_ac=dsETA.latAC.values

pathETA='/home/athelandersson/NETCDFs/' + str(coast) + '/ETANAC.nc'
dsETA= xr.open_dataset(pathETA)
pathVEL='/home/athelandersson/NETCDFs/' + str(coast) + '/WVELAC.nc'
dsVEL= xr.open_dataset(pathVEL)

if coast == 'smooth':
	eta=dsETA.ValfiltAll.values[72:]
	WVEL=dsVEL.Valfilt.values[72:]
	TIME=dsETA.time.values[72:]/60
else: 
	eta=dsETA.ValfiltAll.values
	WVEL=dsVEL.Valfilt.values
	TIME=dsETA.time.values

distETA=dsETA.dist.values
lon_ac=dsETA.lonAC
lat_ac=dsETA.latAC

WVEL=dsVEL.ValfiltAll.values
distVEL=dsVEL.dist.values
lat_acVEL=dsVEL.latAC.values
lon_acVEL=dsVEL.lonAC.values

varname='phiHyd'
i=0
pathn='/home/athelandersson/NETCDFs/smooth_NO/'+ str(varname)+'noSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
pathw='/home/athelandersson/NETCDFs/smooth/'+ str(varname)+'withSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
        
dsw  = xr.open_dataset(pathw)
dsn = xr.open_dataset(pathn)


LAT=dsw.YC
LON=dsw.XC-360
Z = dsw.Z
hFacC = dsw.hFacC

hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)
        
depth=dsw.Depth
depthno=dsn.Depth

params = {'font.size': 26,
          'figure.figsize': (25, 8),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300





fig = plt.figure()
gs = GridSpec(nrows=1, ncols=3)

vmin=-5
vmax=5
cbarall=0

ax0 = fig.add_subplot(gs[0, 0])

plot_HOVMOLLER(ax0,distVEL,TIME,WVEL*1e6,'','Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,1,cbarall,'(a)')


vmin=-0.2
vmax=0.2
ax1 = fig.add_subplot(gs[0, 1])
plot_HOVMOLLER(ax1,dist,TIME,eta*1000,'','SSH [mm]',vmin,vmax,fig,lat_ac,lon_ac,0,cbarall,'(b)')


ax2 = fig.add_subplot(gs[0,2])
plotMap(ax2,LON,LAT,depth,mask,fig,'(c)')


plt.savefig('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + '/Fig3_HOV.png')

