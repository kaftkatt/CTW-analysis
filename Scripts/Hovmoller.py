import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as pl
import scipy.signal as sig
from scipy.io import loadmat

from SVBfuncPlotting import plot_HOVMOLLER, closest

coasto='original'
pathVELo='/home/athelandersson/NETCDFs/' + str(coasto) + '/WVELAC.nc'
dsVELo= xr.open_dataset(pathVELo)
WVELo=dsVELo.Valfilt.values
TIMEVELo=dsVELo.time.values*60 #To make it in seconds

coast='smooth'
pathVEL='/home/athelandersson/NETCDFs/' + str(coast) + '/WVELAC.nc'
dsVEL= xr.open_dataset(pathVEL)

tstart=72
WVEL=dsVEL.Valfilt.values[tstart:]

#To make zero right outside of the bay and not inside the bay
distVEL=dsVEL.dist.values-dsVEL.dist.values[0]

TIMEVEL=dsVEL.time.values[tstart:]/60

lat_acVEL=dsVEL.latAC.values
lon_acVEL=dsVEL.lonAC.values

matfile=loadmat( '/home/athelandersson/CTW-analysis/Files/smooth/BT_PALL_MovAv2025.mat')
latout=matfile['lat'][0][:11]
lonout=matfile['lon'][0][:11]
dep=matfile['d'][0][:11]

#Define the latitude at which the crossection has been taken, corresponding to locations along the coast
lonin=[]
latin=[]
for i in range(len(lonout)):
    lonin.append(lonout[i][0][0])
    latin.append(latout[i][0][0])

#Define latitudes for citites along the coast
ind_lon_cities = [ -115.939167, -116.605833, -117.1625]
ind_lat_cities = [ 30.556389, 31.857778, 32.715]

ind_lat,ind_lon=closest(lat_acVEL,ind_lat_cities,lon_acVEL,ind_lon_cities)

#Find the actual latitude for the current depth corresponding to the latitudes of the chosen crossections
latloc=[]
lonloc=[]

for i in range(len(dep)):
    ind=np.where(dep[i][0]>=500)[0][0]
    latloc.append(latout[i][0][ind])
    lonloc.append(lonout[i][0][ind])

#Find the corresponding index for the values along the coast
loclatIn,loclonIn=closest(lat_acVEL,latloc,lon_acVEL,lonloc)

#Find the last element before the islands
end=np.where(lat_acVEL>latin[-1])[0][0]

fig = plt.figure()
gs = GridSpec(nrows=1, ncols=2)

ax = fig.add_subplot(gs[0, 0])
vmin=-5
vmax=5
cbarall=1

plot_HOVMOLLER(ax,distVEL[:end],TIMEVEL,WVEL[:,:end]*1e6,'Smoothened','Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,1,cbarall,'(a)')

for ll, lab in zip(ind_lat,
                       ['San Quintín', 'Ensenada', 'San Diego']):
    ax.plot(0, distVEL[ll], "_",markersize=50, color='white',zorder=5)
    ax.text(0, distVEL[ll]+5, lab, fontsize=25,color='white')

ax = fig.add_subplot(gs[0, 1])

plot_HOVMOLLER(ax,distVEL[:end],TIMEVELo,WVELo[:,:end]*1e6,'Original','Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,0,cbarall,'(b)')
for ll, lab in zip(ind_lat,
                       ['San Quintín', 'Ensenada', 'San Diego']):
    ax.plot(0, distVEL[ll], "_",markersize=50, color='white',zorder=5)
    ax.text(0, distVEL[ll]+5, lab, fontsize=25,color='white')

fig.tight_layout()
plt.savefig('/home/athelandersson/CTW-analysis/Figures/Article/Hovmoller.png')
