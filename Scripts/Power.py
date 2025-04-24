import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as pl
import scipy.signal as sig
from scipy.io import loadmat

from SVBfuncPlotting import FFRQ, closest

coasto='original'
pathVELo='/home/athelandersson/NETCDFs/' + str(coasto) + '/WVELAC.nc'
dsVELo= xr.open_dataset(pathVELo)
WVELo=dsVELo.Valfilt.values
TIMEVELo=dsVELo.time.values*60 #To make it in seconds

coast='smooth'
pathVEL='/home/athelandersson/NETCDFs/' + str(coast) + '/WVELAC.nc'
dsVEL= xr.open_dataset(pathVEL)

WVEL=dsVEL.Valfilt.values

#To make zero right outside of the bay and not inside the bay
distVEL=dsVEL.dist.values-dsVEL.dist.values[0]

TIMEVEL=dsVEL.time.values

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

#Perform the fast fourier transform 
psdfilt, freqfilt=FFRQ(WVEL,TIMEVEL,distVEL)
psdfilto, freqfilto=FFRQ(WVELo,TIMEVELo,distVEL)

# Define earths rotation
omega=(2*np.pi)/(23*60*60 + 56*60 + 4.1)

#### Begin Plotting ####
params = {'font.size': 24,
          'figure.figsize': (20, 16),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

fig = plt.figure()
gs = GridSpec(nrows=2, ncols=5, height_ratios=[1, 0.8])

const=1e6
vmin=0
vmax=3

title='PSD [10$^{-6}$ ms$^{-1}$Hz$^{-1}$]'
xlab='Frequency [cpd]'
ylab='Distance [km]'

ax = fig.add_subplot(gs[0, 0:2])


ax.text(-0.08, 1.02, '(a)', fontweight='bold', color='k', 
        transform=ax.transAxes)

cax = ax.pcolormesh(freqfilt[1,:]*(24*3600),distVEL[:end],psdfilt[:end,:]*const,vmin=vmin,vmax=vmax, cmap=cmocean.cm.balance)

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.9)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label(title)
ax.set_xlim((0,0.000016*(24*3600)))
ax.set(ylabel=ylab)

for ll, lab in zip(ind_lat,
                       ['San Quintín', 'Ensenada', 'San Diego']):
    ax.plot(0, distVEL[ll], "_",markersize=50, color='white',zorder=5)
    ax.text(0, distVEL[ll]+5, lab, fontsize=25,color='white')
    

ax.set_title('Smoothened')
    
ax1 = fig.add_subplot(gs[1, 0:2])
ax1.set(xlabel=xlab, ylabel=ylab)
ax1.text(-0.08, 1.05, '(b)', fontweight='bold', color='k', 
        transform=ax1.transAxes)

cax1 = ax1.pcolormesh(freqfilto[1,:]*(24*3600),distVEL[:end],psdfilt[:end,:]*const,vmin=vmin,vmax=vmax, cmap=cmocean.cm.balance)

ax1.set_xlim((0,0.000016*(24*3600)))

ax1.set(xlabel=xlab, ylabel=ylab)

ax1.set_title('Original')

for ll, lab in zip(ind_lat,
                       ['San Quintín', 'Ensenada', 'San Diego']):
    ax1.plot(0, distVEL[ll], "_",markersize=50, color='white',zorder=5)
    ax1.text(0, distVEL[ll]+5, lab, fontsize=25,color='white')

axin = fig.add_subplot(gs[:,2:])

for nr in np.arange(0,len(loclatIn),2):
    axin.axhline(nr,color='red',zorder=1)
    if nr == 0:
        axin.plot((freqfilt[1])*(24*3600),(psdfilt[loclatIn[nr]]*const)+nr,c='k',linewidth=2,zorder=20,label='Smoothened') 
        axin.plot((freqfilto[1])*(24*3600),(psdfilto[loclatIn[nr]]*const)+nr,alpha=0.6,c='k',linewidth=2,linestyle='dashed',zorder=10,label='Original') 
        axin.scatter((2*omega*np.sin(np.deg2rad(lat_acVEL[loclatIn[nr]]))*24*3600)/(2*np.pi),nr,color='blue',marker="|",s=1000,label='Inertial Frequency')
    else:
        axin.plot((freqfilt[1])*(24*3600),(psdfilt[loclatIn[nr]]*const)+nr,c='k',linewidth=2,zorder=20) 
        axin.scatter((2*omega*np.sin(np.deg2rad(lat_acVEL[loclatIn[nr]]))*24*3600)/(2*np.pi),nr,color='blue',marker="|",s=1000)
        axin.plot((freqfilto[1])*(24*3600),(psdfilto[loclatIn[nr]]*const)+nr,alpha=0.6,c='k',linewidth=2,linestyle='dashed',zorder=10) 
        ax1.axhline(distVEL[loclatIn[nr]],color='red',linestyle='dashed',linewidth=4)
        ax.axhline(distVEL[loclatIn[nr]],color='red',linestyle='dashed',linewidth=4)
    

axin.axvline((freqfilt[1][np.argmax(psdfilt[loclatIn[nr]])])*(24*3600),linestyle='dotted',color='k')
axin.legend()
xlab='Frequency [cpd]'
ylab='PSD [$10^{-4}$  ms$^{-1}$Hz$^{-1}$]'

axin.set_xlim((0,0.000016*(24*3600)))

axin.set(xlabel=xlab, ylabel=title)
axin.set_xticks(freqfilt[1][:15]*(24*3600))

axin.minorticks_on()

axin.grid(which='minor',linestyle='--', alpha=0.5)
axin.grid(which='major')

axin.text(-0.08, 1, '(c)', fontweight='bold', color='k', 
        transform=axin.transAxes)

fig.tight_layout()
plt.savefig('/home/athelandersson/CTW-analysis/Figures/Article/Power.png')

