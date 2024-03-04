from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import xarray as xr
import cmocean
import pylab as pl
from scipy.io import loadmat
from math import radians, cos, sin, asin, sqrt
import warnings
warnings.filterwarnings("ignore")
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import SVBfunc


# In[2]:

hej=[58, 85, 205, 227]
pathETA='/media/amelia/Trillian/SVB/ETA.nc'
dsETA= xr.open_dataset(pathETA)
pathVEL='/media/amelia/Trillian/SVB/coastVEL.nc'
dsVEL= xr.open_dataset(pathVEL)
dswAC= xr.open_dataset("/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01_SVB_febTS/ETAwithSVBACall.nc")


# In[3]:


eta=dsETA.ETAfiltcoast.values
distETA=dsETA.dist.values
TIME=dsETA.time.values
lon_ac=dswAC.lonAC
lat_ac=dswAC.latAC


# In[4]:


WVEL=dsVEL.Wfilt.values
dist=dsVEL.dist.values
lat_acVEL=dsVEL.lat_ac.values
lon_acVEL=dsVEL.lon_ac.values


# In[5]:


varname='PHIHYD'
i=0
pathn='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL/01b_noSVB_febTS/'+ str(varname)+'noSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
pathw='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01b_SVB_febTS/'+ str(varname)+'withSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
        
dsw  = xr.open_dataset(pathw)
dsn = xr.open_dataset(pathn)


# In[6]:


LAT=dsw.YC
LON=dsw.XC-360
Z = dsw.Z
hFacC = dsw.hFacC

hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)
        
depth=dsw.Depth
depthno=dsn.Depth

matfile=loadmat('BT_Perp.mat')
x,dep,indXlon,indYlat=matfile['dist'],matfile['d'],matfile['indexXlon'],matfile['indexYlat']




params = {'font.size': 16,
          'figure.figsize': (11, 12),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


# In[15]:


fig = plt.figure()
gs = GridSpec(nrows=2, ncols=2)#, height_ratios=[1, 1,1])

vmin=-5
vmax=5
cbarall=0

ax0 = fig.add_subplot(gs[0, 0])
SVBfunc.plot_HOVMOLLER(ax0,dist,TIME,WVEL*1e6,'',
               'Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,1,cbarall,'(a)')




vmin=-0.2
vmax=0.2
ax1 = fig.add_subplot(gs[0, 1])
SVBfunc.plot_HOVMOLLER(ax1,distETA,TIME,eta*1000,'',
               'SSH [mm]',vmin,vmax,fig,lat_acVEL,lon_acVEL,0,cbarall,'(b)')




ax2 = fig.add_subplot(gs[1,0])
ax2.set_facecolor('tan')

pc = ax2.contourf(LON,LAT,np.ma.masked_array(depth, mask=mask[0,:,:]),50,
                 vmin=0, vmax=5000, cmap=cmocean.cm.deep)


cn = ax2.contour(LON,LAT,depth, colors=['0.2','0.4'], 
                levels=[0,500])
divider = make_axes_locatable(ax2)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(pc, cax=axdiv,orientation='horizontal',ticks=np.arange(0,np.max(depth),1000))
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Depth [m]')


ax2.set_xlabel('Lon [°]')
ax2.set_ylabel('Lat [°]')
ax2.set_xlim(238-360, 246-360)
ax2.set_ylim(27,35.3)
ax2.set_aspect(1)
ax2.text(-0.1, 1.2, '(c)', fontweight='bold', color='k', 
        transform=ax2.transAxes)

markers=Line2D.filled_markers
markers=np.delete(markers,np.arange(2,5,1))
colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999']


for i in np.arange(0,6,1):
    #ax0.plot(TIME[0],dist[107*i],'o',color='#dede00',markersize=15)
    #ax1.plot(TIME[0],distETA[107*i],'o',color='#984ea3',markersize=15)
    ax2.plot(LON[lon_acVEL[107*i]],LAT[lat_acVEL[107*i]],'ok')
    ax2.plot(LON[lon_ac[107*i]],LAT[lat_ac[107*i]],'or')
    
    #ax0.annotate('Where SVB is',xy=(0.93,0.4), xytext=(0.7, 0.15),xycoords='axes fraction')
    #ax2.plot(LON[lon_acVEL[53*(1+2*i)]],LAT[lat_acVEL[53*(1+2*i)]],marker=markers[0], 
           #markersize=7, color='k')
    

p=0
for i in hej:
    p=p+1
    ax0.axhline(y=dist[i],color=colors[p-1],linewidth=2)
    ax1.axhline(y=distETA[np.where(lat_ac==lat_acVEL[i])[0][0]],color=colors[p-1],linewidth=2)
    ax2.plot(LON[indXlon[p-1]],LAT[indYlat[p-1]],color=colors[p-1],linewidth=2)


ax3 = fig.add_subplot(gs[1,1])
for i in np.arange(0,4,1):
    ax3.plot(x[i],-dep[i], label=f'{dist[hej[i]]:.0f} km', color=colors[i])


ax3.legend()
ax3.set_xlabel('Distance [km]')
ax3.set_ylabel('Depth [m]')
ax3.text(-0.08, 1.02, '(d)', fontweight='bold', color='k', 
        transform=ax3.transAxes)
fig.tight_layout()


# ## For the presentation

# In[11]:


params = {'font.size': 26,
          'figure.figsize': (25, 8),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


# In[12]:


fig = plt.figure()
gs = GridSpec(nrows=1, ncols=3)#, height_ratios=[1, 1,1])

vmin=-5
vmax=5
cbarall=0

ax0 = fig.add_subplot(gs[0, 0])
plot_HOVMÖLLER(ax0,dist,TIME,WVEL*1e6,'',
               'Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,1,cbarall)

ax0.text(-0.15, 1.2, '(a)', fontweight='bold', color='k', 
        transform=ax0.transAxes)


vmin=-0.2
vmax=0.2
ax1 = fig.add_subplot(gs[0, 1])
plot_HOVMÖLLER(ax1,distETA,TIME,eta*1000,'',
               'SSH [mm]',vmin,vmax,fig,lat_acVEL,lon_acVEL,0,cbarall)

ax1.text(-0.1, 1.2, '(b)', fontweight='bold', color='k', 
        transform=ax1.transAxes)


ax2 = fig.add_subplot(gs[0,2])
ax2.set_facecolor('tan')

pc = ax2.contourf(LON,LAT,np.ma.masked_array(depth, mask=mask[0,:,:]),50,
                 vmin=0, vmax=5000, cmap=cmocean.cm.deep)


cn = ax2.contour(LON,LAT,depth, colors=['0.2','0.4'], 
                levels=[0,500])
divider = make_axes_locatable(ax2)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(pc, cax=axdiv,orientation='horizontal',ticks=np.arange(0,np.max(depth),1000))
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Depth [m]')


ax2.set_xlabel('Lon [°]')
ax2.set_ylabel('Lat [°]')
ax2.set_xlim(238-360, 246-360)
ax2.set_ylim(27,35.3)
ax2.set_aspect(1)
ax2.text(-0.1, 1.2, '(c)', fontweight='bold', color='k', 
        transform=ax2.transAxes)

markers=Line2D.filled_markers
markers=np.delete(markers,np.arange(2,5,1))
colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999']


for i in np.arange(0,6,1):
    #ax0.plot(TIME[0],dist[107*i],'o',color='#dede00',markersize=15)
    #ax1.plot(TIME[0],distETA[107*i],'o',color='#984ea3',markersize=15)
    ax2.plot(LON[lon_acVEL[107*i]],LAT[lat_acVEL[107*i]],'ok')
    ax2.plot(LON[lon_ac[107*i]],LAT[lat_ac[107*i]],'or')


# In[31]:


params = {'font.size': 26,
          'figure.figsize': (25, 9),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


# In[34]:


fig = plt.figure()
gs = GridSpec(nrows=1, ncols=3, width_ratios=[1,1,1])

vmin=-5
vmax=5
cbarall=0

ax0 = fig.add_subplot(gs[0, 0])
plot_HOVMÖLLER(ax0,dist,TIME,WVEL*1e6,'',
               'Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,1,cbarall)

ax0.text(-0.15, 1.2, '(a)', fontweight='bold', color='k', 
        transform=ax0.transAxes)



ax2 = fig.add_subplot(gs[0,1])
ax2.set_facecolor('tan')

pc = ax2.contourf(LON,LAT,np.ma.masked_array(depth, mask=mask[0,:,:]),50,
                 vmin=0, vmax=5000, cmap=cmocean.cm.deep)


cn = ax2.contour(LON,LAT,depth, colors=['0.2','0.4'], 
                levels=[0,500])
divider = make_axes_locatable(ax2)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(pc, cax=axdiv,orientation='horizontal',ticks=np.arange(0,np.max(depth),1000))
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Depth [m]')


ax2.set_xlabel('Lon [°]')
ax2.set_ylabel('Lat [°]')
ax2.set_xlim(238-360, 246-360)
ax2.set_ylim(27,35.3)
ax2.set_aspect(1)
ax2.text(-0.1, 1.2, '(b)', fontweight='bold', color='k', 
        transform=ax2.transAxes)

markers=Line2D.filled_markers
markers=np.delete(markers,np.arange(2,5,1))
colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999']


for i in np.arange(0,6,1):
    #ax0.plot(TIME[0],dist[107*i],'o',color='#dede00',markersize=15)
    #ax1.plot(TIME[0],distETA[107*i],'o',color='#984ea3',markersize=15)
    ax2.plot(LON[lon_acVEL[107*i]],LAT[lat_acVEL[107*i]],'ok')
    ax2.plot(LON[lon_ac[107*i]],LAT[lat_ac[107*i]],'or')
p=0
for i in hej:
    p=p+1
    ax0.axhline(y=dist[i],color=colors[p-1],linewidth=2)
    ax2.axhline(y=LAT[lat_acVEL[i]],color=colors[p-1],linewidth=2)


ax3 = fig.add_subplot(gs[0,2])
for i in np.arange(0,4,1):
    ax3.plot(distcross[i,depcross[i,:]!=0],-depcross[i,depcross[i,:]!=0], label=f'{dist[hej[i]]:.0f} km', color=colors[i])
ax3.legend()
ax3.set_xlabel('Distance [km]')
ax3.set_ylabel('Depth [m]')

ax3.text(-0.1, 1.05, '(c)', fontweight='bold', color='k', 
        transform=ax3.transAxes)


fig.tight_layout()


# In[62]:


params = {'font.size': 22,
          'figure.figsize': (10,10),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


# In[63]:


ax3 = plt.subplot()
for i in np.arange(0,4,1):
    ax3.plot(distcross[i,depcross[i,:]!=0],-depcross[i,depcross[i,:]!=0], label=f'{LAT.values[lat_acVEL[hej[i]]]:.2f} °N', color=colors[i])
ax3.legend()
ax3.set_xlabel('Distance [km]')
ax3.set_ylabel('Depth [m]')

ax3.text(-0.1, 1.05, '(c)', fontweight='bold', color='k', 
        transform=ax3.transAxes)
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")

fig.tight_layout()


# In[ ]:




