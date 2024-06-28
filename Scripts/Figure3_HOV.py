from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import xarray as xr
import pylab as pl
import SVBfunc
from matplotlib.lines import Line2D
from scipy.io import loadmat

pathETA='ETAAC.nc'
dsETA= xr.open_dataset(pathETA)
pathVEL='WVELAC.nc'
dsVEL= xr.open_dataset(pathVEL)


eta=dsETA.ValfiltAll.values
dist=dsETA.dist.values
TIME=dsETA.time2.astype(int).values*1e-9
lon_ac=dsETA.lonAC.values
lat_ac=dsETA.latAC.values


# In[4]:


WVEL=dsVEL.ValfiltAll.values
distVEL=dsVEL.dist.values
TIMEVEL=dsVEL.time2.values
lat_acVEL=dsVEL.latAC.values
lon_acVEL=dsVEL.lonAC.values

varname='PHIHYD'
i=0
pathn='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL/01b_noSVB_febTS/'+ str(varname)+'noSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
pathw='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01b_SVB_febTS/'+ str(varname)+'withSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
        
dsw  = xr.open_dataset(pathw)
dsn = xr.open_dataset(pathn)

matfile=loadmat('BT_P.mat')
indXlon,indYlat=matfile['indexXlon'],matfile['indexYlat']
# In[6]:


LAT=dsw.YC
LON=dsw.XC-360
Z = dsw.Z
hFacC = dsw.hFacC

hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)
        
depth=dsw.Depth
depthno=dsn.Depth

params = {'font.size': 16,
          'figure.figsize': (11, 12),
         'font.family':'sans'}
pl.rcParams.update(params)

fig = plt.figure()
gs = GridSpec(nrows=2, ncols=2)

vmin=-5
vmax=5
cbarall=0

ax0 = fig.add_subplot(gs[0, 0])

SVBfunc.plot_HOVMOLLER(ax0,distVEL,TIMEVEL,WVEL*1e6,'','Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,1,cbarall,'(a)')



vmin=-0.2
vmax=0.2
ax1 = fig.add_subplot(gs[0, 1])
SVBfunc.plot_HOVMOLLER(ax1,dist,TIME,eta*1000,'','SSH [mm]',vmin,vmax,fig,lat_acVEL,lon_acVEL,0,cbarall,'(b)')


ax2 = fig.add_subplot(gs[1,0])
SVBfunc.plotMap(ax2,LON,LAT,depth,mask,fig,'(c)')

markers=Line2D.filled_markers
markers=np.delete(markers,np.arange(2,5,1))
colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999']

hej=[35,54,79,120,154,194,219]  
corrind=[30.49,30.77,31.13,31.69,32.11,32.65,33.02]  #30.706,31.1276,32.5,32.95
hejVEL=np.zeros(len(hej))
for i in range(len(hej)):
	hejVEL[i]=np.where(lat_acVEL==lat_ac[hej[i]])[0][0]
	#hej[i]=np.where(lat_ac==indYlat[i][0]+1)[0][0]

hejVEL=hejVEL.astype(int)

for i in range(len(hej)):
    ax0.axhline(y=distVEL[hejVEL[i]],color=colors[i],linewidth=2)
    ax1.axhline(y=dist[hej[i]],color=colors[i],linewidth=2)
    ax2.plot(LON[indXlon[i]].values,LAT[indYlat[i]].values,color=colors[i],linewidth=2)

for i in range(5):
    ax2.plot(LON[lon_acVEL[107*i]],LAT[lat_acVEL[107*i]],'ok')
    ax2.plot(LON[lon_ac[99*i]],LAT[lat_ac[99*i]],'or')


ax3 = fig.add_subplot(gs[1,1])
SVBfunc.plot_batylines(lat_ac,dist,LAT,ax3,'(d)',hej)


plt.show()
