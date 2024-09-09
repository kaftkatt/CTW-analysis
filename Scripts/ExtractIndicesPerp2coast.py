import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import pylab as pl
import scipy.signal as sig
import scipy.io as sio
import scipy.interpolate as sciint
from xmitgcm import open_mdsdataset
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
from os.path import exists

from scipy.io import savemat
import SVBfunc

from math import radians, cos, sin, asin, sqrt, atan, degrees, log

coast='smooth'
if coast == 'original':
	tstart=0
	indstart=2
elif coast == 'smooth':
	tstart=72
	indstart=1

dirn = '/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
dirw = '/home/athelandersson/NETCDFs/' + str(coast) + '/'
dsw, dsn = SVBfunc.loadNetCDFs(dirw, dirn, 'phiHyd',indstart)

pathETA='/home/athelandersson/NETCDFs/' + str(coast) + '/ETANAC.nc'
ds= xr.open_dataset(pathETA)

lon_ac=ds.lonAC.values
lat_ac=ds.latAC.values
distAC=ds.dist.values

hej2=[35,54,79,120,154,194,219]
#hej2=np.arange(1,len(lon_ac)-1,1)

LAT = dsw[0].YC
LON = dsw[0].XC - 360

Z = dsw[0].Z
hFacC = dsw[0].hFacC

hfa = np.ma.masked_values(hFacC[0, :, :], 0)
mask = np.ma.getmask(hfa)

depth = dsw[0].Depth
depthno = dsn[0].Depth

iniX=[]
iniY=[]
dist=[]
dep=[]

for ind in hej2:
	
	lon1=LON[lon_ac[ind-2]]
	lat1=LAT[lat_ac[ind-2]]
	lon2=LON[lon_ac[ind-2]]
	lat2=LAT[lat_ac[ind+2]]
	a=SVBfunc.haversine(lon1, lat1, lon2, lat2)
	print(LON[lon_ac[ind-2]].values)
	print(LON[lon_ac[ind+2]].values)
	print(LAT[lat_ac[ind-2]].values)
	print(LAT[lat_ac[ind+2]].values)
	lon3=LON[lon_ac[ind+2]]
	lat3=LAT[lat_ac[ind+2]]
	
	b=SVBfunc.haversine(lon2, lat2, lon3, lat3)
	
	deg=-atan(a/b)
	R1=LON*cos(deg)-LAT*sin(deg)
	R2=LON*sin(deg)+LAT*cos(deg)
	
	startX=R1.sel(XC=LON[lon_ac[ind]]+360,YC=LAT[lat_ac[ind]])
	startY=R2.sel(XC=LON[lon_ac[ind]]+360,YC=LAT[lat_ac[ind]])
	
	indexX,indexY=np.where(np.isnan(R2.where(R2==startY).values)==False)
	
	Rcoast = R1[indexX,indexY]
	
	LONIN=R1[:indexX[0],indexY]
	
	LATIN=R2[(np.ones(len(R1[:indexX[0],indexY]))*indexX).astype(int),(indexY).astype(int)]
	
	R1Back=np.flip(LONIN.values*cos(-deg)-LATIN.values*sin(-deg))
	R2Back=np.flip(LONIN.values*sin(-deg)+LATIN.values*cos(-deg))
	
	indlon=np.flip(np.where(np.logical_and(LON>min(R1Back),LON<max(R1Back)))[0])
	indlat=np.flip(np.where(np.logical_and(LAT>min(R2Back),LAT<max(R2Back)))[0])                 	
	
	Lati=LAT[indlat].values
	Loni=LON[indlon].values
	
	dist_array = np.zeros(len(Lati)-1)
	
	p=0
	for jj,ii in zip(range(len(Lati)-1),range(len(Loni)-1)):
	    lat1 = Lati[jj]
	    lon1 = Loni[ii]
	    lat2 = Lati[jj+1]
	    lon2 = Loni[ii+1]
	    dist_array[p]=  SVBfunc.haversine(lat1, lon1, lat2, lon2)
	    p=p+1
	
	
	dist_rot = np.cumsum(dist_array)
	dist_rot = np.insert(dist_rot,0,0)
	
	if np.any(dist_rot>=100):
	    hunKm=np.where(dist_rot>=100)[0][0]
	else:
	    hunKm=len(dist_rot)
	
	iniX.append(indlon[:hunKm])
	iniY.append(indlat[:hunKm])
	dist.append(dist_rot[:hunKm])
	deppre=depth.values[indlat[:hunKm],indlon[:hunKm]]
	dep.append(deppre[np.argpartition(deppre,range(10))])




mdic = {"dist": dist, "d":dep, 'indexXlon':iniX,'indexYlat':iniY }
if len(hej2)>20:
	savemat("/home/athelandersson/CTW-analysis/Files/" + str(coast) + "/BT_PALL.mat", mdic)
else:
	savemat("/home/athelandersson/CTW-analysis/Files/" + str(coast) + "/BT_P.mat", mdic)


pathVEL='/home/athelandersson/NETCDFs/' + str(coast) + '/WVELAC.nc'
dsVEL= xr.open_dataset(pathVEL)

WVEL=dsVEL.ValfiltAll.values[tstart:]
distVEL=dsVEL.dist.values
if coast == 'original':
	TIMEVEL=dsVEL.time2.values
elif coast == 'smooth':
	TIMEVEL=dsVEL.time2.values[tstart:]/60

lat_acVEL=dsVEL.latAC.values
lon_acVEL=dsVEL.lonAC.values

ETA=ds.ValfiltAll.values[tstart:]

params = {'font.size': 16,
          'figure.figsize': (11, 12),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

fig = plt.figure()
gs = GridSpec(nrows=2, ncols=2)

ax = fig.add_subplot(gs[0, 0])
ax.set_facecolor('tan')
pc = ax.contourf(LON, LAT, np.ma.masked_array(depth, mask=mask), 50,
	         vmin=0, vmax=5000, cmap=cmocean.cm.deep) 

cn = ax.contour(LON, LAT, depth, colors=['0.2', '0.4', '0.6', '0.8'],
	        levels=[200,500, 1000, 2000])
ax.contour(LON, LAT, depthno[:, :], levels=[0], colors='brown', linestyles=':', linewidths=2.5)
divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(pc, cax=axdiv,orientation='horizontal',ticks=[0,1000,2000,3000,4000])
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Depth [m]')
ax.set_xlabel('Lon [°]')
ax.set_ylabel('Lat [°]')

ax.set_aspect(1)
ax.text(-0.1, 1.2, '(a)', fontweight='bold', color='k',transform=ax.transAxes)

colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999','tab:blue']

ax1 = fig.add_subplot(gs[0, 1])
ax1.set_xlabel('Distance from coast [km]')
ax1.set_ylabel('Depth [m]')

for i in range(len(dep)):
	if len(hej2)>20:
		ax.scatter(LON[iniX[i]].values,LAT[iniY[i]].values,linewidth=2)
		ax1.plot(dist[i],-dep[i],linewidth=2)
	else:
                ax.scatter(LON[iniX[i]].values,LAT[iniY[i]].values,color=colors[i],linewidth=2)
                ax1.plot(dist[i],-dep[i],color=colors[i],linewidth=2)


ax1.text(-0.1, 1.02, '(b)', fontweight='bold', color='k', 
        transform=ax1.transAxes)

ax = fig.add_subplot(gs[1, 0])

vmin=-5
vmax=5
cbarall=0
SVBfunc.plot_HOVMOLLER(ax,distVEL,TIMEVEL,WVEL*1e6,'','Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,1,cbarall,'(c)')

for i in range(len(hej2)):
	if len(hej2)>20:
		ax.axhline(y=distAC[hej2[i]],linewidth=2,alpha=0.7)
	else:
		ax.axhline(y=distAC[hej2[i]],color=colors[i],linewidth=2,alpha=0.7)

ax = fig.add_subplot(gs[1, 1])

vmin=-0.2
vmax=0.2
cbarall=0
SVBfunc.plot_HOVMOLLER(ax,distAC,TIMEVEL,ETA*1e3,'','SSH  [mm]',vmin,vmax,fig,lat_ac,lon_ac,1,cbarall,'(d)')

for i in range(len(hej2)):
	if len(hej2)>20:
		ax.axhline(y=distAC[hej2[i]],linewidth=2,alpha=0.7)
	else:
		ax.axhline(y=distAC[hej2[i]],color=colors[i],linewidth=2,alpha=0.7)

fig.tight_layout()

if len(hej2)>20:
	plt.savefig('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + '/indsperpALL.png')	
else:
	plt.savefig('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + '/indsperp.png')	
