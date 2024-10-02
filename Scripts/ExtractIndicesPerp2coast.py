import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import pylab as pl
import scipy.signal as sig
import scipy.io as sio
import scipy.interpolate as sciint
from scipy.interpolate import RegularGridInterpolator
from xmitgcm import open_mdsdataset
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
from os.path import exists

from scipy.io import savemat, loadmat
import SVBfunc

from math import radians, cos, sin, asin, sqrt, atan, degrees, log

plotALL=0
coast='original'
nr=1 #The larger the number the larger the distance over which the angle is calculated
#Thus the less accurate the angle. 

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

lon_ac=ds.lonAC.values-1
lat_ac=ds.latAC.values-1
distAC=ds.dist.values

#hej2=[35,54,79,120,154,194,235]
hej2=np.arange(nr,len(lon_ac)-nr,1)

LAT = dsw[0].YC
LON = dsw[0].XC - 360

Z = dsw[0].Z
hFacC = dsw[0].hFacC

hfa = np.ma.masked_values(hFacC[0, :, :], 0)
mask = np.ma.getmask(hfa)

depth = dsw[0].Depth
depthno = dsn[0].Depth

interp = RegularGridInterpolator((LAT.values,LON.values), depth.values)

shift=30
lonWeight30=SVBfunc.weightedmovingaverage(LON[lon_ac], shift)
latWeight30=SVBfunc.weightedmovingaverage(LAT[lat_ac], shift)

hej2ind=np.array((35,54,79,120,154,194,240))
loc=[]
for i in hej2ind: 
    loc.append(np.where(latWeight30>=LAT[lat_ac[i]].values)[0][0])


lonNew=[]
latNew=[]
dist=[]
dep=[]
degree=[]
finNR=[]

lonNewSHORT=[]
latNewSHORT=[]
distSHORT=[]
depSHORT=[]
degreeSHORT=[]

for ind in range(1,len(lonWeight30)-1,1):
    lon1=lonWeight30[ind-nr]
    lat1=latWeight30[ind-nr]
    lon2=lonWeight30[ind-nr]
    lat2=latWeight30[ind+nr]
    
    a=SVBfunc.haversine(lon1, lat1, lon2, lat2) 
    
    lon3=lonWeight30[ind+nr]
    lat3=latWeight30[ind+nr]
    
    b=SVBfunc.haversine(lon2, lat2, lon3, lat3)
    
    if lon1 == lon3: 
        deg=0
    elif lat1==lat3:
        deg=radians(270)
    else:
        deg=-atan(a/b)
    
    R1=LON*cos(deg)-LAT*sin(deg)
    R2=LON*sin(deg)+LAT*cos(deg)
    indexX,indexY=np.where(np.logical_and(dsw[0].XC==LON[lon_ac[ind+shift-1]]+360,dsw[0].YC==LAT[lat_ac[ind+shift-1]]))
   
    Rcoast = R1[indexX,indexY]
    
    LONIN=np.arange(np.min(R1),R1[indexX,indexY],dsw[0].XC[1].values-dsw[0].XC[0].values)
    LATIN=np.ones(len(LONIN))*R2[indexX,indexY].values[0,0]
    
    R1Back=np.flip(LONIN*cos(-deg)-LATIN*sin(-deg))
    R2Back=np.flip(LONIN*sin(-deg)+LATIN*cos(-deg))
    
    
    dist_array = np.zeros(len(R1Back)-1)
    
    p=0
    for jj in range(len(R1Back)-1):
        lat1 = R2Back[jj]
        lon1 = R1Back[jj]
        lat2 = R2Back[jj+1]
        lon2 = R1Back[jj+1]
        dist_array[p]=  SVBfunc.haversine(lat1, lon1, lat2, lon2)
        p=p+1
    
    dist_rot = np.cumsum(dist_array)
    dist_rot = np.insert(dist_rot,0,0)
    
    if np.any(dist_rot>=100):
        hunKm=np.where(dist_rot>=100)[0][0]
    else:
        hunKm=len(dist_rot)
    if hunKm>=len(R1Back):
        hunKm=len(R1Back)
    elif hunKm>=len(R1Back):
        hunKm=len(R1Back)
    
    if max(R2Back[:hunKm])>max(LAT):
        hunKm=np.where(R2Back<max(LAT).values)[0][-1]
    
    lonNew.append(R1Back[:hunKm])
    latNew.append(R2Back[:hunKm])
    dist.append(dist_rot[:hunKm])
    
    deppre=interp((R2Back[:hunKm],R1Back[:hunKm]))
    dep.append(deppre)
    degree.append(deg)

    if np.any(lat_ac[hej2ind]==lat_ac[ind+shift-1]):
        print(LAT[lat_ac[ind+shift-1]])
        lonNewSHORT.append(R1Back[:hunKm])
        latNewSHORT.append(R2Back[:hunKm])
        distSHORT.append(dist_rot[:hunKm])
        depSHORT.append(deppre)
        degreeSHORT.append(deg)


mdicALL = {"dist": dist, "d":dep, 'lon':lonNew,'lat':latNew,'degree':degree }

savemat('/home/athelandersson/CTW-analysis/Files/' + str(coast) + "/BT_PALL_MovAv.mat", mdicALL)

mdichej2 = {"dist": distSHORT, "d":depSHORT, 'lon':lonNewSHORT,'lat':latNewSHORT,'degree':degreeSHORT }

savemat('/home/athelandersson/CTW-analysis/Files/' + str(coast) + "/BT_P_MovAv.mat", mdichej2)

if plotALL==1:
	filenam="/BT_PALL_MovAv.mat"
else:
	filenam="/BT_P_MovAv.mat"

matfile=loadmat( '/home/athelandersson/CTW-analysis/Files/' + str(coast) + str(filenam))
x,dep,lon,lat,deg=matfile['dist'][0],matfile['d'][0],matfile['lon'][0],matfile['lat'][0],matfile['degree'][0]


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

colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999','tab:blue']

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

ax.scatter(LON[lon_ac[hej2ind]],LAT[lat_ac[hej2ind]])
ax.scatter(lonWeight30[loc],latWeight30[loc])

ax.set_aspect(1)
ax.text(-0.1, 1.2, '(a)', fontweight='bold', color='k',transform=ax.transAxes)


ax1 = fig.add_subplot(gs[0, 1])
ax1.set_xlabel('Distance from coast [km]')
ax1.set_ylabel('Depth [m]')

for i in range(len(dep)):
	if len(hej2ind)>20:
		ax.scatter(lon[i][0],lat[i][0],linewidth=2)
		ax1.plot(x[i][0],-dep[i][0],linewidth=2)
	else:
                ax.scatter(lon[i][0],lat[i][0],color=colors[i],linewidth=2)
                ax1.plot(x[i][0],-dep[i][0],color=colors[i],linewidth=2)


ax1.text(-0.1, 1.02, '(b)', fontweight='bold', color='k', 
        transform=ax1.transAxes)

ax = fig.add_subplot(gs[1, 0])

vmin=-5
vmax=5
cbarall=0
SVBfunc.plot_HOVMOLLER(ax,distVEL,TIMEVEL,WVEL*1e6,'','Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,1,cbarall,'(c)')

for i in range(len(hej2ind)):
	if len(hej2ind)>20:
		ax.axhline(y=distAC[hej2ind[i]],linewidth=2,alpha=0.7)
	else:
		ax.axhline(y=distAC[hej2ind[i]],color=colors[i],linewidth=2,alpha=0.7)

ax = fig.add_subplot(gs[1, 1])

vmin=-0.2
vmax=0.2
cbarall=0
SVBfunc.plot_HOVMOLLER(ax,distAC,TIMEVEL,ETA*1e3,'','SSH  [mm]',vmin,vmax,fig,lat_ac,lon_ac,1,cbarall,'(d)')

for i in range(len(hej2ind)):
	if len(hej2ind)>20:
		ax.axhline(y=distAC[hej2ind[i]],linewidth=2,alpha=0.7)
	else:
		ax.axhline(y=distAC[hej2ind[i]],color=colors[i],linewidth=2,alpha=0.7)


if len(dep)>20:
	plt.savefig('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + '/indsperpALL_res' + str(nr) + '.png')	
else:
	plt.savefig('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + '/indsperp_res' + str(nr) + '.png')	
