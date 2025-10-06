from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import xarray as xr
import cmocean
import pylab as pl

from SVBfunc import loadNetCDFs
from SVBfuncPlotting import haversine

coast = 'smooth'
tstart=2
i=0
dirn = '/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
dirw = '/home/athelandersson/NETCDFs/' + str(coast) + '/'

dsw, dsn = loadNetCDFs(dirw, dirn, 'dynVars',tstart)

var=np.zeros((len(dsw[0].WVEL[0,55,:,0]),len(dsw[0].WVEL[0,55,0,:])))
indsIn=[]
indsIn.append(0)
for l in range(len(dsw[0].WVEL[0,55,:,0])):
	vel=np.zeros((576,len(dsw[0].WVEL[0,55,0,:])))
	for k in range(0,len(dsw)):
		Ww=dsw[k].WVEL[:,55,l,:].values
		Wn=dsn[k].WVEL[:,55,l,:].values
		indsIn.append(len(dsw[k].WVEL[:,55,0,0]))
		vel[(72*k):(72*(k+1)),:]=Ww-Wn
	var[l,:]=np.var(vel,axis=0)
	print(l)


pathETA='/home/athelandersson/NETCDFs/' + str(coast) + '/ETA.nc'
dsETA= xr.open_dataset(pathETA)

etafiltall=dsETA.VALfilt.values
LON=dsETA.x.values
LAT=dsETA.y.values
TIME=dsETA.time.values/60
hFacC = dsn[0].hFacC.values

hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)

depth=dsw[0].Depth.values

coast = 'original'

i=0
dirn = '/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
dirw = '/home/athelandersson/NETCDFs/' + str(coast) + '/'

dsw, dsn = loadNetCDFs(dirw, dirn, 'dynVars',tstart)

varo=np.zeros((len(dsw[0].WVEL[0,55,:,0]),len(dsw[0].WVEL[0,55,0,:])))
indsIn=[]
indsIn.append(0)
for l in range(len(dsw[0].WVEL[0,55,:,0])):
        vel=np.zeros((576,len(dsw[0].WVEL[0,55,0,:])))
        for k in range(0,len(dsw)):
                Ww=dsw[k].WVEL[:,55,l,:].values
                Wn=dsn[k].WVEL[:,55,l,:].values
                vel[72*k:(72*(k+1)),:]=Ww-Wn
        varo[l,:]=np.var(vel,axis=0)
        print(l)

pathETA='/home/athelandersson/NETCDFs/' + str(coast) + '/ETA.nc'
dsETA= xr.open_dataset(pathETA)

etafiltallo=dsETA.VALfilt.values
LONo=dsETA.x.values
LATo=dsETA.y.values
TIMEo=dsETA.time.values
hFacCo = dsn[0].hFacC.values

hfao = np.ma.masked_values(hFacCo, 0)
masko = np.ma.getmask(hfao)

deptho=dsw[0].Depth.values
coast='smooth'
dir = '/home/athelandersson/NETCDFs/' + str(coast) + '/'
pathVEL=str(dir) + 'WVELAC.nc'
dsVEL= xr.open_dataset(pathVEL)

lat_acVEL=dsVEL.latAC.values
lon_acVEL=dsVEL.lonAC.values

maxpow=[20,75,204,248]
maxpowOrg=[60,130,183,229]

params = {'font.size': 16,
          'figure.figsize': (10, 10),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

fig = plt.figure()
gs = GridSpec(nrows=2, ncols=2)

vmin= 0 #np.max(np.ma.masked_array(np.var(etafiltall[:,:,:],axis=0), mask=mask[0,:,:]))/100
vmax=np.max(np.ma.masked_array(np.var(etafiltall[:,:,:],axis=0), mask=mask[0,:,:]))/5e1


xlab='Longitude [°]'
ylab='Latitude [°]'

ax = fig.add_subplot(gs[0,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(np.var(etafiltall[:,:,:],axis=0), mask=mask[0,:,:]),cmap=cmocean.cm.amp,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth,  colors=['0.2','0.6'], 
                levels=[0,500])

ind_lon = [ -115.939167, -116.605833, -117.1625]
ind_lat = [30.556389, 31.857778, 32.715]

for kk, ll, lab in zip(ind_lon, ind_lat,
			['San \nQuintín', 'Ensenada', 'San Diego']):
	ax.plot(kk, ll, 'o',markersize=4, color='tab:green', markeredgecolor='k',zorder=1)
	if lab =='San \nQuintín':
		ax.text(kk-0.23 , ll + 0.2, lab, fontsize=14)
	else:
		ax.text(kk , ll + 0.18, lab, fontsize=14)
ax.scatter(lon_acVEL[maxpowOrg],lat_acVEL[maxpowOrg],color='red',s=4,zorder=4,label='Original')
ax.scatter(lon_acVEL[maxpow],lat_acVEL[maxpow],color='blue',s=4,zorder=4,label='Smooth')
ax.legend()

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Variance in SSH')

ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)
ax.set( ylabel='Latitude [°]')

ax.text(-0.1,1.2, '(a)', transform=ax.transAxes)
ax.text(0.4,0.93, f'Smoothened \nSurface', transform=ax.transAxes,horizontalalignment='left')
ax.set_xlim((-118.5,-115.7))
ax.set_ylim((29,34))

vmin= 0 #np.max(np.ma.masked_array(np.var(etafiltall[:,:,:],axis=0), mask=mask[0,:,:]))/100
vmax=np.max(np.ma.masked_array(np.var(etafiltall[:,:,:],axis=0), mask=mask[0,:,:]))/100

ax = fig.add_subplot(gs[1,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LONo,LAT,np.ma.masked_array(np.var(etafiltallo[:,:,:],axis=0), mask=masko[0,:,:]),cmap=cmocean.cm.amp,vmin=vmin,vmax=vmax)
ax.contour(LONo,LATo,deptho,  colors=['0.2','0.6'], 
                levels=[0,500])
ind_lon = [ -115.939167, -116.605833, -117.1625]
ind_lat = [30.556389, 31.857778, 32.715]

for kk, ll, lab in zip(ind_lon, ind_lat,
			['San \nQuintín', 'Ensenada', 'San Diego']):
	ax.plot(kk, ll, 'o',markersize=4, color='tab:green', markeredgecolor='k',zorder=1)
	if lab =='San \nQuintín':
		ax.text(kk-0.23 , ll + 0.2, lab, fontsize=14)
	else:
		ax.text(kk , ll + 0.18, lab, fontsize=14)
ax.scatter(lon_acVEL[maxpowOrg],lat_acVEL[maxpowOrg],color='red',s=4,zorder=4,label='Original')
ax.scatter(lon_acVEL[maxpow],lat_acVEL[maxpow],color='blue',s=4,zorder=4,label='Smooth')

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Variance in SSH')

ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=True)
ax.tick_params(axis='y',which='both', left=False, labelleft=True)
ax.set(xlabel='Longitude [  ]',ylabel='Latitude [°]')

ax.text(-0.1,1.2, '(c)', transform=ax.transAxes)
ax.text(0.4,0.93, f'Original \nSurface', transform=ax.transAxes,horizontalalignment='left')
ax.set_xlim((-118.5,-115.7))
ax.set_ylim((29,34))


vmin=0
vmax=6e-6


ax = fig.add_subplot(gs[0,1]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(var, mask=mask[55,:,:])*1e6,cmap=cmocean.cm.amp,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth, colors=['0.2','0.6'], 
                levels=[0,500])
ind_lon = [ -115.939167, -116.605833, -117.1625]
ind_lat = [30.556389, 31.857778, 32.715]

for kk, ll, lab in zip(ind_lon, ind_lat,
			['San \nQuint  n', 'Ensenada', 'San Diego']):
	ax.plot(kk, ll, 'o',markersize=4, color='tab:green', markeredgecolor='k',zorder=1)
	if lab =='San \nQuint  n':
		ax.text(kk-0.23 , ll + 0.2, lab, fontsize=14)
	else:
		ax.text(kk , ll + 0.18, lab, fontsize=14)
ax.scatter(lon_acVEL[maxpowOrg],lat_acVEL[maxpowOrg],color='red',s=4,zorder=4,label='Original')
ax.scatter(lon_acVEL[maxpow],lat_acVEL[maxpow],color='blue',s=4,zorder=4,label='Smooth')

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Variance in Vertical Velocity')

ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False) 
ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)

ax.text(0.4,0.93, f'Smoothened \n480 m depth', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.2, '(b)', transform=ax.transAxes)
ax.set_xlim((-118.5,-115.7))
ax.set_ylim((29,34))

ax = fig.add_subplot(gs[1,1]) 
vmin=0
vmax=2e-6
ax.set_facecolor('tan')
cax = ax.pcolormesh(LONo,LATo,np.ma.masked_array(varo, mask=masko[55,:,:])*1e6,cmap=cmocean.cm.amp,vmin=vmin,vmax=vmax)
ax.contour(LONo,LATo,deptho, colors=['0.2','0.6'], 
                levels=[0,500])
ind_lon = [ -115.939167, -116.605833, -117.1625]
ind_lat = [30.556389, 31.857778, 32.715]

for kk, ll, lab in zip(ind_lon, ind_lat,
		['San \nQuint  n', 'Ensenada', 'San Diego']):
	ax.plot(kk, ll, 'o',markersize=4, color='tab:green', markeredgecolor='k',zorder=1)
	if lab =='San \nQuint  n':
		ax.text(kk-0.23 , ll + 0.2, lab, fontsize=14)
	else:
		ax.text(kk , ll + 0.18, lab, fontsize=14)
ax.scatter(lon_acVEL[maxpowOrg],lat_acVEL[maxpowOrg],color='red',s=4,zorder=4,label='Original')
ax.scatter(lon_acVEL[maxpow],lat_acVEL[maxpow],color='blue',s=4,zorder=4,label='Smooth')

ax.set(xlabel='Longitude [°]')

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Variance in Vertical Velocity')

ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False) 
ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=True)

ax.text(0.4,0.93, f'Original \n480 m depth', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.2, '(d)', transform=ax.transAxes)

ax.set_xlim((-118.5,-115.7))
ax.set_ylim((29,34))

plt.tight_layout()
plt.savefig('/home/athelandersson/CTW-analysis/Figures/Variance.png')
