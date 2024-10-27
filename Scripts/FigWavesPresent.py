from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import xarray as xr
import cmocean
import pylab as pl

from SVBfunc import haversine

coast = 'original'

i=0
varname='dynVars'
pathn='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/'+ str(varname)+'noSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
pathw='/home/athelandersson/NETCDFs/'+ str(coast) +'/'+ str(varname)+'withSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
        
dsw  = xr.open_dataset(pathw)
dsn = xr.open_dataset(pathn)

Ww=dsw.WVEL[0,55,:,:].values
Wn=dsn.WVEL[0,55,:,:].values
time=dsw.time.values.astype(int)
TIMEvel=time*1e-9


i2=2
varname='dynVars'
pathn2='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/'+ str(varname)+'noSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
pathw2='/home/athelandersson/NETCDFs/'+ str(coast) +'/'+ str(varname)+'withSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
        
dsw2  = xr.open_dataset(pathw2)
dsn2 = xr.open_dataset(pathn2)

Ww2=dsw2.WVEL[0,55,:,:].values
Wn2=dsn2.WVEL[0,55,:,:].values
time2=dsw2.time.values.astype(int)
TIMEvel2=time2*1e-9

i2=6
varname='dynVars'
pathn3='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/'+ str(varname)+'noSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
pathw3='/home/athelandersson/NETCDFs/'+ str(coast) + '/'+ str(varname)+'withSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
        
dsw3  = xr.open_dataset(pathw3)
dsn3 = xr.open_dataset(pathn3)

Ww3=dsw3.WVEL[0,55,:,:].values
Wn3=dsn3.WVEL[0,55,:,:].values
time3=dsw3.time.values.astype(int)
TIMEvel3=time3*1e-9


pathETA='/home/athelandersson/NETCDFs/' + str(coast) + '/ETA.nc'
dsETA= xr.open_dataset(pathETA)

etafiltall=dsETA.VALfilt.values
LON=dsETA.x.values
LAT=dsETA.y.values
TIME=dsETA.time.values
hFacC = dsn.hFacC.values

hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)

depth=dsw.Depth.values

params = {'font.size': 16,
          'figure.figsize': (10, 20),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

fig = plt.figure()
gs = GridSpec(nrows=1, ncols=2,hspace=0.01)

vmin=-0.009*10
vmax=0.009*10
ind=np.where(TIME==TIMEvel[0])[0][0]
ind2=np.where(TIME==TIMEvel2[0])[0][0]
ind3=np.where(TIME==TIMEvel3[0])[0][0]


xlab='Longitude [°]'
ylab='Latitude [°]'

ax = fig.add_subplot(gs[0,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(etafiltall[ind2,:,:]*1000, mask=mask[0,:,:]),cmap=cmocean.cm.curl,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth,  colors=['0.2','0.6'], 
                levels=[0,500])

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('SSH [mm]')

ax.text(-0.1,1.2, '(a)', transform=ax.transAxes)
ax.text(0.4,0.87, f'Surface \nDay {TIME[ind2]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')

ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)
ax.set(xlabel=xlab, ylabel=ylab)

vmin=-0.000002*1e6
vmax=0.000002*1e6


ax = fig.add_subplot(gs[0,1]) 
ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(Ww2-Wn2, mask=mask[55,:,:])*1e6,cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth, colors=['0.2','0.6'], 
                levels=[0,500])

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Vertical velocity [$10^{-6}$ ms$^{-1}$]')

ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)

ax.set(ylabel=ylab)

ax.text(0.4,0.87, f'480 m depth \nDay {TIMEvel2[0]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.2, '(b)', transform=ax.transAxes)


ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)



plt.savefig('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + 'fig2Waves.png', transparent=True,bbox_inches='tight')
