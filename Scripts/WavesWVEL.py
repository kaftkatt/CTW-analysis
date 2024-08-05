from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import xarray as xr
import cmocean
import pylab as pl

from SVBfunc import haversine

coast = 'smooth'

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


coast = 'original'

i=0
varname='dynVars'
pathn='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/'+ str(varname)+'noSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
pathw='/home/athelandersson/NETCDFs/'+ str(coast) +'/'+ str(varname)+'withSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
        
dswo  = xr.open_dataset(pathw)
dsno = xr.open_dataset(pathn)

Wwo=dswo.WVEL[0,55,:,:].values
Wno=dsno.WVEL[0,55,:,:].values
timeo=dswo.time.values.astype(int)
TIMEvelo=timeo*1e-9


i2=2
varname='dynVars'
pathn2='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/'+ str(varname)+'noSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
pathw2='/home/athelandersson/NETCDFs/'+ str(coast) +'/'+ str(varname)+'withSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
        
dsw2o  = xr.open_dataset(pathw2)
dsn2o = xr.open_dataset(pathn2)

Ww2o=dsw2o.WVEL[0,55,:,:].values
Wn2o=dsn2o.WVEL[0,55,:,:].values
time2o=dsw2o.time.values.astype(int)
TIMEvel2o=time2o*1e-9

i2=6
varname='dynVars'
pathn3='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/'+ str(varname)+'noSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
pathw3='/home/athelandersson/NETCDFs/'+ str(coast) + '/'+ str(varname)+'withSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
        
dsw3o  = xr.open_dataset(pathw3)
dsn3o = xr.open_dataset(pathn3)

Ww3o=dsw3o.WVEL[0,55,:,:].values
Wn3o=dsn3o.WVEL[0,55,:,:].values
time3o=dsw3o.time.values.astype(int)
TIMEvel3o=time3o*1e-9

LON=dsw.XC - 360
LAT=dsw.YC

hFacC = dsn.hFacC.values

hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)

deptho=dswo.Depth.values
depth=dsw.Depth.values

params = {'font.size': 16,
          'figure.figsize': (10, 20),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

fig = plt.figure()
gs = GridSpec(nrows=3, ncols=2,hspace=0.01)

vmin=-0.000002*1e6
vmax=0.000002*1e6

xlab='Longitude [°]'
ylab='Latitude [°]'
ctitle='Vertical velocity [$10^{-6}$ ms$^{-1}$]'

ax = fig.add_subplot(gs[0,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(Ww-Wn, mask=mask[55,:,:])*1e6,cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth, colors=['0.2','0.6'], 
                levels=[0,500])


ax.set(ylabel=ylab)
ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)

ax.text(0.4,0.87, f'480 m depth \nDay {TIMEvel[0]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.2, '(b)', transform=ax.transAxes)

ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

ax = fig.add_subplot(gs[1,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(Ww2-Wn2, mask=mask[55,:,:])*1e6,cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth, colors=['0.2','0.6'], 
                levels=[0,500])

ax.set(ylabel=ylab)
ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)

ax.text(0.4,0.87, f'480 m depth \nDay {TIMEvel2[0]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.05, '(d)', transform=ax.transAxes)


ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

ax = fig.add_subplot(gs[2,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(Ww3-Wn3, mask=mask[55,:,:])*1e6,cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth, colors=['0.2','0.6'], 
                levels=[0,500])

#ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=True) 
ax.set(xlabel=xlab)
ax.set(ylabel=ylab)

ax.text(0.4,0.87, f'480 m depth \nDay {TIMEvel3[0]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.05, '(f)', transform=ax.transAxes)


ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)


vmin=-0.000002*1e6
vmax=0.000002*1e6


ax = fig.add_subplot(gs[0,1]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(Wwo-Wno, mask=mask[55,:,:])*1e6,cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,deptho, colors=['0.2','0.6'], 
                levels=[0,500])

#divider = make_axes_locatable(ax)
#axdiv = divider.new_vertical(size = '5%', pad = 0.5)
#fig.add_axes(axdiv)
#cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
#cbar_ax.ax.xaxis.set_label_position("top")
#cbar_ax.set_label('Vertical velocity [$10^{-6}$ ms$^{-1}$]')

ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False) 
ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)

ax.text(0.4,0.87, f'480 m depth \nDay {TIMEvelo[0]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.2, '(b)', transform=ax.transAxes)

ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

ax = fig.add_subplot(gs[1,1]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(Ww2o-Wn2o, mask=mask[55,:,:])*1e6,cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,deptho, colors=['0.2','0.6'], 
                levels=[0,500])

ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False) 
ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)

ax.text(0.4,0.87, f'480 m depth \nDay {TIMEvel2o[0]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.05, '(d)', transform=ax.transAxes)


ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

ax = fig.add_subplot(gs[2,1]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(Ww3o-Wn3o, mask=mask[55,:,:])*1e6,cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,deptho, colors=['0.2','0.6'], 
                levels=[0,500])

ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False) 
ax.set(xlabel=xlab)

ax.text(0.4,0.87, f'480 m depth \nDay {TIMEvel3o[0]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.05, '(f)', transform=ax.transAxes)


ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

cbar_ax = fig.add_axes([1, 0.15, 0.03, 0.7])
fig.colorbar(cax, cax=cbar_ax)
cbar_ax.set_ylabel(ctitle)

plt.savefig('/home/athelandersson/CTW-analysis/Figures/both/fig2WWEL.png')
