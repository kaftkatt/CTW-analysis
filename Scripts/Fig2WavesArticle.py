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

Ww=dsw.WVEL[0,51,:,:].values
Wn=dsn.WVEL[0,51,:,:].values
time=dsw.time.values.astype(int)
TIMEvel=time*1e-9


i2=2
varname='dynVars'
pathn2='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/'+ str(varname)+'noSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
pathw2='/home/athelandersson/NETCDFs/'+ str(coast) +'/'+ str(varname)+'withSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
        
dsw2  = xr.open_dataset(pathw2)
dsn2 = xr.open_dataset(pathn2)

Ww2=dsw2.WVEL[0,51,:,:].values
Wn2=dsn2.WVEL[0,51,:,:].values
time2=dsw2.time.values.astype(int)
TIMEvel2=time2*1e-9


i=0
varname='eta'
pathn='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/'+ str(varname)+'noSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
pathw='/home/athelandersson/NETCDFs/'+ str(coast) +'/'+ str(varname)+'withSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
        
dsw  = xr.open_dataset(pathw)
dsn = xr.open_dataset(pathn)

ETAw1=dsw.ETAN[0,:,:].values
ETAn1=dsn.ETAN[0,:,:].values
time=dsw.time.values.astype(int)
TIME=time*1e-9


i2=2
varname='eta'
pathn2='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/'+ str(varname)+'noSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
pathw2='/home/athelandersson/NETCDFs/'+ str(coast) +'/'+ str(varname)+'withSVB'+ str(2+i2)+'_'+ str(3+i2) +'.nc'
        
dsw2  = xr.open_dataset(pathw2)
dsn2 = xr.open_dataset(pathn2)

ETAw2=dsw2.ETAN[0,:,:].values
ETAn2=dsn2.ETAN[0,:,:].values
time2=dsw2.time.values.astype(int)
TIME2=time2*1e-9

LON= dsn.XC.values-360
LAT= dsn.YC.values

hFacC = dsn.hFacC.values

hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)

depth=dsw.Depth.values

params = {'font.size': 16,
          'figure.figsize': (20, 20),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

fig = plt.figure()
gs = GridSpec(nrows=2, ncols=2,hspace=0.1)

if coast == 'smooth': 
	vminSSH=-0.01
	vmaxSSH=0.01
	vminvel=-0.000002*1e6
	vmaxvel=0.000002*1e6
elif coast == 'original': 
        vminSSH=-0.0001
        vmaxSSH=0.0001
        vminvel=-0.000002*1e6
        vmaxvel=0.000002*1e6

xlab='Longitude [°]'
ylab='Latitude [°]'

if coast== 'original':
	timelab=TIME[0]/(60*24)
else:
	timelab=TIME2[0]/(60*60*24)
ax = fig.add_subplot(gs[0,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array((ETAw1-ETAn1)*1000, mask=mask[0,:,:]),cmap=cmocean.cm.curl,vmin=vminSSH,vmax=vmaxSSH)
ax.contour(LON,LAT,depth,  colors=['0.2','0.4','0.6','0.8'], 
                levels=[0,500,700,1000])

ax.tick_params(axis='x',which='both', top=False, bottom=True, labelbottom=False)

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('SSH [mm]')

ax.text(-0.1,1.2, '(a)', transform=ax.transAxes)
ax.text(0.4,0.87, f'Surface \nDay {timelab}', transform=ax.transAxes,horizontalalignment='left')

ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)
ax.set(ylabel=ylab)

ax = fig.add_subplot(gs[1,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array((ETAw2-ETAn2)*1000, mask=mask[0,:,:]),cmap=cmocean.cm.curl,vmin=vminSSH,vmax=vmaxSSH)
ax.contour(LON,LAT,depth,  colors=['0.2','0.4','0.6','0.8'], 
                levels=[0,500,700,1000])

ax.text(-0.1,1.2, '(c)', transform=ax.transAxes)
ax.text(0.4,0.87, f'Surface \nDay {timelab}', transform=ax.transAxes,horizontalalignment='left')

ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)
ax.set(ylabel=ylab,xlabel=xlab)


ax = fig.add_subplot(gs[0,1]) 
ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(Ww-Wn, mask=mask[51,:,:])*1e6,cmap=cmocean.cm.balance,vmin=vminvel,vmax=vmaxvel)
ax.contour(LON,LAT,depth, colors=['0.2','0.4','0.6','0.8'], 
                levels=[0,500,700,1000])

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Vertical velocity [$10^{-6}$ ms$^{-1}$]')

ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False)

ax.tick_params(axis='x',which='both', top=False, bottom=True, labelbottom=False)

ax.text(0.4,0.87, f'{dsw.Z.values[51]:.1f} m depth \nDay {timelab}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.2, '(b)', transform=ax.transAxes)


ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

ax = fig.add_subplot(gs[1,1]) 
ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(Ww2-Wn2, mask=mask[51,:,:])*1e6,cmap=cmocean.cm.balance,vmin=vminvel,vmax=vmaxvel)
ax.contour(LON,LAT,depth, colors=['0.2','0.4','0.6','0.8'], 
                levels=[0,500,700,1000])

ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False)

ax.set(xlabel=xlab)

ax.text(0.4,0.87, f'{dsw.Z.values[51]:.1f} m depth \nDay {timelab}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.15, '(d)', transform=ax.transAxes)


ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

plt.savefig('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + '/fig2Waves.png',bbox_inches='tight')
