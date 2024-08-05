from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import xarray as xr
import cmocean
import pylab as pl

from SVBfunc import haversine

coast = 'smooth'

pathETA='/home/athelandersson/NETCDFs/'+ str(coast) + '/ETA.nc'
pathmask='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/etanoSVB2_3.nc'
dsn=xr.open_dataset(pathmask)
dsETA= xr.open_dataset(pathETA)

etafiltall=dsETA.VALfilt.values
LON=dsETA.x.values
LAT=dsETA.y.values
TIME=dsETA.time.values

hFacC = dsn.hFacC.values

hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)

depth=dsn.Depth.values

coast = 'original'

pathETA='/home/athelandersson/NETCDFs/'+ str(coast) + '/ETA.nc'
dsETA= xr.open_dataset(pathETA)

pathmask='/home/athelandersson/NETCDFs/'+ str(coast) + '_NO/etanoSVB2_3.nc'
dsn=xr.open_dataset(pathmask)

etafiltallo=dsETA.VALfilt.values
TIMEo=dsETA.time.values
hFacCo = dsn.hFacC.values

hfao = np.ma.masked_values(hFacCo, 0)
masko = np.ma.getmask(hfao)

deptho=dsn.Depth.values

params = {'font.size': 16,
          'figure.figsize': (10, 20),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

fig = plt.figure()
gs = GridSpec(nrows=3, ncols=2,hspace=0.01)

vmin=-0.009*10
vmax=0.009*10
ind=np.where(TIME==2*60*60*24)[0][0]
ind2=np.where(TIME==4*60*60*24)[0][0]
ind3=np.where(TIME==8*60*60*24)[0][0]

indo=np.where(TIMEo==2*60*24)[0][0]
ind2o=np.where(TIMEo==4*60*24)[0][0]
ind3o=np.where(TIMEo==8*60*24)[0][0]


xlab='Longitude [°]'
ylab='Latitude [°]'

ax = fig.add_subplot(gs[0,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(etafiltall[ind,:,:]*1000, mask=mask[0,:,:]),cmap=cmocean.cm.curl,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth,  colors=['0.2','0.6'], 
                levels=[0,500])



ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)
ax.set( ylabel=ylab)

ax.text(-0.1,1.2, '(a)', transform=ax.transAxes)
ax.text(0.4,0.87, f'Surface \nDay {TIME[ind]/(60*60*24)}', transform=ax.transAxes,horizontalalignment='left')

ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

ax = fig.add_subplot(gs[1,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(etafiltall[ind2,:,:]*1000, mask=mask[0,:,:]),cmap=cmocean.cm.curl,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth,  colors=['0.2','0.6'], 
                levels=[0,500])

ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)
ax.set(ylabel=ylab)

ax.text(-0.1,1.05, '(c)', transform=ax.transAxes)
ax.text(0.4,0.87, f'Surface \nDay {TIME[ind2]/(60*60*24)}', transform=ax.transAxes,horizontalalignment='left')

ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

ax = fig.add_subplot(gs[2,0]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(etafiltall[ind3,:,:]*1000, mask=mask[0,:,:]),cmap=cmocean.cm.curl,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,depth,  colors=['0.2','0.6'], 
                levels=[0,500])


ax.set(xlabel=xlab, ylabel=ylab)

ax.text(-0.1,1.05, '(e)', transform=ax.transAxes)
ax.text(0.4,0.87, f'Surface \nDay {TIME[ind3]/(60*60*24)}', transform=ax.transAxes,horizontalalignment='left')

ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)


ax = fig.add_subplot(gs[0,1]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(etafiltallo[indo,:,:]*1000, mask=masko[0,:,:]),cmap=cmocean.cm.curl,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,deptho,  colors=['0.2','0.6'], 
                levels=[0,500])

ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False) 
ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)

ax.text(-0.1,1.2, '(b)', transform=ax.transAxes)
ax.text(0.4,0.87, f'Original \nDay {TIMEo[indo]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')


ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

ax = fig.add_subplot(gs[1,1]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(etafiltallo[ind2o,:,:]*1000, mask=masko[0,:,:]),cmap=cmocean.cm.curl,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,deptho,  colors=['0.2','0.6'], 
                levels=[0,500])

ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False) 
ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)

ax.text(0.4,0.87, f'Original \nDay {TIMEo[ind2o]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.05, '(d)', transform=ax.transAxes)


ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

ax = fig.add_subplot(gs[2,1]) 

ax.set_facecolor('tan')
cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(etafiltallo[ind3o,:,:]*1000, mask=masko[0,:,:]),cmap=cmocean.cm.curl,vmin=vmin,vmax=vmax)
ax.contour(LON,LAT,deptho,  colors=['0.2','0.6'], 
                levels=[0,500])

ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False) 
ax.set(xlabel=xlab)

ax.text(0.4,0.87, f'Original \nDay {TIMEo[ind3o]/(60*24)}', transform=ax.transAxes,horizontalalignment='left')
ax.text(-0.1,1.05, '(f)', transform=ax.transAxes)


ax.set_xlim(-122,-114) 
ax.set_ylim(27,35.3)
ax.set_aspect(1)

cb_ax = fig.add_axes([.91,.124,.04,.754])
fig.colorbar(cax,orientation='vertical',cax=cb_ax)
cb_ax.set_label('SSH [mm]')

plt.savefig('/home/athelandersson/CTW-analysis/Figures/both/fig2SSH.png')
