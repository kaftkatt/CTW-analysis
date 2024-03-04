import SVBfunc
import numpy as np
import xarray as xr
from matplotlib.gridspec import GridSpec
import pylab as pl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#dirw='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01_SVB_febTS/'
#dirn='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL/01_noSVB_febTS/'

#dsw,dsn=loadNetCDFs(dirw,dirn,'DYNVARS')


pathw='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL/01_noSVB_febTS/DYNVARSnoSVB'+ str(2)+'_'+ str(3) +'.nc'
dsw  = xr.open_dataset(pathw)

TIME=2

moment=TIME*60*24
t=(moment-2880)/20
VAL=get_snapshot_at_level(t,dep,dsw,dsn)


pathETA='/media/amelia/Trillian/SVB/ETA.nc'
dsETA= xr.open_dataset(pathETA)

etafiltall=dsETA.ETAfiltall.values
LON=dsETA.x.values
LAT=dsETA.y.values
TIME=dsETA.time.values

hFacC = dsw.hFacC.values

hfa = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfa)

depth=dsw.Depth.values

params = {'font.size': 16 ,
          'figure.figsize': (20, 30),
         'font.family':'sans'}
pl.rcParams.update(params)
#plt.rcParams['figure.dpi'] = 300
#plt.rcParams['savefig.dpi'] = 300

fig = plt.figure()
gs = GridSpec(nrows=3, ncols=2)#,hspace=0.01)

vmin=-0.009*10
vmax=0.009*10

ax = fig.add_subplot(gs[0,0]) 

label='Surface'
nr='(a)'
dep=0
TIME=2
moment=TIME*60*24
t=int((moment-2880)/20)

VAL=etafiltall[t,:,:]*1000
cax=SVBfunc.plotsnapshot(ax,VAL,dep,LON,LAT,vmin,vmax,depth,label,nr,TIME,mask)

ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.1)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('SSH [mm]')

ax = fig.add_subplot(gs[1,0]) 



label='Surface'
nr='(c)'
dep=0
TIME=4
moment=TIME*60*24
t=int((moment-2880)/20)

VAL=etafiltall[t,:,:]*1000
cax=SVBfunc.plotsnapshot(ax,VAL,dep,LON,LAT,vmin,vmax,depth,label,nr,TIME,mask)

ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)


ax = fig.add_subplot(gs[2,0]) 

label='Surface'
nr='(e)'
dep=0
TIME=8
moment=TIME*60*24
t=int((moment-2880)/20)

VAL=etafiltall[t,:,:]*1000
cax=SVBfunc.plotsnapshot(ax,VAL,dep,LON,LAT,vmin,vmax,depth,label,nr,TIME,mask)
