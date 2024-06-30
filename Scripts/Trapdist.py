import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cmocean
import pylab as pl

from SVBfunc import haversine

i=0
varname='DYNVARS'
pathn='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL/01_noSVB_febTS/'+ str(varname)+'noSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
pathw='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01_SVB_febTS/'+ str(varname)+'withSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
        
dsw  = xr.open_dataset(pathw)
dsn = xr.open_dataset(pathn)

Ww=dsw.WVEL[:,55,220,:].values
Wn=dsn.WVEL[:,55,220,:].values
time=dsw.time.values.astype(int)
TIMEvel=time*1e-9

dist_array = np.zeros(len(LON))

p=0
for ii in np.arange(len(LON),0,-1):
    lat1 = LAT[220]
    lon1 = LON[len(LON)-1]
    lat2 = LAT[220]
    lon2 = LON[ii-1]
    p=p+1
    dist_array[p-1]=  haversine(lat1, lon1, lat2, lon2)



params = {'font.size': 16,
          'figure.figsize': (10, 6),
         'font.family':'sans'}
pl.rcParams.update(params)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


plt.pcolormesh(TIME/(60*24),dist_array[350:400],etafiltall[:,220,350:400].T*1000,cmap=cmocean.cm.curl)
cbar= plt.colorbar()
cbar.set_label('SSH [mm]')
plt.xlabel('Time [days]')
plt.ylabel('Distance [km]')



