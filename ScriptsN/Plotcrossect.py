import xarray as xr
import numpy as np
import SVBfunc
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import cmocean
from scipy.io import loadmat

hej=[35,54,79,120,154,194,219]   #use one of these indices to create pressure plots in order: 30.84,31.20,32.68,32.98°N
corrind=[30.49,30.77,31.13,31.69,32.11,32.65,33.02]
varlist=['PHIHYD','UVEL','VVEL']


i=0
var=varlist[0]
filt='filt'

FILENAMEfilt='Locations/' + str(var) + str(corrind[i]) + str(filt) +'.nc'
ds=xr.open_dataset(FILENAMEfilt)


dirw='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01b_SVB_febTS/'
varname='PHIHYD'
pathw = dirw + str(varname) + 'withSVB' + str(2 + i) + '_' + str(3 + i) + '.nc'
dsw = xr.open_dataset(pathw)

matfile=loadmat('BT_P.mat')
x,dep,indXlon,indYlat=matfile['dist'],matfile['d'],matfile['indexXlon'],matfile['indexYlat']

VALmit=ds.VAL
grid_X=ds.x.values
grid_Z=ds.z.values

hFacC = dsw.hFacC
hfac = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfac)
maskin=mask[:,indYlat[i],indXlon[i]]

t=200
vmin=-np.nanmax(abs(VALmit[t]))
vmax=np.nanmax(abs(VALmit[t]))

levels=np.linspace(vmin,vmax,15)

maxZ=np.max(np.where(maskin==False)[0])

varbrink=VALmit[t]

xlab='Cross-shore distance [km]'
ylab='Depth [m]'

fig,ax=plt.subplots()
ax.set_facecolor('tan')
cax=ax.contourf(grid_X,grid_Z[:maxZ],np.ma.masked_array(VALmit[t,:maxZ,:],mask=maskin[:maxZ,:]) ,cmap=cmocean.cm.delta,levels=levels)
ax.contour(grid_X,grid_Z[:maxZ],np.ma.masked_array(VALmit[t,:maxZ,:],mask=maskin[:maxZ,:]) , levels=[0], linewidths=2, 
        linestyles='-', colors='k', zorder=2)
ax.set(xlabel=xlab,ylabel=ylab)

plt.colorbar(cax, extend='max',label=var)
