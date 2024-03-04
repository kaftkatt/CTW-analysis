import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from scipy.io import loadmat

corrinds=[30.77,31.13,32.65,33.02]
VAR=['PHIHYD','VVEL','UVEL'] #,'WVEL']
FILT=['no','filt']

filt=FILT[1]
var=VAR[0]
ds=[]

for lat in corrinds:
	ds.append(xr.open_dataset(str(var)+'/LinReg' + str(lat) + str(filt)+ '.nc'))

ik=0
grid_X=ds[ik].gridX.values
grid_Z=ds[ik].gridZ.values
VALfit=ds[ik].valfit.values
VAL=ds[ik].valmit.values
TIME=ds[ik].time.values

matfile=loadmat('BT_PERP.mat')
indXlon,indYlat=matfile['indexXlon'][ik],matfile['indexYlat'][ik]
	
dirw='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01b_SVB_febTS/'
dirn='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL/01b_noSVB_febTS/'

dsw,dsn=SVBfunc.loadNetCDFs(dirw,dirn,'PHIHYD')
hFacC = dsw[0].hFacC

hfa = np.ma.masked_values(hFacC, 0)
maskin = np.ma.getmask(hfa)
mask=maskin[:,indYlat,indXlon]

vals=ds[ik].varbrink.values
modes=len(vals[:,1,1])


vminb=-np.nanmax(abs(vals))
vmaxb=np.nanmax(abs(vals))

levelsb=np.linspace(vminb,vmaxb,15)

day1=4
day2=4.6

lat=32.98

t=np.where(TIME>=day1*60*24)[0][0]

t2=np.where(TIME>=day2*60*24)[0][0]

vmin=-np.nanmax(abs(VAL[t]))*1e4
vmax=np.nanmax(abs(VAL[t]))*1e4

levels=np.linspace(vmin,vmax,15)

hfac = np.ma.masked_values(VALfit[t], 0)
mask = np.ma.getmask(hfac)

xlab='Cross-shore distance [km]'
ylab='Depth [m]'

fig = plt.figure()
gs = GridSpec(nrows=4, ncols=2, height_ratios=[1,1,1,1],hspace=0.35)
ax = fig.add_subplot(gs[0, 0])

ax.text(0.88, 1.34, 'At' + str(corrinds[ik])+'°N', fontweight='bold', fontsize=22,color='k', 
                transform=ax.transAxes)
                
cax1=SVBfunc.plotbrink(ax,grid_X,grid_Z,levels,xlab,ylab,-1,np.ma.masked_array(VAL[t]*1e4,mask=mask),modes,lat,day1)
ax.set_ylim([np.min(np.ma.masked_array(grid_Z,mask=mask)),0])
ax.text(0.96, 1.06, '(a)', fontweight='bold', color='k', 
            transform=ax.transAxes)
            
            
ax = fig.add_subplot(gs[0, 1],sharey=ax)
cax1=SVBfunc.plotbrink(ax,grid_X,grid_Z,levels,xlab,ylab,-2,np.ma.masked_array(VALfit[t]*1e4,mask=mask),modes,lat,day2)

cbar_ax1 = fig.add_axes([0.95, 0.74, 0.03, 0.13])
fig.colorbar(cax1, cax=cbar_ax1)
cbar_ax1.set_ylabel('Pressure [$10^{-4}$ m²s$^{-2}$]')
cbar_ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))

ax.text(0.96, 1.06, '(b)', fontweight='bold', color='k', 
            transform=ax.transAxes)

for i in np.arange(0,modes-1,1):
        if i<=1:
            ax = fig.add_subplot(gs[1, i],sharey=ax)
            if i==1:
                ax.text(0.93, 1.06, '(d)', fontweight='bold', color='k', 
                transform=ax.transAxes)
            else:
                ax.text(0.93, 1.06, '(c)', fontweight='bold', color='k', 
                transform=ax.transAxes)
        elif i<=3:
            ax = fig.add_subplot(gs[2, i-2],sharey=ax)
            if i ==2:
                ax.text(0.93, 1.06, '(e)', fontweight='bold', color='k', 
                transform=ax.transAxes)
            else: 
                ax.text(0.93, 1.06, '(f)', fontweight='bold', color='k', 
                transform=ax.transAxes)
        elif i<=5:
            ax = fig.add_subplot(gs[3, i-4],sharey=ax)
            if i==4:
                ax.text(0.93, 1.06, '(g)', fontweight='bold', color='k', 
                transform=ax.transAxes)
            else:
                ax.text(0.93, 1.06, '(h)', fontweight='bold', color='k', 
                transform=ax.transAxes)
        elif i<=6:
            ax = fig.add_subplot(gs[4, i-5],sharey=ax)
        cax2=SVBfunc.plotbrink(ax,grid_X,grid_Z,levelsb,xlab,ylab,i,np.ma.masked_array(vals[i+1],mask=mask),modes,lat,day1)

cbar_ax = fig.add_axes([0.95, 0.35, 0.03, 0.3])
fig.colorbar(cax2, cax=cbar_ax)
cbar_ax.set_ylabel('Pressure [Pa]')
cbar_ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))


