import xarray as xr
import numpy as np
import SVBfunc
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import cmocean
from os.path import exists
from matplotlib.ticker import FormatStrFormatter


hej=[35,54,79,120,154,194,219]
corrinds=[30.49,30.77,31.13,31.69,32.11,32.65,33.02]
ik=3
day1=0

lat=corrinds[ik]

var='PHIHYD'
u,v,w,r,p,k,omega,epe,eke,xgr,zgr = SVBfunc.openBrink(corrinds[ik])

if var == 'UVEL':
	valsin=u
	clabelB='u-Velocity [cms$^{-1}$]' 
elif var == 'VVEL':
	valsin=v
	clabelB='v-Velocity [cms$^{-1}$]' 
elif var == 'WVEL':
	valsin=w
	clabelB='w-Velocity [cms$^{-1}$]'
elif var == 'PHIHYD':
	valsin=p
	clabelB='Pressure [Pa]'

modes=0
modenr=[]
vals=[]
for i in range(len(valsin)):
	if np.any(valsin[i]!=0):
		modenr.append(i)
		modes=modes+1
		if np.logical_or(var == 'WVEL',var== 'UVEL'):
			vals.append(valsin[i].imag)
		else:
			vals.append(valsin[i])


vals=np.array(vals)
		
vminb=-np.nanmax(abs(vals))
vmaxb=np.nanmax(abs(vals))

levelsb=np.linspace(vminb,vmaxb,15)

xlab='Cross-shore distance [km]'
ylab='Depth [m]'


fig = plt.figure()
if modes<=4:
    gs = GridSpec(nrows=2, ncols=2, height_ratios=[1,1],hspace=0.35)
elif modes<=6:
    gs = GridSpec(nrows=3, ncols=2, height_ratios=[1,1,1],hspace=0.35)
elif modes<=8:
    gs = GridSpec(nrows=4, ncols=2, height_ratios=[1,1,1,1],hspace=0.35)



    
for i in np.arange(0,modes-1,1):
        if i<=1:
            ax = fig.add_subplot(gs[0, i])
            if i==1:
                ax.text(0.93, 1.06, '(b)', fontweight='bold', color='k', 
                transform=ax.transAxes)
            else:
                ax.text(0.93, 1.06, '(a)', fontweight='bold', color='k', 
                transform=ax.transAxes)
        elif i<=3:
            ax = fig.add_subplot(gs[1, i-2])
            if i ==2:
                ax.text(0.93, 1.06, '(d)', fontweight='bold', color='k', 
                transform=ax.transAxes)
            else: 
                ax.text(0.93, 1.06, '(e)', fontweight='bold', color='k', 
                transform=ax.transAxes)
        elif i<=5:
            ax = fig.add_subplot(gs[2, i-4])
            if i==4:
                ax.text(0.93, 1.06, '(f)', fontweight='bold', color='k', 
                transform=ax.transAxes)
            else:
                ax.text(0.93, 1.06, '(g)', fontweight='bold', color='k', 
                transform=ax.transAxes)
        elif i<=6:
            ax = fig.add_subplot(gs[3, i-5])
        cax2=SVBfunc.plotbrink(ax,xgr,zgr,levelsb,xlab,ylab,modenr[i+1],i,vals[i+1],modes,lat,day1)

cbar_ax = fig.add_axes([0.92, 0.4, 0.03, 0.45])
fig.colorbar(cax2, cax=cbar_ax)
cbar_ax.set_ylabel(clabelB)
cbar_ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))

plt.show()
