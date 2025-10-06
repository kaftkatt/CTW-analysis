from SVBfunc import loadNetCDFs

import xarray as xr
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
coast='smooth'
tstart=2

dirn = '/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
dirw = '/home/athelandersson/NETCDFs/' + str(coast) + '/'
dswSurf, dsnSurf = loadNetCDFs(dirw, dirn, 'eta',tstart)
LON=dswSurf[0].XC-360
LAT=dswSurf[0].YC
indy=np.where(LAT>28)[0][0]
indx=np.where(LON>-121)[0][0]
etaplot=[]
TIME=[]
for i in range(len(dswSurf)):
    for k in range(len(dswSurf[i].ETAN[:,0,0])):
        etaplot.append(dswSurf[i].ETAN[k,indy,indx].values-dsnSurf[i].ETAN[k,indy,indx].values)
        TIME.append(dswSurf[i].time[k].astype(int).values*1e-9)

plt.plot(np.asarray(TIME)/(60*60*24),etaplot)
plt.xlabel('Time [days]')
plt.ylabel('SSH [m]')
plt.show()
plt.savefig('/home/athelandersson/CTW-analysis/Figures/figSSHsig.png')
