import SVBfunc
import xarray as xr
import numpy as np

hej=[35,54,79,120,154,194,219]  
corrind=[30.49,30.77,31.13,31.69,32.11,32.65,33.02] #[30.49,31.69,32.11] # 30.4884, 31.6852, 32.1068
varlist=['PHIHYD','UVEL','VVEL','WVEL']

i=0
var=varlist[0]
FILENAMEIN='/home/athelandersson/CTW-analysis/Files/Locations/' + str(var) + str(corrind[i]) + 'no.nc'

ds=xr.open_dataset(FILENAMEIN)

VALMITpre=ds.VAL.values

fs=1/1200
fs2=0

Z=dsw[0].Zl.values

matfile=loadmat('/home/athelandersson/CTW-analysis/Files/BT_P2.mat')
x,dep,indXlon,indYlat=matfile['dist'],matfile['d'],matfile['indexXlon'],matfile['indexYlat']

print('Filtering begins')

for d in np.arange(np.size(VALMITpre,2)):
    VALdif,VALfiltout,VALfiltAll,inds = SVBfunc.FiltDetrend(VALMITpre[:,:,d],filt,detrend,fs,fs2)
    VALfilt[:,:,d]=VALfiltAll

FILENAMEfilt='/home/athelandersson/CTW-analysis/Files/Locations/' + str(var) + str(corrind[i]) + 'filt.nc'

dsf = xr.Dataset({'VAL': (("time","z","x"), VALfilt)
        },
    coords ={
        "x" : dist,
        "z" : Z,
        "time": times
    },
    )
dsf.to_netcdf(FILENAMEfilt)
