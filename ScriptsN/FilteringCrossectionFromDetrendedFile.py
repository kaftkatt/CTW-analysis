import SVBfunc
import xarray as xr
import numpy as np

hej=[35,54,79,120,154,194,219]  
corrind=[30.49,30.77,31.13,31.69,32.11,32.65,33.02] #[30.49,31.69,32.11] # 30.4884, 31.6852, 32.1068
varlist=['PHIHYD','UVEL','VVEL','WVEL']
for var in varlist:
    for i in range(len(corrind)):
        
        FILENAMEIN='/home/athelandersson/CTW-analysis/Files/Locations/' + str(var) + str(corrind[i]) + 'no.nc'
        
        ds=xr.open_dataset(FILENAMEIN)
        
        VALMITpre=ds.VAL.values
        
        fs=1/1200
        fs2=0
        
        filt=1
        detrend=0
        VALfilt=np.zeros(np.shape(VALMITpre))
        
        Z=ds.z.values
        dist=ds.x.values
        times=ds.time.values
        
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
