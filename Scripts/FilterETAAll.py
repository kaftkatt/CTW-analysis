import xarray as xr
import numpy as np
import SVBfunc

dirn='/home/athelandersson/NETCDFs/smooth_NO/'
dirw='/home/athelandersson/NETCDFs/smooth/'

dsw,dsn=SVBfunc.loadNetCDFs(dirw,dirn,'eta',1)

LAT = dsw[0].YC
LON = dsw[0].XC - 360
hFacC = dsw[0].hFacC
hfac = np.ma.masked_values(hFacC, 0)
mask = np.ma.getmask(hfac)

time12=dsw[0].time.values.astype(int)
time23=dsw[1].time.values.astype(int)
time34=dsw[2].time.values.astype(int)
time45=dsw[3].time.values.astype(int)
time56=dsw[4].time.values.astype(int)
time67=dsw[5].time.values.astype(int)
time78=dsw[6].time.values.astype(int)
time89=dsw[7].time.values.astype(int)
time910=dsw[8].time.values.astype(int)

Time=np.concatenate((time12,time23, time34, time45, time56,time67, time78,time89, time910), axis=0)#, time910), axis=0)

times=Time*1e-9

VALMITpre=np.zeros((len(Time),np.size(mask[0,:,:]),np.size(mask[0,:,:])))

for tt in np.arange(0,9,1):
	
	VALMIT=np.zeros(np.shape(dsw[tt].ETAN))
		    
	for t in np.arange(0,len(dsw[tt].ETAN[:,1,1]),1):
		VALb=dsw[tt].ETAN[t,:,:].values
		VALn=dsn[tt].ETAN[t,:,:].values
		VALmit=VALb-VALn
		VALMIT[t,:,:]=VALmit
	
	VALMITpre[len(dsw[tt-1].ETAN[:,1,1])*tt:len(VALMIT[:,1,1])*(tt+1),:,:]=VALMIT
	print('Day '+str(tt+2))


VALfilt=np.zeros(np.shape(VALMITpre))

fs=1/1200
fs2=0

filt=1
detrend=1

print('Filtering begins')
for d in np.arange(np.size(VALMITpre,2)):
	VALdif,VALfiltout,VALfiltAll,inds = SVBfunc.FiltDetrend(VALMITpre[:,:,d],filt,detrend,fs,fs2)
	VALfilt[:,:,d]=VALfiltAll

r



FILENAMEfilt='/home/athelandersson/NETCDFs/ETA.nc'
dsf = xr.Dataset({"VALfilt": (("time","y","x"), VALfilt),
		"VAL": (("time","y","x"), VALMITpre)
		    },
		coords ={
		    "x" : LON.values,
		    "y" : LAT.values,
		    "time": times
		},
		)
dsf.to_netcdf(FILENAMEfilt)
