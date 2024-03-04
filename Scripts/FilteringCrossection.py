import SVBfunc
import xarray as xr

hej=[35,54,79,120,154,194,219]  
corrind=[30.49,30.77,31.13,31.69,32.11,32.65,33.02] #[30.49,31.69,32.11] # 30.4884, 31.6852, 32.1068
varlist=['PHIHYD','UVEL','VVEL','WVEL']

#for var in varlist:
var='WVEL'
for i in range(len(hej)):
	if var=='PHIHYD':
		dirn='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL/01b_noSVB_febTS/'
		dirw='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01b_SVB_febTS/'
		
		dsw,dsn=SVBfunc.loadNetCDFs(dirw,dirn,'PHIHYD')
	else:	
		dirn='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL/01_noSVB_febTS/'
		dirw='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01_SVB_febTS/'
		
		dsw,dsn=SVBfunc.loadNetCDFs(dirw,dirn,'DYNVARS')
		
	VALfilt,VALfilttwe,VALMITpre,dist,Z,times,timestwe=SVBfunc.ExtractAndFiltCrossect(i,dsw,dsn,1,1,var,corrind)
	
	
	FILENAME='Locations/' + str(var) + str(corrind[i]) + 'no.nc'
	ds = xr.Dataset({'VAL': (("time","z","x"), VALMITpre)
			    },
			coords ={
			    "x" : dist,
			    "z" : Z,
			    "time": times
			},
			)
			
	ds.to_netcdf(FILENAME)
	
	
	FILENAMEfilt='Locations/' + str(var) + str(corrind[i]) + 'filt.nc'
	dsf = xr.Dataset({'VAL': (("time","z","x"), VALfilt),
			'VAL2': (("time2","z","x"), VALfilttwe)
			    },
			coords ={
			    "x" : dist,
			    "z" : Z,
			    "time": times,
			    "time2": timestwe
			},
			)
	dsf.to_netcdf(FILENAMEfilt)



