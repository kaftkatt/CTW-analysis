import SVBfunc
import xarray as xr
COASTS=['original','smooth']

for coast in COASTS:
	
	hej=[35,54,79,120,154,194,219]  
	corrind=[30.49,30.77,31.13,31.69,32.11,32.65,33.02] #[30.49,31.69,32.11] # 30.4884, 31.6852, 32.1068
	varlist= ['ashore','cshore'] #['PHIHYD','UVEL','VVEL','WVEL','ashore','cshore']
	
	dirn='/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
	dirw='/home/athelandersson/NETCDFs/' + str(coast) + '/'
	 
	if coast == 'smooth':
		startday=1
	elif coast == 'original':
		startday=2

	for var in varlist:
		for i in range(len(hej)):
			if var=='PHIHYD':			
				dsw,dsn=SVBfunc.loadNetCDFs(dirw,dirn,'phiHyd',startday)
			else:	
				dsw,dsn=SVBfunc.loadNetCDFs(dirw,dirn,'dynVars',startday)
				
			VALfilt,VALMITpre,dist,Z,times=SVBfunc.ExtractAndFiltCrossectNEW(i,dsw,dsn,1,1,var,corrind,coast)
			
			
			FILENAME='/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/Locations/' + str(var) + str(corrind[i]) + 'no.nc'
			ds = xr.Dataset({'VAL': (("time","z","x"), VALMITpre)
					    },
					coords ={
				    	"x" : dist,
				    	"z" : Z,
				   	 "time": times
					},
					)
					
			ds.to_netcdf(FILENAME)
			
			
			FILENAMEfilt='/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/Locations/' + str(var) + str(corrind[i]) + 'filt.nc'
			dsf = xr.Dataset({'VAL': (("time","z","x"), VALfilt)
					    },
					coords ={
				    	"x" : dist,
				    	"z" : Z,
				    	"time": times
					},
					)
			dsf.to_netcdf(FILENAMEfilt)
	
	
	
