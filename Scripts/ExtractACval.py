# ### To save all velocity along the coast
import numpy as np
import SVBfunc
import xarray as xr
# In[141]:

dep=51
coast='original'

if coast== 'original':
  tstart=2
  ts=8
elif coast == 'smooth':
  tstart=1
  ts=9

for var in ['phiHyd', 'rhoAnoma', 'dynVars','eta']:
	dirw = '/home/athelandersson/NETCDFs/' +  str(coast) + '/'
	dirn = '/home/athelandersson/NETCDFs/' +  str(coast) + '_NO/'
	
	dsw, dsn = SVBfunc.loadNetCDFs(dirw, dirn, var,tstart)
	
	name='Val'
	name2='lonAC'
	name3='latAC'
	
	LAT = dsw[0].YC
	LON = dsw[0].XC - 360
	
	Zdyn = dsw[0].Z
	hFacCw = dsw[0].hFacC
	hFacCn = dsn[0].hFacC
	hFacCusew=hFacCw.values
	hFacCusen=hFacCn.values

	if var == 'dynVars': 
		for VAR in ['UVEL','VVEL','WVEL','SALT','THETA']:
			pathn= '/home/athelandersson/NETCDFs/' +  str(coast) + '_NO/' + str(VAR)+ 'ACnoSVBPREFILT.nc'
			pathw= '/home/athelandersson/NETCDFs/' +  str(coast) + '/' + str(VAR)+'ACwithSVBPREFILT.nc'
			coordPathW='/home/athelandersson/CTW-analysis/Files/smooth/PAC/pac505mSVB.nc'
			coordPathN='/home/athelandersson/CTW-analysis/Files/smooth/PAC/pac505mNoSVB.nc'
			for l in np.arange(0,ts,1):
				
				exec (f'W=dsw[l].{VAR}[:,51,:,:].values')
				exec (f'N=dsn[l].{VAR}[:,51,:,:].values')
				
				TIMEdyn=dsw[l].time.astype(int).values*1e-9
				ntdyn = np.size(TIMEdyn)
				
				val, lon_fix,lat_fix,dist_cummul=SVBfunc.varsalongcoasts(LAT,LON,Zdyn,W,ntdyn,dep,l,VAR,coordPathW)
				valn,lon_fixn,lat_fixn,dist_cummuln=SVBfunc.varsalongcoasts(LAT,LON,Zdyn,N,ntdyn,dep,l,VAR,coordPathN)
				    
				if l == 0: valOUTw = val;valOUTn = valn;time = TIMEdyn
				else: valOUTw=np.concatenate((valOUTw,val),axis=0); valOUTn=np.concatenate((valOUTn,valn),axis=0); time=np.concatenate((time,TIMEdyn),axis=0)
				print(str(l))
			
			
			
			dsnALL = xr.Dataset({name: (("time","x"), valOUTn),
					 name2:(("x"), lon_fixn),
					 name3:(("x"), lat_fixn)
					    },
				   coords ={
					 "x" : dist_cummuln,
					 "time": time,
				      },
				       )
			dswALL = xr.Dataset({name: (("time","x"), valOUTw),
					 name2:(("x"), lon_fix),
					 name3:(("x"), lat_fix)
					    },
				   coords ={
					 "x" : dist_cummul,
					 "time": time,
				      },
				       )
			
			dsnALL.to_netcdf(pathn)
			dswALL.to_netcdf(pathw)
	elif var == 'eta': 
		VAR='ETAN'
		pathn= '/home/athelandersson/NETCDFs/' +  str(coast) + '_NO/' + str(VAR)+ 'ACnoSVBPREFILT.nc'
		pathw= '/home/athelandersson/NETCDFs/' +  str(coast) + '/' + str(VAR)+'ACwithSVBPREFILT.nc'
		coordPathW='/home/athelandersson/CTW-analysis/Files/smooth/PAC/pac20mSVB.nc'
		coordPathN='/home/athelandersson/CTW-analysis/Files/smooth/PAC/pac20mNoSVB.nc'
		for l in np.arange(0,ts,1):
			
			exec (f'W=dsw[l].{VAR}[:,:,:].values')
			exec (f'N=dsn[l].{VAR}[:,:,:].values')
			
			TIMEdyn=dsw[l].time.astype(int).values*1e-9
			nt = np.size(TIMEdyn)
			
			val, lon_fix,lat_fix,dist_cummul=SVBfunc.varsalongcoasts(LAT,LON,Zdyn,W,nt,dep,l,VAR,coordPathW)
			valn,lon_fixn,lat_fixn,dist_cummuln=SVBfunc.varsalongcoasts(LAT,LON,Zdyn,N,nt,dep,l,VAR,coordPathN)
			    
			if l == 0: valOUTw = val;valOUTn = valn;time = TIMEdyn
			else: valOUTw=np.concatenate((valOUTw,val),axis=0); valOUTn=np.concatenate((valOUTn,valn),axis=0); time=np.concatenate((time,TIMEdyn),axis=0)
			print(str(l))
			
		
		dsnALL = xr.Dataset({name: (("time","x"), valOUTn),
				 name2:(("x"), lon_fixn),
				 name3:(("x"), lat_fixn)
				    },
			   coords ={
				 "x" : dist_cummuln,
				 "time": time,
			      },
			       )
		dswALL = xr.Dataset({name: (("time","x"), valOUTw),
				 name2:(("x"), lon_fix),
				 name3:(("x"), lat_fix)
				    },
			   coords ={
				 "x" : dist_cummul,
				 "time": time,
			      },
			       )
		
		dsnALL.to_netcdf(pathn)
		dswALL.to_netcdf(pathw)

	else: 
		if var == 'phiHyd':
			VAR='PHIHYD'
		elif var == 'rhoAnoma': 
			VAR='RHOAnoma'
			
		pathn= '/home/athelandersson/NETCDFs/' +  str(coast) + '_NO/' + str(VAR)+ 'ACnoSVBPREFILT.nc'
		pathw= '/home/athelandersson/NETCDFs/' +  str(coast) + '/' + str(VAR)+'ACwithSVBPREFILT.nc'
		coordPathW='/home/athelandersson/CTW-analysis/Files/smooth/PAC/pac505mSVB.nc'
		coordPathN='/home/athelandersson/CTW-analysis/Files/smooth/PAC/pac505mNoSVB.nc'
		for l in np.arange(0,ts,1):
			
			exec (f'W=dsw[l].{VAR}[:,51,:,:].values')
			exec (f'N=dsn[l].{VAR}[:,51,:,:].values')
			
			TIMEdyn=dsw[l].time.astype(int).values*1e-9
			ntdyn = np.size(TIMEdyn)
			
			val, lon_fix,lat_fix,dist_cummul=SVBfunc.varsalongcoasts(LAT,LON,Zdyn,W,ntdyn,dep,l,VAR,coordPathW)
			valn,lon_fixn,lat_fixn,dist_cummuln=SVBfunc.varsalongcoasts(LAT,LON,Zdyn,N,ntdyn,dep,l,VAR,coordPathN)
			    
			if l == 0: valOUTw = val;valOUTn = valn;time = TIMEdyn
			else: valOUTw=np.concatenate((valOUTw,val),axis=0); valOUTn=np.concatenate((valOUTn,valn),axis=0); time=np.concatenate((time,TIMEdyn),axis=0)
			print(str(l))
			
		
		dsnALL = xr.Dataset({name: (("time","x"), valOUTn),
				 name2:(("x"), lon_fixn),
				 name3:(("x"), lat_fixn)
				    },
			   coords ={
				 "x" : dist_cummuln,
				 "time": time,
			      },
			       )
		dswALL = xr.Dataset({name: (("time","x"), valOUTw),
				 name2:(("x"), lon_fix),
				 name3:(("x"), lat_fix)
				    },
			   coords ={
				 "x" : dist_cummul,
				 "time": time,
			      },
			       )
		
		dsnALL.to_netcdf(pathn)
		dswALL.to_netcdf(pathw)
	
