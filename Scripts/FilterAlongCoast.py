import xarray as xr
import SVBfunc
import numpy as np

coast='smooth'

for varname in ['PHIHYD', 'RHOAnoma', 'SALT','THETA','UVEL','VVEL','WVEL','ETAN']:

  pathn= '/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'  + str(varname)+'ACnoSVBPREFILT.nc'
  pathw= '/home/athelandersson/NETCDFs/' + str(coast) + '/' + str(varname)+'ACwithSVBPREFILT.nc'

  dswALL = xr.open_dataset(pathw)
  dsnALL = xr.open_dataset(pathn)
  lonac = dswALL.lonAC.values
  latac = dswALL.latAC.values

  lonacn = dsnALL.lonAC.values
  latacn = dsnALL.latAC.values
  
  if coast == 'original':
     start=np.where(latacn>29.91)[0][0]
  else:
     start=np.where(latacn>29.9)[0][0]
  
  startSVB=np.where(latac==latacn[start])[0][-1]
  
  valw=dswALL.Val.values
  valn=dsnALL.Val.values
  fs=1/1200
  fs2=0
  
  FILENAMEw='/home/athelandersson/NETCDFs/' + str(coast) + '/' + str(varname)+ 'AC_SVB.nc'
  FILENAMEn='/home/athelandersson/NETCDFs/' + str(coast) + '_NO/' + str(varname)+ 'AC_NoSVB.nc'
  FILENAME='/home/athelandersson/NETCDFs/' + str(coast) + '/' + str(varname)+ 'AC.nc'

  SVBfunc.SavingFilteredValues(valn,dsnALL,FILENAMEn,1,1,fs,fs2,0)
  SVBfunc.SavingFilteredValues(valw,dswALL,FILENAMEw,1,1,fs,fs2,0) # startSVB is 0 becasue we want the whole section of the values without and with bay separately
  SVBfunc.SavingFilteredValues(valw[:,startSVB:]-valn[:,start:],dswALL,FILENAME,1,1,fs,fs2,startSVB)
  print(varname)
