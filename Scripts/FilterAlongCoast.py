import xarray as xr
import SVBfunc

coast='original'

for varname in ['PHIHYD', 'RHOAnoma', 'SALT','THETA','UVEL','VVEL','WVEL','ETAN']:

  pathn= '/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'  + str(varname)+'ACnoSVBPREFILT.nc'
  pathw= '/home/athelandersson/NETCDFs/' + str(coast) + '/' + str(varname)+'ACwithSVBPREFILT.nc'

  dswALL = xr.open_dataset(pathw)
  dsnALL = xr.open_dataset(pathn)

  valw=dswALL.Val
  valn=dsnALL.Val
  fs=1/1200
  fs2=0
  
  FILENAMEw='/home/athelandersson/NETCDFs/' + str(coast) + '/' + str(varname)+ 'AC_SVB.nc'
  FILENAMEn='/home/athelandersson/NETCDFs/' + str(coast) + '_NO/' + str(varname)+ 'AC_NoSVB.nc'

  SVBfunc.SavingFilteredValues(valn,dsnALL,FILENAMEn,1,1,fs,fs2)
  SVBfunc.SavingFilteredValues(valw,dswALL,FILENAMEw,1,1,fs,fs2)


