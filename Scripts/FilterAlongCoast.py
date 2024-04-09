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
  
  FILENAME='/home/athelandersson/NETCDFs/' + str(coast) + '/' + str(varname)+ 'AC.nc'

  SVBfunc.SavingFilteredValues(valw,valn,dswALL,FILENAME,1,1,fs,fs2)


