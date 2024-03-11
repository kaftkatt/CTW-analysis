import xarray as xr
import SVBfunc

for varname in ['PHIHYD', 'RHOAnoma', 'SALT','THETA','UVEL','VVEL','WVEL','ETAN']:

  pathn= '/home/athelandersson/NETCDFs/' + str(varname)+'ACnoSVBPREFILT.nc'
  pathw= '/home/athelandersson/NETCDFs/' + str(varname)+'ACwithSVBPREFILT.nc'

  dswALL = xr.open_dataset(pathw)
  dsnALL = xr.open_dataset(pathn)

  valw=dswALL.Val
  valn=dsnALL.Val

  FILENAME='/home/athelandersson/NETCDFs/' + str(varname)+ 'AC.nc'

  SVBfunc.SavingFilteredValues(valw,valn,dswALL,FILENAME,1,1)


