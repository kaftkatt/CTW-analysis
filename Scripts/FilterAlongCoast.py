import xarray as xr
import SVBfunc

for varname in ['phiHyd', 'rhoAnoma', 'dynVars']:

  pathn= str(varname)+'ACdep55noSVB.nc'
  pathw= str(varname)+'ACdep55withSVB.nc'

  dswALL = xr.open_dataset(pathw)
  dsnALL = xr.open_dataset(pathn)

  valw=dswALL.Val
  valn=dsnALL.Val

  FILENAME='/home/athelandersson/NETCDFs/' + str(varname)+ 'AC.nc'

  SVBfunc.SavingFilteredValues(valw,valn,dswALL,FILENAME,1,1)


