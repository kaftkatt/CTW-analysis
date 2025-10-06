import SVBfunc
import numpy as np

coast='original'

if coast == 'smooth':
	startday=1
elif coast == 'original':
	startday=2

dirn='/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
dirw='/home/athelandersson/NETCDFs/' + str(coast) + '/'

dsw,dsn=SVBfunc.loadNetCDFs(dirw,dirn,'dynVars',startday)

for i in np.arange(len(dsw)):
	print(f'starting on {i}')
	ashorew,ashoren,cshorew,cshoren,t,Z,LON,LAT=SVBfunc.rotatedvelocities(dsw[i],dsn[i])
	print(f'Finished rotation {i}')
	varname='velrot'
	FILENAMEw='/home/athelandersson/NETCDFs/' + str(coast) + '/' + str(varname) + 'withSVB' + str(i + 1) + '_' + str(2 + i) + '.nc'
	FILENAMEn='/home/athelandersson/NETCDFs/' + str(coast) + '_NO/' + str(varname) + 'noSVB' + str(i + 1) + '_' + str(2 + i) + '.nc'
	SVBfunc.createDataSet(ashorew,cshorew,t,Z,LON,LAT,FILENAMEw)
	SVBfunc.createDataSet(ashoren,cshoren,t,Z,LON,LAT,FILENAMEn)

