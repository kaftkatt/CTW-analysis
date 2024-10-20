import SVBfunc
import xarray as xr
COASTS=['original','smooth']
NUM=['1','30']
coast='original'
nr=1
#for nr in NUM: 
#for coast in COASTS:
hej=[35,54,79,120,154,194,219]  
corrind=[30.49,30.77,31.13,31.69,32.11,32.65,33.02] #[30.49,31.69,32.11] # 30.4884, 31.6852, 32.1068
varlist= ['ashore','PHIHYD','WVEL'] #['ashore','cshore']
varlistLONGNAME= ['Crosshore velocity','Hydrostatic Pressure Pot.(p/rho) Anomaly','U-Velocity','V-Velocity','Vertical Velocity'] #['Alongshore velocity','Crosshore velocity']

dirn='/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
dirw='/home/athelandersson/NETCDFs/' + str(coast) + '/'

if coast == 'smooth':
	startday=1
elif coast == 'original':
	startday=2
k=0
for var in varlist:
	varname=varlistLONGNAME[k]
	k=k+1
	if var=='PHIHYD':			
		dsw,dsn=SVBfunc.loadNetCDFs(dirw,dirn,'phiHyd',startday)
		units=dsw[0].PHIHYD.units
	else:	
		dsw,dsn=SVBfunc.loadNetCDFs(dirw,dirn,'dynVars',startday)
		units='m/s'
	
	for i in range(len(hej)):
		VALfilt,VALMITpre,x,dep,lon,lat,deg,Z,times=SVBfunc.CrossectExctraction(i,dsw,dsn,1,1,var,corrind,coast)
	
		FILENAME='/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/Locations/' + str(var) + str(corrind[i]) + '.nc'
		
		title = 'Crossection of' + varname
		description = 'Extracted crossection from MITgcm output. Using the provided longitude and latitude. If crosshore or alongshore velocity is specified they have been calculated from U and V velocity using the provided angle.' 
		SVBfunc.create_descriptive_file(times, Z, x, dep,lon,lat,deg, VALMITpre, VALfilt, varname, var, units, FILENAME, title, description)
	
	
	
