import SVBfunc
import xarray as xr
from scipy.io import loadmat
import numpy as np

COASTS=['original','smooth']
NUM=['1','30']
coast='smooth'
nr=1
all=1
#for nr in NUM: 
#for coastin COASTS:
pathDIST= '/home/athelandersson/NETCDFs/' + str(coast) + '/ETANAC.nc'
dsDIST= xr.open_dataset(pathDIST)

if all==1:  
	matfile=loadmat( '/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/BT_PALL_MovAv2025.mat')
	x=matfile['dist'][0]
	hej=range(len(x))
	latout=matfile['lat'][0]
	lonout=matfile['lon'][0]
	
	lonin=[]
	latin=[]
	for i in range(len(x)):
		lonin.append(lonout[i][0][0])
		latin.append(latout[i][0][0])
	
	distAC=np.zeros(len(x)) 
	
	p=0
	dist_array = np.zeros(len(x)-1)
	for ii in np.arange(1,len(x),1):
		lat1 = latin[ii-1]
		lon1 = lonin[ii-1]
		lat2 = latin[ii]
		lon2 = lonin[ii]
		p=p+1
		dist_array[p-1]=  SVBfunc.haversine(lat1, lon1, lat2, lon2)

	distAC[1:] = np.cumsum(dist_array)
else: 
	hej=[35,54,79,120,154,194,219] 
	distAC=dsDIST.dist[hej].values

varlist=['ashore','PHIHYD','WVEL'] #['ashore','cshore']
varlistLONGNAME= ['Alongshore velocity','Hydrostatic Pressure Pot.(p/rho) Anomaly','Vertical Velocity'] #['Alongshore velocity','Crosshore velocity']

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
	VALMIT=[]
	VALFILT=[]
	VALMITw=[]
	VALFILTw=[]
	VALMITn=[]
	VALFILTn=[]
	X=[]
	DEP=[]
	LON=[]
	LAT=[]
	DEG=[]
	Zout=[]

	for i in range(len(hej)):
		
		VALfilt,VALMITpre,VALfiltn,VALMITpren,VALfiltw,VALMITprew,x,dep,lon,lat,deg,Z,times=SVBfunc.CrossectExctraction(i,dsw,dsn,1,1,var,all,coast)
		
		VALMIT.append(VALMITpre[:,:,:58])
		VALFILT.append(VALfilt[:,:,:58])
		VALMITw.append(VALMITprew[:,:,:58])
		VALFILTw.append(VALfiltw[:,:,:58])
		VALMITn.append(VALMITpren[:,:,:58])
		VALFILTn.append(VALfiltn[:,:,:58])
		X.append(x[:58])
		DEP.append(dep[:58])
		LON.append(lon[:58])
		LAT.append(lat[:58])
		DEG.append(deg)
		Zout.append(Z)
	
	FILENAME='/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/Locations/' + str(var) + '.nc'
	
	title = 'Crossection of ' + varname
	description = 'Extracted crossection from MITgcm output. Using the provided longitude and latitude. If crosshore or alongshore velocity is specified they have been calculated from U and V velocity using the provided angle.' 
	SVBfunc.create_descriptive_file(times, Zout, X, DEP,LON,LAT,DEG, VALMIT, VALFILT,VALMITn,VALFILTn,VALMITw,VALFILTw,distAC,varname, var, units, FILENAME, title, description)

	
	
