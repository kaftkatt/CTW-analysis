import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import pylab as pl
import scipy.signal as sig
import scipy.io as sio
import scipy.interpolate as sciint
from xmitgcm import open_mdsdataset
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
from os.path import exists
from scipy.io import loadmat

from math import radians, cos, sin, asin, sqrt, atan, degrees

def loadNetCDFs(dirw,dirn,varname,startday):
    dsw = []
    dsn = []
    for i in np.arange(startday, 10, 1):
	
        pathn = dirn + str(varname) + 'noSVB' + str(i) + '_' + str(1 + i) + '.nc'
        pathw = dirw + str(varname) + 'withSVB' + str(i) + '_' + str(1 + i) + '.nc'

        dswin = xr.open_dataset(pathw)
        dsnin = xr.open_dataset(pathn)

        dsw.append(dswin)
        dsn.append(dsnin)

    return dsw, dsn

def loadWVEL(dsw,dsn):
    var23=dsw[0].Ww.values
    var34=dsw[1].Ww.values
    var45=dsw[2].Ww.values
    var56=dsw[3].Ww.values
    var67=dsw[4].Ww.values
    var78=dsw[5].Ww.values
    var89=dsw[6].Ww.values
    var910=dsw[7].Ww.values

    Ww=np.concatenate((var23, var34, var45, var56,var67,var78,var89, var910), axis=0)

    var23n=dsn[0].Wn.values
    var34n=dsn[1].Wn.values
    var45n=dsn[2].Wn.values
    var56n=dsn[3].Wn.values                
    var67n=dsn[4].Wn.values
    var78n=dsn[5].Wn.values
    var89n=dsn[6].Wn.values
    var910n=dsn[7].Wn.values

    Wn=np.concatenate((var23n, var34n, var45n, var56n,var67n,var78n,var89n, var910n), axis=0) 
   
    time23=dsw[0].time.values.astype(int)
    time34=dsw[1].time.values.astype(int)
    time45=dsw[2].time.values.astype(int)
    time56=dsw[3].time.values.astype(int)
    time67=dsw[4].time.values.astype(int)
    time78=dsw[5].time.values.astype(int)
    time89=dsw[6].time.values.astype(int)
    time910=dsw[7].time.values.astype(int)
    
    Time=np.concatenate((time23, time34, time45, time56,time67, time78,time89, time910), axis=0)#, time910), axis=0)
    
    times=Time*1e-9
      
    return Ww,Wn,times


def createNetCDF(coast,prefix, varname):
	
	levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
		11, 12, 13, 14, 15, 16, 17,
		18, 19, 20, 21, 22, 23, 24, 25,
		26, 27, 28, 29, 30, 31,
		32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
		45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
		58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
		74, 79, 84, 89, 94, 99, ]
	
	if coast == 'smooth':
		pathw = '/data/SO2/sio-kramosmusalem/exp11_512x612x100_smooth_SVB/01_febTS_1000x'
		pathn = '/data/SO2/sio-kramosmusalem/exp11_512x612x100_smooth/01_febTS_1000x'
		pathnN = '/home/athelandersson/NETCDFs/smooth_NO/'
		pathwN = '/home/athelandersson/NETCDFs/smooth/'
		for whatdaystart,whatdayfinish in zip(np.arange(1,10,1),np.arange(2,11,1)):	
			dayarr = np.arange(whatdaystart * 24 * 60 * 2, whatdayfinish * 24 * 60 * 2, 40)
			day = dayarr.tolist()	
			dsw = open_mdsdataset(pathw, pathw, prefix=[prefix], default_dtype='>f4', levels=levels, iters=day, delta_t=30)
			dsn = open_mdsdataset(pathn, pathn, prefix=[prefix], default_dtype='>f4', levels=levels, iters=day, delta_t=30)
			
			pathwNEW = pathwN + str(varname) + 'withSVB' + str(whatdaystart) + '_' + str(whatdayfinish) + '.nc'
			pathnNEW = pathnN + str(varname) + 'noSVB' + str(whatdaystart) + '_' + str(whatdayfinish) + '.nc'
			
			dsw.to_netcdf(path=pathwNEW)
			print('Done with SVB, day ' + str(whatdaystart) + ' - ' + str(whatdayfinish))
			dsn.to_netcdf(path=pathnNEW)
			print('Done without SVB, day ' + str(whatdaystart) + ' - ' + str(whatdayfinish))
	elif coast == 'originial':
		if np.logical_or(varname == 'phiHyd',varname == 'rhoAnoma'):
			pathw = '/data/SO2/sio-kramosmusalem/exp06_512x612x100_ORL_SVB/01b_SVB_febTS_output/'
			pathn = '/data/SO2/sio-kramosmusalem/exp06_512x612x100_ORL/01b_noSVB_febTS/'
		else:
			pathw = '/data/SO2/sio-kramosmusalem/exp06_512x612x100_ORL_SVB/01_SVB_febTS_output/'
			pathn = '/data/SO2/sio-kramosmusalem/exp06_512x612x100_ORL/01_noSVB_febTS/'
		pathnN = '/home/athelandersson/NETCDFs/original_NO/'
		pathwN = '/home/athelandersson/NETCDFs/original/'
		for whatdaystart,whatdayfinish in zip(np.arange(2,10,1),np.arange(3,11,1)):	
			dayarr = np.arange(whatdaystart * 24 * 60, whatdayfinish * 24 * 60, 20)
			day = dayarr.tolist()	
			dsw = open_mdsdataset(pathw, pathw, prefix=[prefix], default_dtype='>f4', levels=levels, iters=day)
			dsn = open_mdsdataset(pathn, pathn, prefix=[prefix], default_dtype='>f4', levels=levels, iters=day)
			
			pathwNEW = pathwN + str(varname) + 'withSVB' + str(whatdaystart) + '_' + str(whatdayfinish) + '.nc'
			pathnNEW = pathnN + str(varname) + 'noSVB' + str(whatdaystart) + '_' + str(whatdayfinish) + '.nc'
			
			dsw.to_netcdf(path=pathwNEW)
			print('Done with SVB, day ' + str(whatdaystart) + ' - ' + str(whatdayfinish))
			dsn.to_netcdf(path=pathnNEW)
			print('Done without SVB, day ' + str(whatdaystart) + ' - ' + str(whatdayfinish))


def findlonlat(hFacCuse,d,var):
 # ## Doing the picking out 
    nx = 512
    ny = 612
    ind30=170+50 #30° N
    indlon = 50 # Lon of land at N boundary

    lon_inds_off = np.argmax(np.squeeze(hFacCuse[d,:,::-1]), axis=1)


    lon_inds = np.ones_like(lon_inds_off[ind30:])*nx - lon_inds_off[ind30:]
    lat_inds = np.ones_like(lon_inds)*ind30 + np.arange(len(lon_inds))

    lat_inds_off = np.argmax(np.squeeze(hFacCuse[d,::-1,:]), axis=0)



    lat_inds_2 = np.ones_like(lat_inds_off[indlon:])*ny - lat_inds_off[indlon:]
    lon_inds_2 = np.ones_like(lat_inds)*indlon + np.arange(len(lat_inds))

    p=0

    k=np.ones_like(lon_inds)
    l=np.ones_like(lat_inds)

    for i in np.arange(len(lon_inds_2)):
        if lon_inds_2[i] not in lon_inds:
            p=p+1
            k[p]=lon_inds_2[i]
            l[p]=lat_inds_2[i]
        
    k=k[k!=1]
    l=l[l!=1]
    
    # Fill the holes
    holearray=[]
    for i in np.arange(len(l)):
        if l[i] not in holearray:
            if l[i] in lat_inds:
                indexpts=np.where(lat_inds==l[i])[0]
                indexs=np.where(l==lat_inds[indexpts[0]])
                lon_inds=np.insert(lon_inds,indexpts[0],k[indexs])
                lat_inds=np.insert(lat_inds,indexpts[0],l[indexs])
                holearray.append(l[i+1])
    
    
    #Organize where the values are in the wrong order
    if var==1:
        lon_lat=np.zeros((2,len(lon_inds[332:436])))
        lon_lat[0,:]=lat_inds[332:436]
        lon_lat[1,:]=lon_inds[332:436]

        inds=lon_lat[1, :].argsort()
        indices=np.flip(inds)

        latis=lon_lat[0,indices]
        lonis=lon_lat[1,indices]

        lat_inds[332:436]=latis
        lon_inds[332:436]=lonis

    print("Done")
    
    return(lat_inds, lon_inds)


def haversine(lon1, lat1, lon2, lat2):

# ######   Calculate the great circle distance in kilometers between two points on the earth (specified in decimal degrees)
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r

def loadingcoastpts(hFacCuse,Z,LAT,LON,prevalue,nt,t,var):

    dep=0
    lat_fix, lon_fix=findlonlat(hFacCuse,dep,0)


    dist_array = np.zeros(len(lon_fix))
    p=0
    for jj,ii in zip(lon_fix, lat_fix):
        lat1 = LAT[ii-2]
        lon1 = LON[jj-2]
        lat2 = LAT[ii-1]
        lon2 = LON[jj-1]
        p=p+1
        dist_array[p-1]=  haversine(lat1, lon1, lat2, lon2)
    
    dist_cummul = np.cumsum(dist_array)
    
    value=prevalue.values
    val = np.zeros((nt,len(lon_fix)))
    
    p=0
    
    for ii,jj in zip(lon_fix, lat_fix):
        p=p+1
        val[:,p-1] = value[:,jj-1,ii-1]

    
    if t==7:
       plotpointsAC(LON,LAT,lon_fix,lat_fix,var)
    
    return val[:,val[0,:]!=0],lon_fix[val[0,:] != 0],lat_fix[val[0,:] != 0],dist_cummul[val[0,:] != 0]


#For dynvars and pressure, just don't forget to change the depth
def varsalongcoasts(hFacCuse,LAT,LON,Z,valvar,nt,dep,plot,var):

    
    lat_fix, lon_fix=findlonlat(hFacCuse,dep,1)
    
    dist_array = np.zeros(len(lon_fix[5:-5]))
    p=0
    for ii,jj in zip(lon_fix[5:-5],lat_fix[5:-5]):
        lat1 = LAT[jj-2]
        lon1 = LON[ii-2]
        lat2 = LAT[jj-1]
        lon2 = LON[ii-1]
        p=p+1
        dist_array[p-1]=  haversine(lat1, lon1, lat2, lon2)

    dist_cummul = np.cumsum(dist_array)
    rangedep=np.arange(0,71,10)
    val = np.zeros((nt,len(lon_fix[5:-5]))) 
    dd=0

    for h in range(nt):
        p=0
        for n in np.arange(5,len(lon_fix)-5,1):
            p=p+1
            
            value=np.mean(valvar[h-1,lat_fix[n-5:n+5]-1,lon_fix[n-5:n+5]-1])
            val[h-1,p-1] = value

    lon_fix=lon_fix[5:-5]
    lat_fix=lat_fix[5:-5]
   
    if plot==7:
    	plotpointsAC(LON,LAT,lon_fix,lat_fix,var)
    
    return val[:,val[0,:]!=0],lon_fix[val[0,:]!=0],lat_fix[val[0,:]!=0],dist_cummul[val[0,:]!=0]
    
# When you want a diagonal across the area    
def varsalldepths(hFacCuse,LAT,LON,Z,valvar,name,name2,name3,nt,t,FILENAME):
    Zdyn=Z[:71].values
    
    
    latss=LAT.values
    leo=np.arange(4,590,4)
    LATSpre=np.delete(latss,np.arange(0,10,1))
    LATSpre=np.delete(LATSpre,leo)

    lat_ind=np.zeros(455)
    for i in np.arange(0,455,1):
        indexx=np.where(latss==LATSpre[i])
        lat_ind[i]=indexx[0]
    
    lon_ind=np.flip(np.arange(10,445,1))
    lat_ind=lat_ind[10:-10].astype(int)


    dist_array = np.zeros(len(lon_ind))
    p=0
    for ii,jj in zip(lon_ind,lat_ind):
        lat1 = LAT[jj-2]
        lon1 = LON[ii-2]
        lat2 = LAT[jj-1]
        lon2 = LON[ii-1]
        p=p+1
        dist_array[p-1]=  haversine(lat1, lon1, lat2, lon2)

    dist_cummul = np.cumsum(dist_array)
    rangedep=np.arange(0,71,10)
    val = np.zeros((nt,len(Zdyn[rangedep]),len(lon_ind))) 
    dd=0
    
    for dep in rangedep:
        dd=dd+1
        for h in range(nt):
            p=0
            for ii,jj in zip(lon_ind,lat_ind):
                p=p+1
                value=valvar[h-1,dep,jj-1,ii-1]
                val[h-1,dd-1,p-1] = value
            

    ds = xr.Dataset({name: (("time","z","x"), val),                 
                     name2:(("x"), lon_ind),
                     name3:(("x"), lat_ind)
                    },
                coords ={
                    "x" : dist_cummul,
                    "z" : Zdyn[rangedep],
                    "time": t,
                },
                )
    ds.to_netcdf(FILENAME)

def sgolay2d ( z, window_size, order, derivative=None):
    """
    """
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0
    
    if  window_size % 2 == 0:
        raise ValueError('window_size must be odd')
    
    if window_size**2 < n_terms:
        raise ValueError('order is too high for the window size')
    
    
    half_size = window_size // 2
    
    # exponents of the polynomial. 
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
    # this line gives a list of two item tuple. Each tuple contains 
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]
    
    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat( ind, window_size )
    dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )
    
    # build matrix of system of equation
    A = np.empty( (window_size**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:,i] = (dx**exp[0]) * (dy**exp[1])
        
    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = np.zeros( (new_shape) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band -  np.abs( np.flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -band ) 
    # left band
    band = np.tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z
    
    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band ) 
    
    # top right corner
    band = Z[half_size,-half_size:]
    Z[:half_size,-half_size:] = band - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band ) 
    # bottom left corner
    band = Z[-half_size:,half_size].reshape(-1,1)
    Z[-half_size:,:half_size] = band - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band ) 
    
    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -c, mode='valid')        
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid')        
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid'), scipy.signal.fftconvolve(Z, -c, mode='valid')       
  
def butter_lowpass(cutoff, fs, order):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = sig.butter(order, normal_cutoff, btype='bandpass', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = sig.filtfilt(b, a, data)
    return y

#Output every 10 min (600 seconds) until location nr 432 (day 5 in minutes 7200) then to 433 its 20 mins. So to filter the first 5 days with 1/600 fs we write VALdif[:433] and then VALdif[433:]. 
def FiltDetrend(VAL,filt,detrend,fs,fs2):
	if detrend==1:
		VALdif=sig.detrend(VAL,0)
	elif detrend==0:
		VALdif=VAL
	if filt==1:
		VALfilt=np.zeros(np.shape(VALdif))

		order = 3      
		cutoff =np.array([1/432000, 1/43200])

		if fs2 != 0:
			inds=np.append(np.arange(0,433,2),np.arange(433,792,1))
			VALfiltAll=np.zeros(np.shape(VALdif[inds,:]))

			for d in np.arange(np.size(VALdif,1)):
				data = VALdif[:433,d] 
				VALfilt[:433,d] = butter_lowpass_filter(data, cutoff, fs2, order)


			for d in np.arange(np.size(VALdif,1)):
				data = VALdif[433:,d] 
				VALfilt[433:,d] = butter_lowpass_filter(data, cutoff, fs, order)
			
			
			for d in np.arange(np.size(VALdif,1)):
				data = VALdif[inds,d] 
				VALfiltAll[:,d] = butter_lowpass_filter(data, cutoff, fs, order)
	
		else:
			VALfiltAll=np.zeros(np.shape(VALdif))
			for d in np.arange(np.size(VALdif,1)):
				data = VALdif[:,d] 
				VALfiltAll[:,d] = butter_lowpass_filter(data, cutoff, fs, order)
		
		
			Valfilt=np.zeros(np.shape(VALfiltAll))
			inds=VALfilt[:,0]==0
		
	return(VALdif,VALfilt,VALfiltAll,inds)
       	
def SavingFilteredValues(valw,valn,dsw, FILENAME,filt,detrend,fs,fs2):
	
	dist = dsw.x.values
	TIME=dsw.time.values #TIME = dsw.time.astype(int).values*1e-9
	VAL = valw.values-valn.values
	lon_ac=dsw.lonAC.values
	lat_ac=dsw.latAC.values
	
	valdet,valfilt,valfiltall,inds = FiltDetrend(VAL,filt,detrend,fs,fs2)

	ds = xr.Dataset({'Valdet': (("time","dist"),valdet),
                 	'Valfilt':(("time","dist"),valfilt),
                 	'ValfiltAll':(("time2","dist"),valfiltall),
                 	'Valorig':(("time","dist"),VAL),
                 	'lonAC':(("dist"),lon_ac),
                 	'latAC':(("dist"),lat_ac)
                 
                 	   },
                	coords ={
                    	"time": TIME,
                    	"time2": TIME[inds],
                   	 "dist":dist
                	},
                	)
	ds.to_netcdf(FILENAME)

def FFRQ(Wdif,Wfilt,timemin,dist):

    times = (timemin)*60 #in s

    t0 = 172800 # start is day 2 

    dt = 1200 # 20 min 

        
    nx = len(dist)
    nt = int(len(times)/2)
        
    psd = np.zeros((nx,nt))*np.nan
    phase = np.zeros((nx,nt))*np.nan
    
    psdfilt = np.zeros((nx,nt))*np.nan
    phasefilt = np.zeros((nx,nt))*np.nan
    
    freq =  np.zeros((nx,nt))*np.nan
    freqfilt = np.zeros((nx,nt))*np.nan
    

    for ii in np.arange(nx): #nx
        signalFFT = np.fft.rfft(Wdif[:,ii])
        signalFFTfilt = np.fft.rfft(Wfilt[:,ii])
    
        ## Get Power Spectral Density
        signalPSD = np.abs(signalFFT)
        signalPSDfilt= np.abs(signalFFTfilt)
        
        ## Get frequencies corresponding to signal 
        fftFreq = np.fft.rfftfreq(len(Wdif[:,ii]), dt)
        fftFreqfilt = np.fft.rfftfreq(len(Wfilt[:,ii]), dt)
        
        psd[ii,:] = signalPSD[1:]
   
        psdfilt[ii,:] = signalPSDfilt[1:]
        
        freq[ii,:] =  fftFreq[1:]

        freqfilt[ii,:] = fftFreqfilt[1:]
        
        
    return psd, freq,psdfilt,freqfilt




## FOR LINEAR REGRESSION -------------------------------------------------------------------------
def create_descriptive_file(t, Z, X,dep,lon,lat,deg, var, varfilt, nameLong, nameShort, units, filename, title, description): 
    
    """ This function creates a netCDF4 file for
    the given variable given the filename and 
    the time, x and z grid.
    
    :arg time: Time in seconds
    :arg Z: Depth 
    :arg X: Cross-shore distance
    :arg var: Variable in question
    :arg varfilt: Filtered version of variable in question
    :arg filename: Directory and name of netcdf file
    :arg title: Title of version
    :arg description: Details about version
    """
    
    dataset = Dataset(filename, 'w')
    file_time = dataset.createDimension('time', len(t))
    file_z = dataset.createDimension('z', len(Z) )
    file_x = dataset.createDimension('x', len(X) )

    file_X = dataset.createVariable('X', 'f8', ('x'))
    file_lon = dataset.createVariable('LON', 'f8', ('x'))
    file_lat = dataset.createVariable('LAT', 'f8', ('x'))
    file_depth = dataset.createVariable('Depth', 'f8', ('x'))
    file_degree = dataset.createVariable('Degree Rot', 'f8', ('time'))
	
    file_Z = dataset.createVariable('Y', 'f8', ('z'))
    file_TIME = dataset.createVariable('TIME', 'f8', ('time'))
	
    VAR = dataset.createVariable(str(nameShort), 'f8', ('time','z','x'))
    VARFILT = dataset.createVariable('Filt' + str(nameShort), 'f8', ('time','z','x'))

    dataset.title = title
    dataset.author = 'Amelia Thelandersson'
    dataset.institution = 'Departamento de Oceanografía Física, Centro de Investigación Científica y de Educación Superior de Ensenada'
    #dataset.source = 'bitbucket.org/CanyonsUBC/BuildCanyon/Bathymetry/GenerateTankBathymetry_Inserts.ipynb'
    dataset.description = description
    dataset.timeStamp = time.ctime(time.time())
    file_X.standard_name = 'Crosshelf Distance'
    file_X.units = 'km'
    file_lon.standard_name = 'Longitude'
    file_lon.units = '°W'	
    file_lat.standard_name = 'Latitude'
    file_lat.units = '°N'
    file_depth.standard_name = 'Depth with distance from coast'
    file_lon.units = 'm'	
    file_degree.standard_name = 'Degree of rotation'
    file_lat.units = 'rad'
	
    file_Z.standard_name = 'Depth'
    file_Y.units = 'm'
    file_time.standard_name = 'Time'
    file_time.units = 'timedelta64[ns]'
    file_time.calendar = 'gregorian'	
    VAR.standard_name = str(nameLong)
    VAR.units = str(units)
	
    VARFILT.standard_name = 'Filtered' + str(nameLong)
    VARFILT.units = str(units)
    #VAR.positive = 'upward'

    file_X[:] = X[:]
    file_Z[:] = Z[:]
    file_TIME[:] = t[:]
    file_lon[:] = lon[:]
    file_lat[:] = lat[:]
    file_degree[:] = deg[:]
    file_depth[:] = dep[:]
    VAR[:] = var[:]
    VARFILT[:] = varfilt[:]

    dataset.close()

def recenter(vel,Z,LON,LAT,lon,lat):
	Recent=np.zeros((len(Z),len(lat)))
	
	for d in range(len(Z)):
		interp=sciint.RegularGridInterpolator((LAT,LON),vel.values[d])
		Recent[d]=interp((lat,lon))
	
	return Recent

def CrossectExctraction(i,dsw,dsn,filt,detrend,var,corrind,coast,nr):
	
	Z=dsw[0].Z.values
	LAT = dsw[0].YC.values
	LON = dsw[0].XC.values - 360
	
	if coast == 'smooth':
		day=9
		time12=dsw[0].time.values.astype(int)
		time23=dsw[1].time.values.astype(int)
		time34=dsw[2].time.values.astype(int)
		time45=dsw[3].time.values.astype(int)
		time56=dsw[4].time.values.astype(int)
		time67=dsw[5].time.values.astype(int)
		time78=dsw[6].time.values.astype(int)
		time89=dsw[7].time.values.astype(int)
		time910=dsw[8].time.values.astype(int)
		
		Time=np.concatenate((time12, time23, time34, time45, time56,time67, time78,time89, time910), axis=0)#, time910), axis=0)
		times=Time*1e-9
	else:
		day=8
		time23=dsw[0].time.values.astype(int)
		time34=dsw[1].time.values.astype(int)
		time45=dsw[2].time.values.astype(int)
		time56=dsw[3].time.values.astype(int)
		time67=dsw[4].time.values.astype(int)
		time78=dsw[5].time.values.astype(int)
		time89=dsw[6].time.values.astype(int)
		time910=dsw[7].time.values.astype(int)
		
		Time=np.concatenate((time23, time34, time45, time56,time67, time78,time89, time910), axis=0)#
		#inds=np.append(np.arange(0,433,2),np.arange(433,792,1))
		times=Time*1e-9
	
	
	matfile=loadmat('/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/BT_P_res' + str(nr) + '.mat')
	x,dep,lon,lat,deg=matfile['dist'][0][i][0],matfile['d'][0][i][0],matfile['lon'][0][i][0],matfile['lat'][0][i][0],matfile['degree'][0][i]

	
	
	VALMITpre=np.zeros((len(times), len(Z),len(lon)))
	
	if np.logical_or(var=='ashore',var=='cshore'):
		print(var)
		pathETA='/home/athelandersson/NETCDFs/' + str(coast) + '/ETANAC.nc'
		ds= xr.open_dataset(pathETA)
			
		for tt in np.arange(0,day,1):
	    		
			VALMIT=np.zeros((len(times),len(Z),len(lon[i])))
				    
			for t in np.arange(0,len(times),1):					

				Ub=recenter(dsw[tt].UVEL[t],Z,LON,LAT,lon,lat)
				Vb=recenter(dsw[tt].VVEL[t],Z,LON,LAT,lon,lat)
				Un=recenter(dsn[tt].UVEL[t],Z,LON,LAT,lon,lat)
				Vn=recenter(dsn[tt].VVEL[t],Z,LON,LAT,lon,lat)
				
				VALb=(sin(deg) * Ub +  Vb * cos(deg))
				VALn=(sin(deg) * Un +  Vb * cos(deg))
				VALmit=VALb-VALn
				VALMIT[t,:,:]=VALmit
			
			VALMITpre[len(dsw[tt-1].time)*tt:len(VALMIT[:,1,1])*(tt+1),:,:]=VALMIT
			print('Day '+str(tt+2))
			
	else: 
		print(var)
		for tt in np.arange(0,day,1):
	    		
			VALMIT=np.zeros((len(times),len(Z),len(lon)))
				    
			for t in np.arange(0,times,1):
				exec (f'VALb=recenter(dsw[tt].{VAR}[t],Z,LON,LAT,lon,lat)')
				exec (f'VALn=recenter(dsn[tt].{VAR}[t],Z,LON,LAT,lon,lat)')
				VALmit=VALb-VALn
				VALMIT[t,:,:]=VALmit
			
			VALMITpre[len(dsw[tt-1].time)*tt:len(VALMIT[:,1,1])*(tt+1),:,:]=VALMIT
			print('Day '+str(tt+2))
			
	return(VALfilt,VALMITpre,x,dep,lon,lat,deg,Z,times)
		


def ExtractAndFiltCrossectNEW(i,dsw,dsn,filt,detrend,var,corrind,coast):

	Z=dsw[0].Zl.values
	hFacC = dsw[0].hFacC
	hfac = np.ma.masked_values(hFacC, 0)
	mask = np.ma.getmask(hfac)
	
	if coast == 'smooth':
		day=9
		time12=dsw[0].time.values.astype(int)
		time23=dsw[1].time.values.astype(int)
		time34=dsw[2].time.values.astype(int)
		time45=dsw[3].time.values.astype(int)
		time56=dsw[4].time.values.astype(int)
		time67=dsw[5].time.values.astype(int)
		time78=dsw[6].time.values.astype(int)
		time89=dsw[7].time.values.astype(int)
		time910=dsw[8].time.values.astype(int)
		
		Time=np.concatenate((time12, time23, time34, time45, time56,time67, time78,time89, time910), axis=0)#, time910), axis=0)
		times=Time*1e-9
	else:
		day=8
		time23=dsw[0].time.values.astype(int)
		time34=dsw[1].time.values.astype(int)
		time45=dsw[2].time.values.astype(int)
		time56=dsw[3].time.values.astype(int)
		time67=dsw[4].time.values.astype(int)
		time78=dsw[5].time.values.astype(int)
		time89=dsw[6].time.values.astype(int)
		time910=dsw[7].time.values.astype(int)
		
		Time=np.concatenate((time23, time34, time45, time56,time67, time78,time89, time910), axis=0)#
		inds=np.append(np.arange(0,433,2),np.arange(433,792,1))
		times=Time[inds]*1e-9
	
	
	matfile=loadmat('/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/BT_P.mat')
	x,dep,indXlon,indYlat=matfile['dist'],matfile['d'],matfile['indexXlon'],matfile['indexYlat']

	maskin=mask[:,indYlat[i],indXlon[i]]
	if coast == 'original': 
		VALMITpre=np.zeros((72*5+144*3,np.size(maskin[:,0]),np.size(maskin[0,:])))
	elif coast == 'smooth':
		VALMITpre=np.zeros((len(Time),np.size(maskin[:,0]),np.size(maskin[0,:])))
	if var=='UVEL':
		print('u')
		for tt in np.arange(0,9,1):
	    		
			VALMIT=np.zeros((len(dsw[tt].UVEL[:,1,1,1]),len(Z),len(indXlon[i])))
				    
			for t in np.arange(0,len(dsw[tt].UVEL[:,1,1,1]),1):
				VALb=dsw[tt].UVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]]
				VALn=dsn[tt].UVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]]
				VALmit=VALb-VALn
				VALMIT[t,:,:]=VALmit
			
			VALMITpre[len(dsw[tt-1].UVEL[:,1,1,1])*tt:len(VALMIT[:,1,1])*(tt+1),:,:]=VALMIT
			print('Day '+str(tt+2))
			
	elif var=='VVEL':
		print('v')
		for tt in np.arange(0,9,1):
			VALMIT=np.zeros((len(dsw[tt].VVEL[:,1,1,1]),len(Z),len(indXlon[i])))
				    
			for t in np.arange(0,len(dsw[tt].VVEL[:,1,1,1]),1):
				VALb=dsw[tt].VVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]]
				VALn=dsn[tt].VVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]]
				VALmit=VALb-VALn
				VALMIT[t,:,:]=VALmit
			
			VALMITpre[len(dsw[tt-1].VVEL[:,1,1,1])*tt:len(VALMIT[:,1,1])*(tt+1),:,:]=VALMIT
			print('Day '+str(tt+2))
			
	elif var=='WVEL':
		print('w')
		for tt in np.arange(0,9,1):
			VALMIT=np.zeros((len(dsw[tt].WVEL[:,1,1,1]),len(Z),len(indXlon[i])))
				    
			for t in np.arange(0,len(dsw[tt].WVEL[:,1,1,1]),1):
				VALb=dsw[tt].WVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]]
				VALn=dsn[tt].WVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]]
				VALmit=VALb-VALn
				VALMIT[t,:,:]=VALmit
			

			VALMITpre[len(dsw[tt-1].WVEL[:,1,1,1])*tt:len(VALMIT[:,1,1])*(tt+1),:,:]=VALMIT
			print('Day '+str(tt+2))

	
	elif var=='PHIHYD':
		print('phi')
		for tt in np.arange(0,9,1):
			VALMIT=np.zeros((len(dsw[tt].PHIHYD[:,1,1,1]),len(Z),len(indXlon[i])))
			for t in np.arange(0,len(dsw[tt].PHIHYD[:,1,1,1]),1):
				VALb=dsw[tt].PHIHYD[t,:,:,:].values[:,indYlat[i],indXlon[i]]
				VALn=dsn[tt].PHIHYD[t,:,:,:].values[:,indYlat[i],indXlon[i]]
				VALmit=VALb-VALn
				VALMIT[t,:,:]=VALmit
		

			VALMITpre[len(dsw[tt-1].PHIHYD[:,1,1,1])*tt:len(VALMIT[:,1,1])*(tt+1),:,:]=VALMIT
			print('Day '+str(tt+2))
	if var=='cshore':
		print('crosshore')
		pathETA='/home/athelandersson/NETCDFs/' + str(coast) + '/ETANAC.nc'
		ds= xr.open_dataset(pathETA)
		
		lon_ac=ds.lonAC.values
		lat_ac=ds.latAC.values
		distAC=ds.dist.values
	
		
		LAT = dsw[0].YC
		LON = dsw[0].XC - 360
		
		ind=i

		lon1=LON[lon_ac[ind]]
		lat1=LAT[lat_ac[ind]]
		lon2=LON[-1]
		lat2=LAT[lat_ac[ind]]
		a=haversine(lon1, lat1, lon2, lat2)
		
		lon3=LON[-1]
		lat3=LAT[0]
		
		b=haversine(lon2, lat2, lon3, lat3)
		
		deg=-atan(a/b)

			
		for tt in np.arange(0,day,1):
	    		
			VALMIT=np.zeros((len(dsw[tt].UVEL[:,1,1,1]),len(Z),len(indXlon[i])))
				    
			for t in np.arange(0,len(dsw[tt].UVEL[:,1,1,1]),1):
				VALb=(cos(deg) * dsw[tt].UVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]] - dsw[tt].VVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]]*sin(deg))
				VALn=(cos(deg) * dsn[tt].UVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]] - dsn[tt].VVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]]*sin(deg))
				VALmit=VALb-VALn
				VALMIT[t,:,:]=VALmit
			
			VALMITpre[len(dsw[tt-1].UVEL[:,1,1,1])*tt:len(VALMIT[:,1,1])*(tt+1),:,:]=VALMIT
			print('Day '+str(tt+2))
			
	if var=='ashore':
		print('alongshore')
		pathETA='/home/athelandersson/NETCDFs/' + str(coast) + '/ETANAC.nc'
		ds= xr.open_dataset(pathETA)
		
		lon_ac=ds.lonAC.values
		lat_ac=ds.latAC.values
		distAC=ds.dist.values
	
		
		LAT = dsw[0].YC
		LON = dsw[0].XC - 360
		
		ind=i

		lon1=LON[lon_ac[ind]]
		lat1=LAT[lat_ac[ind]]
		lon2=LON[-1]
		lat2=LAT[lat_ac[ind]]
		a=haversine(lon1, lat1, lon2, lat2)
		
		lon3=LON[-1]
		lat3=LAT[0]
		
		b=haversine(lon2, lat2, lon3, lat3)
		
		deg=-atan(a/b)
		R1=LON*cos(deg)-LAT*sin(deg)
		R2=LON*sin(deg)+LAT*cos(deg)
			
		for tt in np.arange(0,day,1):
	    		
			VALMIT=np.zeros((len(dsw[tt].UVEL[:,1,1,1]),len(Z),len(indXlon[i])))
				    
			for t in np.arange(0,len(dsw[tt].UVEL[:,1,1,1]),1):
				VALb=(sin(deg) * dsw[tt].UVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]] +  dsw[tt].VVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]] * cos(deg))
				VALn=(sin(deg) * dsn[tt].UVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]] +  dsn[tt].VVEL[t,:,:,:].values[:,indYlat[i],indXlon[i]] * cos(deg))
				VALmit=VALb-VALn
				VALMIT[t,:,:]=VALmit
			
			VALMITpre[len(dsw[tt-1].UVEL[:,1,1,1])*tt:len(VALMIT[:,1,1])*(tt+1),:,:]=VALMIT
			print('Day '+str(tt+2))
	
	if coast == 'original': 
		VALMITpre=VALMITpre[inds]
		
	VALfilt=np.zeros(np.shape(VALMITpre))
	
	if filt==1:
		fs=1/1200
		fs2=0
		
		print('Filtering begins')
		for d in np.arange(np.size(VALMITpre,2)):
	    		VALdif,VALfiltout,VALfiltAll,inds = FiltDetrend(VALMITpre[:,:,d],filt,detrend,fs,fs2)
	    		VALfilt[:,:,d]=VALfiltAll
	else:
		print('No filtering')
		VALfilt = 0
		VALfilttwe=0
		inds=0	

	return(VALfilt,VALMITpre,x[i],Z,times)
      
def get_Brink(file_fig):#,file_h): #, file_ratio):
    # Brink mode
    file = sio.loadmat(file_fig)
    z, xpl, xxx, zzz, xgr, zgr = file['z'][0,:], file['xpl'][0,:], file['xxx'][0,:], file['zzz'][0,:], file['xgr'], file['zgr']
    k, omega,epe,eke = file['wavenumber'][0][0], file['frequency'][0][0], file['epe'][0][0],file['eke'][0][0]

    # (u is cross-shore and v is alongshore in Brink.)
    p0, u0, v0, w0, r0 = file['p'], file['u'],file['v'], file['wvel'], file['rho']

    scale=0.2
    w = w0.transpose() * 0.01 * scale # cms-1 to ms-1 and normalization (?)
    u = u0.transpose() * 0.01 * scale # cms-1 to ms-1 and normalization 
    v = v0.transpose() * 0.01 * scale # cms-1 to ms-1 and normalization 
    r = r0.transpose() * 1.0 * scale # mg/cm³ to kg/m³ and normalization
    p = p0.transpose() * 0.1 * scale # dyn/cm² to 0.1 Pa (or kg m-1 s-2) and normalization
    

    
    return(u,v,w,r,p,z,k,omega,xpl, xxx, zzz,zgr.transpose(),xgr.transpose(), epe, eke)
    
def interpolate(VALmit,valinBrink,dist,xpl,Z,z,zgr,xgr,mask):
    
    #finding the coordinates for all points in the brink variable
    zpib,xpib=np.where(np.array(valinBrink[1])<1000)
    zpim,xpim=np.where(VALmit!=0)
    
    valbrink1d=[]
    #picking out the values using the coordinates and creating a 1D vector
    for i in np.arange(0,len(valinBrink),1):
         valbrink1dPRE=valinBrink[i][zpib,xpib]
         valbrink1d.append(valbrink1dPRE)
    
    valmitin=VALmit[zpim,xpim]
    
    # go back to 2D
    # valmit2d = np.zeros(np.shape(VALmit))
    #for i in np.arange(len(valmitout)):
    #    valmit2d[zpi[i],xpi[i]]=valmitout[i]
    

    #Creating a grid of Z which does not include the values in the coast, 
    #similar to the grid from Brink
    grid_Z=np.zeros(np.shape(VALmit))
    for i in np.arange(0,len(mask[0,:]),1):
	#finding the indices of the coast at one x-location
        coastz=np.where(VALmit[:,i]!=0) 
        if len(Z[coastz])>0:
        #creating a vector of the same size as how many
            places=np.arange(0,76,76/len(Z[coastz]))   
            #depth values we have that aren't the coast
            #Sometimes the vector became too long 
            if len(places)> len(Z[coastz]):           
                places=places[:-1] 
            #making sure the last value in the vector of indices 
            #is 71 so we interpolate over the same start and finish vaues
            places[-1]=76       
            #Interpolating the values that arent the coast
            #to agree with the amount of depth values the 
            # the model outputs. Using the range of indices
            #given from the vector created to show indices
            # that the depth values not on the coast correspond to
            # and then filling the rest with values between these.
             
            grid_Z[:,i]=np.interp(np.arange(0,76,1),places,Z[coastz])
            
    
    grid_Z[grid_Z==0]=-2                                                              
    #Creating a grid of the x-values, which are the same for every Z
    grid_X,grid_no = np.meshgrid(dist,Z)
    
    #find the depth and distance for each of the indices in the Brink code
    points=[zgr[zpib,xpib],xgr[zpib,xpib]] 
    
    #transposing this to so the values are the columns and not the rows for the interpolation
    points=np.asarray(points).transpose()
    
    valbrinkint=[]
    for i in np.arange(0,len(valbrink1d),1):
        valbrinkout  = sciint.griddata(points, np.array(valbrink1d[i]), (grid_Z, grid_X), method='linear')
        valbrinkint.append(valbrinkout)
    
    pointsMIT=[grid_no[zpim,xpim],grid_X[zpim,xpim]] 
    pointsMIT=np.asarray(pointsMIT).transpose() 
    
    valmitint = sciint.griddata(pointsMIT, valmitin, (grid_Z, grid_X), method='linear')
   
    zpi,xpi=np.where(~np.isnan(valbrinkout))
    
    valbrinkR=[]
    for i in np.arange(0,len(valbrink1d),1):
        xOUT= valbrinkint[i][zpi,xpi]
        valbrinkR.append(xOUT)

    valmit1d= valmitint[zpi,xpi]
    
    return valmit1d,valbrinkR,xpi,zpi,valbrinkint,grid_X,grid_Z,xpi,valmitint

def lin_reg(VALmit,valinBrink,dist,xpl,Z,z,zgr,xgr,maskin):

    Ypre,valbrinkR,xpi,zpi,valbrinkint,grid_X,grid_Z,xpi,valmitint=interpolate(VALmit,valinBrink,dist,xpl,
                                                                              Z,z,zgr,xgr,maskin)
    Y=Ypre-np.mean(Ypre)
    
    
    for l in np.arange(0,len(valbrinkR),1):
        if l==0:
           Xpre=np.vstack((np.ones(len(valbrinkR[l])),valbrinkR[l]))
        else:
           Xpre=np.vstack((Xpre,valbrinkR[l]))
		
    X=Xpre.T
    beta_hat = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)
    
    yhat=X.dot(beta_hat) +np.mean(Ypre)
    valout = np.zeros(np.shape(valbrinkint[0]))
    varbrink = np.zeros(np.shape(valbrinkint))

    
    for i in np.arange(len(Y)):
        valout[zpi[i],xpi[i]]=yhat[i]
        for k in np.arange(0,len(varbrink),1):
        	varbrink[k][zpi[i],xpi[i]]=valbrinkR[k][i]

    
    return beta_hat,yhat,X,valout,varbrink,grid_X,grid_Z,Y,xpi,Ypre,valmitint

def fitmodes(dsw,dsn,valinBrink,xpl,Z,z,indXlon,indYlat,dist,zgr,xgr,ds,filt,time,coast):
    
    hFacC = dsw[0].hFacC
    hfac = np.ma.masked_values(hFacC, 0)
    mask = np.ma.getmask(hfac)


    maskin=mask[:,indYlat,indXlon]
    VALMIT=np.zeros((len(time),len(Z),len(dist)))
    VALfit=np.zeros((len(time),len(Z),len(dist)))
    betas=np.zeros((len(time),len(valinBrink)+1))
    fit=np.zeros((len(time)))
    RMSE=np.zeros((len(time)))
    for t in np.arange(0,len(time),1):
            print(str(time[t]))
            VALmit=ds.VAL[t,:,:].values
            beta_hat,yhat,xbeta,valout,varbrink,grid_X,grid_Z,Y,xpi,Ypre,valmitint=lin_reg(VALmit,valinBrink,dist,xpl,Z,z,zgr,xgr,maskin)
            VALfit[t,:,:]=valout
            betas[t,:]=beta_hat
            VALMIT[t,:,:]=valmitint
            
            mse=np.mean((Ypre - yhat) ** 2)
            rmse=sqrt(mse)
            RMSE[t]=rmse
            
            yi=VALmit
            sst=np.sum((yi-np.mean(yi))**2)
            yihat=valout
            ssr=np.sum((yihat-np.mean(yihat))**2)
            fit[t]=100*(ssr/sst)
        
    
    return VALfit,betas,xbeta,yhat,dist,VALMIT,varbrink,grid_X,grid_Z,fit,Y,xpi,Ypre,RMSE,maskin,valmitint

def linearregressionSave(filt,var,coast):
	if coast == 'smooth':
		startday=1
	elif coast == 'original':
		startday=2
	
	hej=[35] #,54,79,120,154,194,219]  
	corrinds=[30.49] #,30.77,31.13,31.69,32.11,32.65,33.02] 
	
	
	for ik in np.arange(0,len(hej),1):
		print(str(corrinds[ik]))
		
		u=[]
		v=[]
		w=[]
		r=[]
		p=[]
		k=[]
		omega=[]
		epe=[]
		eke=[]
		
		
		for l in np.arange(0,25,1):
			if exists('/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/CrossectsPerp/dataSVB'+ str(corrinds[ik]) +'mode' + str(l) + '.mat') == True:
				print('Mode ' + str(l))
				uo,vo,wo,ro,po,z,ko,omegao, xpl, xxx, zzz, zgr, xgr, epeo, ekeo = get_Brink('/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/CrossectsPerp/dataSVB'+ str(corrinds[ik]) +'mode' + str(l) + '.mat')	
				u.append(uo.imag) 
				v.append(vo)
				w.append(wo.imag)
				r.append(ro)
				p.append(po)
				k=np.append(k,ko)
				omega=np.append(omega,omegao)
				epe=np.append(epe,epeo)
				eke=np.append(eke,ekeo)
			

			dirn='/home/athelandersson/NETCDFs/' + str(coast) + '_NO/'
			dirw='/home/athelandersson/NETCDFs/' + str(coast) + '/'

			if var == 'PHIHYD':
				dsw,dsn=loadNetCDFs(dirw,dirn,'phiHyd',startday)
			else:
				dsw,dsn=loadNetCDFs(dirw,dirn,'dynVars',startday)
			
		
		ds=xr.open_dataset('/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/Locations/' + str(var) + str(corrinds[ik]) + str(filt)+ '.nc')
		
		matfile=loadmat('/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/BT_P.mat')
		dist,indXlon,indYlat=matfile['dist'][ik],matfile['indexXlon'][ik],matfile['indexYlat'][ik]
		
		Z=dsw[0].Zl.values
		TIME=ds.time.values
		
		if var == 'PHIHYD':
			valinBrink=p
		elif var == 'ashore':
			valinBrink=v
		elif var == 'cshore':
			valinBrink=u
		elif var == 'WVEL':
			valinBrink=w
		VALfit,betas,xbeta,yhat,dist,VALmit,varbrink,grid_X,grid_Z,fit,Y,xpi,Ypre,RMSE,mask,valmitint=fitmodes(dsw,dsn,valinBrink,xpl,Z,z,indXlon,indYlat,dist,zgr,xgr,ds,filt,TIME,coast)
		
		FILENAME='/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/' + str(var) + '/LinReg' + str(corrinds[ik]) + str(filt)+ '.nc'
		ds = xr.Dataset({'valfit': (("time","z","x"), VALfit),
				 'valmit': (("time","z","x"), VALmit),
				 'varbrink': (("nrM","z","x"), varbrink),
				 'fit':(("time"), fit),
				 'gridX':(("z","x"), grid_X),
				 'gridZ':(("z","x"), grid_Z),
				 'betas':(("time","nrB"), betas),
				 'rmse':(("time"), RMSE)
				    },
				coords ={
				    "x" : dist,
				    "z" : Z,
				    "time": TIME,
				    "nrM": np.arange(0,len(varbrink),1),
				    "nrB": np.arange(0,len(betas[0,:]),1)
				},
				)
		ds.to_netcdf(FILENAME)

def LoadLinReg(ds):
    VALFIT=ds.valfit.values
    VALMIT=ds.valmit.values
    VALBRINK=ds.varbrink.values
    FIT=ds.fit.values
    grid_Z=ds.gridZ.values
    grid_X=ds.gridX.values
    BETA=ds.betas.values
    RMSE=ds.rmse.values
    dist=ds.x.values
    Z=ds.z.values
    TIME=ds.time.values
    
    return VALFIT, VALMIT, VALBRINK, FIT, grid_Z, grid_X, BETA, RMSE, dist, Z, TIME
		
def loadalllats(filt,var):
	hej=[58, 85, 205, 227] 
	corrinds=[30.84,31.2,32.68,32.98]

	VALfit=[]
	VALmit=[]
	vals=[]
	fit=[]
	grid_Z=[]
	grid_X=[]
	betas=[]
	RMSE=[]
	dist=[]
	Z=[]
	zgr=[]
	xgr=[]
	
	for ik in np.arange(0,len(hej),1):
		FILENAME='LinReg' + str(var) + str(corrinds[ik]) + str(filt)+ '.nc'
		ds=xr.open_dataset('Locations/'+ str(FILENAME)) 
		
		VALfito, VALmito, valso, fito, grid_Zo, grid_Xo, betaso, RMSEo, disto, Zo, TIME =LoadLinReg(ds)
		
		uo,vo,wo,ro,po,z,ko,omegao, xpl, xxx, zzz, zgro, xgro, epeo, ekeo = get_Brink('PerpendicularCrossections/dataSVB'+ str(corrinds[ik]) +'mode1.mat')
		
		VALfit.append(VALfito)
		VALmit.append(VALmito)
		vals.append(valso)
		fit.append(fito)
		grid_Z.append(grid_Zo)
		grid_X.append(grid_Xo)
		betas.append(betaso)
		RMSE.append(RMSEo)
		dist.append(disto)
		Z.append(Zo)
		zgr.append(zgro)
		xgr.append(xgro)
		
	return(VALfit,VALmit,vals,fit,grid_Z,grid_X,betas,RMSE,dist,Z,TIME,xgr,zgr)

def openBrink(loc):
	u=[]
	v=[]
	w=[]
	r=[]
	p=[]
	k=[]
	omega=[]
	epe=[]
	eke=[]


	for l in np.arange(0,10,1):
		if exists('/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/CrossectsPerp/dataSVB'+ str(loc) +'mode' + str(l) + '.mat') == True:
			uo,vo,wo,ro,po,z,ko,omegao, xpl, xxx, zzz, zgr, xgr, epeo, ekeo = get_Brink('/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/CrossectsPerp/dataSVB'+ str(loc) +'mode' + str(l) + '.mat')	
		else: 
			uo=0
			vo=0
			wo=0
			ro=0
			po=0
			ko=0
			omegao=0
			epeo=0
			ekeo=0
			
		u.append(uo) 
		v.append(vo)
		w.append(wo)
		r.append(ro)
		p.append(po)
		k=np.append(k,ko)
		omega=np.append(omega,omegao)
		epe=np.append(epe,epeo)
		eke=np.append(eke,ekeo)
	
	
	return(u,v,w,r,p,k,omega,epe,eke,xgr,zgr)

## FOR PLOTTING ----------------------------------------------------------------------------------

def get_snapshot_at_level(t,dep,dsw,dsn,var):
    ind=0
    if t>=72 and t <(72*2):
        ind=1
        t=t-72
    elif t>=(72*2) and t<(72*3):
        ind=2
        t=t-(72*2)
    elif t>=(72*3) and t<(72*4):
        ind=3
        t=t-(72*3)
    elif t>=(72*4) and t<(72*5):
        ind=4
        t=t-(72*4)
    elif t>=(72*5) and t<(72*6):
        ind=5
        t=t-(72*5)
    elif t>=(72*6) and t<(72*7):
        ind=6
        t=t-(72*6)
    elif t>=(72*7) and t<(72*8):
        ind=7
        t=t-(72*7)
    elif t>=(72*8) and t<(72*9):
        ind=8
        t=t-(72*8)
    if var=='WVEL':
        Ww=dsw[ind].WVEL[t,dep,:,:].values
        Wn=dsn[ind].WVEL[t,dep,:,:].values
        W = Ww-Wn
    elif var=='ETAn':
        W=dsn[ind].ETAN[t,:,:].values
    else:
        W=dsw[ind].ETAN[t,:,:].values
    return(W)

def plotMap(ax,LON,LAT,depth,mask,fig,nr):
	ax.set_facecolor('tan')

	pc = ax.contourf(LON,LAT,np.ma.masked_array(depth, mask=mask[0,:,:]),50,vmin=0, vmax=5000, cmap=cmocean.cm.deep)


	cn = ax.contour(LON,LAT,depth, colors=['0.2','0.4'],levels=[0,500])
	divider = make_axes_locatable(ax)
	axdiv = divider.new_vertical(size = '5%', pad = 0.5)
	fig.add_axes(axdiv)
	cbar_ax = plt.colorbar(pc,cax=axdiv,orientation='horizontal',ticks=np.arange(0,np.max(depth),1000))
	cbar_ax.ax.xaxis.set_label_position("top")
	cbar_ax.set_label('Depth [m]')


	ax.set_xlabel('Lon [°]')
	ax.set_ylabel('Lat [°]')
	ax.set_xlim(238-360, 246-360)
	ax.set_ylim(27,35.3)
	ax.set_aspect(1)
	ax.text(-0.1, 1.2, nr, fontweight='bold', color='k', transform=ax.transAxes)
        
def plotpointsAC(LON,LAT,lon_inds,lat_inds,var):
	params = {'font.size': 18,
          'figure.figsize': (30, 30),
         'font.family':'serif'}
	pl.rcParams.update(params)
	
	dirw = '/home/athelandersson/NETCDFs/smooth/'
	dirn = '/home/athelandersson/NETCDFs/smooth_NO/'

	dsw,dsn=loadNetCDFs(dirw,dirn,'rhoAnoma')
	depth = dsw[0].Depth
	
	fig, ax = plt.subplots(1,1,figsize=(10,9))
	ax.set_facecolor('tan')
	pc = ax.pcolormesh(LON,LAT,depth, cmap=cmocean.cm.deep)#, extend='max')
	cb = plt.colorbar(pc, extend='max',label='depth / m')
	cn = ax.contour(LON,LAT,depth, colors=['0.3','0.6'], 
                	levels=[0,500])

	for ii,jj in zip(lon_inds,lat_inds):

    		ax.plot(LON[ii-1],LAT[jj-1],'o', 
            		markersize=4, color='r')

	# To show it begins and ends where it should
	ax.plot(LON[lon_inds[0]],LAT[lat_inds[0]],'o', 
           	markersize=10, color='orange') 
	ax.plot(LON[lon_inds[-1]],LAT[lat_inds[-1]],'o', 
           	markersize=10, color='blue') 



	cb.set_label('Depth [m]')
	ax.set_xlabel('Lon [°]')
	ax.set_ylabel('Lat [°]')
	ax.set_xlim(238-360, 246-360)
	ax.set_ylim(27,35.3)
	ax.set_aspect(1)
	plt.savefig('/home/athelandersson/CTW-analysis/Figures/PAC' + str(var) + '.png')
	
	
def plotsnapshot(ax,VAL,dep,LON,LAT,vmin,vmax,depth,label,nr,TIME,mask):
	xlab='Longitude [°]'
	ylab='Latitude [°]'
	
	ax.set_facecolor('tan')
	cax = ax.pcolormesh(LON,LAT,np.ma.masked_array(VAL, mask=mask[dep,:,:]),cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax)
	ax.contour(LON,LAT,depth, colors=['0.2','0.6'],levels=[0,500])
	
	ax.set(xlabel=xlab,ylabel=ylab)

	ax.text(0.4,0.87, f'Day {TIME} \n' + label, transform=ax.transAxes,horizontalalignment='left')
	ax.text(-0.1,1.05, nr, transform=ax.transAxes)
	         
	ax.set_xlim(-122,-114)
	ax.set_ylim(27,35.3)
	ax.set_aspect(1)
	return(cax)
	
def plot_HOVMOLLER(ax,LON,TIME,VAL,title,ctitle,vmin,vmax,fig,lat,lon,lab,cbarall,nr):
    
   
    xlab='Time [days]'
    if lab==1:
        ylab='Distance [km]' #'Depth [m]'
    else:
        ylab=''

    ax.set(xlabel=xlab, ylabel=ylab)
    ax.set_xticks([2880, 4320, 5760, 7200, 8640, 10080, 11520, 12960, 14400])
    ax.set_xticklabels([2, 3, 4, 5, 6, 7, 8, 9, 10])
    ax.set_title(title)

    if ctitle=='SSH [mm]':
        cax = ax.pcolormesh(TIME,LON,np.transpose(VAL),cmap=cmocean.cm.curl,vmin=vmin,vmax=vmax) 
    else:    
        cax = ax.pcolormesh(TIME,LON,np.transpose(VAL),cmap=cmocean.cm.balance,vmin=vmin,vmax=vmax) 
   
    if cbarall==1:
    ##FOR THE SAME COLORBAR FOR ALL OF THE PLOTS
        cbar_ax = fig.add_axes([1, 0.15, 0.03, 0.7])
        fig.colorbar(cax, cax=cbar_ax)
        cbar_ax.set_ylabel(ctitle)
    else:
    ##FOR A COLORBAR FOR EACH PLOT
        divider = make_axes_locatable(ax)
        axdiv = divider.new_vertical(size = '5%', pad = 0.5)
        fig.add_axes(axdiv)
        cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
        cbar_ax.ax.xaxis.set_label_position("top")
        cbar_ax.set_label(ctitle)
    #ax.set_aspect(1./ax.get_data_ratio())
    ax.text(-0.15, 1.2, nr, fontweight='bold', color='k', 
        transform=ax.transAxes)

def plot_batylines(lat_ac,dist,LAT,ax,lab,hej):

	colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999']
	matfile=loadmat('BT_P.mat')
	x,dep,indXlon,indYlat=matfile['dist'],matfile['d'],matfile['indexXlon'],matfile['indexYlat']


	for i in np.arange(0,len(hej),1):
    		ax.plot(x[i],-dep[i], label=f'{dist[hej[i]]:.0f} km \n {LAT[lat_ac[hej[i]]]:.2f} °N', color=colors[i])
    		ax.legend()
    		ax.set_xlabel('Distance [km]')
    		ax.set_ylabel('Depth [m]')
    		ax.text(-0.08, 1.02, lab, fontweight='bold', color='k',transform=ax.transAxes)
    		
    		
def plotbrink(ax,grid_X,grid_Z,levelsb,xlab,ylab,modenr,nr,varbrink,modes,lat,t):
    ax.set_facecolor('tan')
    cax=ax.contourf(grid_X,grid_Z,varbrink ,cmap=cmocean.cm.delta,levels=levelsb)
    ax.contour(grid_X,grid_Z,varbrink , levels=[0], linewidths=2, 
                linestyles='-', colors='k', zorder=2)
   
    if nr<=-1:
        ax.set_title(f'MITgcm cross-section \n Day {t}', fontdict={'fontsize': 15})
        if nr==-1:    
            ax.set(ylabel=ylab)
        else: 
            ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False) 
        ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)
    else:
        ax.set_title(f'Mode {modenr}')
        if nr>=modes-3:
            ax.set(xlabel=xlab)
        else:
            ax.tick_params(axis='x',which='both', bottom=True, top=False, labelbottom=False)
        if (nr % 2) != 0:
            ax.tick_params(axis='y',which='both', left=True, right=False, labelleft=False) # labels along the bottom edge are off
        else: 
            ax.set(ylabel=ylab)
            
    return cax

	    		

