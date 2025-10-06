



import numpy as np
from math import cos, sin, asin
import cmocean
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

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


def closest(lat,indlat,lon,indlon):
    loclat=[]
    loclon=[]
    for i in range(len(indlat)):
        loclat.append(np.where(lat>=indlat[i])[0][0])
        loclon.append(np.where(lon<=indlon[i])[0][0])
    return loclat,loclon


def FFRQ(var,timemin,dist):

    times = (timemin) #in s

    t0 = 172800 # start is day 2 

    dt = 1200 # 20 min, sample frequency 
    fs=1/dt
        
    nx = len(dist)
    nt = int(720/2+1)
        
    psd = np.zeros((nx,nt))
    phase = np.zeros((nx,nt))
    
    freq =  np.zeros((nx,nt))
    nr=720

    if len(var)>576:
        start=72
    else:
        start=144

    for ii in np.arange(nx):
        arrin=np.zeros(nr)
        N = len(var[:,ii])
        
        arrin[start:]=var[:,ii]
        
        xdft = np.fft.rfft(arrin,n=nr) #cycles/second
        xdft = xdft[0:int(nr/2+1)]
        psdx = (1/(fs*N)) * np.abs(xdft)**2
        psdx[1:-1] = 2*psdx[1:-1]
        FFTfreq = np.fft.rfftfreq(nr, d=dt)
    
        psd[ii,:] = psdx
        freq[ii,:] =  FFTfreq

    return psd, freq

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
        cbar_ax = fig.add_axes([1, 0.15, 0.03, 0.6])
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
    ax.text(-0.15, 1.05, nr, fontweight='bold', color='k', 
        transform=ax.transAxes)


