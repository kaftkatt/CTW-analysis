#!/usr/bin/env python
# coding: utf-8

# # Calculate the power spectral density and visualise it

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import cmocean
import xarray as xr
import pylab as pl
from math import radians, cos, sin, asin, sqrt
import SVBfunc
from scipy.io import loadmat

pathVEL='WVELAC.nc'
dsVEL= xr.open_dataset(pathVEL)

Wfilt=dsVEL.ValfiltAll.values
Wdif=dsVEL.Valorig.values
dist=dsVEL.dist.values
TIME=dsVEL.time2.values
lat_ac=dsVEL.latAC.values


# In[5]:


pathETA='ETAAC.nc'
dsETA= xr.open_dataset(pathETA)

etafilt=dsETA.ValfiltAll.values
etadif=dsETA.Valorig.values
disteta=dsETA.dist.values
TIMEeta=dsETA.time2.values
lat_aceta=dsETA.latAC.values


# #### Loading pressure data just for the domain characteristics

# In[6]:


varname='PHIHYD'
i=0
pathn='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL/01b_noSVB_febTS/'+ str(varname)+'noSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
pathw='/media/amelia/Trillian/SVB/exp06_512x612x100_ORL_SVB/01b_SVB_febTS/'+ str(varname)+'withSVB'+ str(2+i)+'_'+ str(3+i) +'.nc'
        
dsw  = xr.open_dataset(pathw)
dsn = xr.open_dataset(pathn)


# In[7]:


LAT=dsw.YC
LON=dsw.XC-360
Z = dsw.Z
hFacC = dsw.hFacC

hfa = np.ma.masked_values(hFacC[:,:,:], 0)
mask = np.ma.getmask(hfa)
        
depth=dsw.Depth
depthno=dsn.Depth

matfile=loadmat('BT_PERP.mat')
indXlon,indYlat=matfile['indexXlon'],matfile['indexYlat']

# In[8]:
inds=np.append(np.arange(0,433,2),np.arange(433,792,1))

psd, freq,psdfilt,freqfilt=SVBfunc.FFRQ(Wdif[inds],Wfilt,TIME,dist) 


# In[9]:


psdeta, freqeta,psdfilteta,freqfilteta=SVBfunc.FFRQ(etadif[inds],etafilt,TIMEeta,disteta) 


# In[10]:


timepsd=np.arange(2880,TIME[-1],40)
timepsd=np.append(timepsd,14400)
periodsplot=(1/freq)/(60*60*24)
periodsplotfilt=(1/freqfilt)/(60*60*24)


# In[11]:
hejeta=[35,55,80,120,154,195,220] #30.706,31.1276,32.5,32.95
hej=np.zeros(len(hejeta))
for i in range(len(hejeta)):
	hej[i]=np.where(lat_ac==lat_aceta[hejeta[i]])[0][0]
	#hej[i]=np.where(lat_ac==indYlat[i][0]+1)[0][0]
	

hej=hej.astype(int)
#hej=[50,75,185,215] #30.706,31.1276,32.5,32.95
#hej=[58, 85, 205, 227]
markers=Line2D.filled_markers
markers=np.delete(markers,np.arange(2,5,1))
colors=colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999']

fig = plt.figure()
gs = GridSpec(nrows=1, ncols=2, width_ratios=[1.1, 0.9])


vmin=0
vmax=9
psdplot=np.transpose(psd[:,1:])
title='PSD [$10^{-4}$  ms$^{-1}$Hz$^{-1}$]'
xlab='Frequency [days$^{-1}$]'
ylab='Distance [km]'

#time= TIME.values.astype(int)
#time=dat.date2num(TIME)
markers=Line2D.filled_markers
markers=np.delete(markers,np.arange(2,5,1))
colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999']

ax = fig.add_subplot(gs[0, 0])
ax.set(xlabel=xlab, ylabel=ylab)
ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
ax.set_xticklabels(['0$^{-1}$', '1$^{-1}$', '2$^{-1}$', '3$^{-1}$','4$^{-1}$', '5$^{-1}$', '6$^{-1}$', '7$^{-1}$', '8$^{-1}$'])
ax.text(-0.08, 1.05, '(a)', fontweight='bold', color='k', 
        transform=ax.transAxes)
#val=VAL.values
cax = ax.pcolormesh(periodsplotfilt[1,:],dist,psdfilt*1e4,vmin=vmin,vmax=vmax, cmap=cmocean.cm.balance)

cb = plt.colorbar(cax)

cb.set_label(title)

ax1 = fig.add_subplot(gs[0, 1])
sort=np.argsort(periodsplotfilt[1,:])
p=0
for i in hej:
    p=p+1
    ax.axhline(y=dist[i],color=colors[p-1],linewidth=2)
    plt.plot(periodsplotfilt[1,sort],psdfilt[i,sort]*1e4, label=f'{dist[i]:.0f} km', color=colors[p-1])


ax1.legend()

ax1.set(xlabel=xlab, ylabel=title)
ax1.yaxis.tick_right()
ax1.yaxis.set_label_position("right")
#ax1.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
#ax1.set_xticklabels(['1/0', '1/1', '1/2', '1/3','1/4', '1/5', '1/6', '1/7', '1/8'])
ax1.text(-0.08, 1.05, '(b)', fontweight='bold', color='k', 
        transform=ax1.transAxes)
plt.show()


# In[15]:
'''

fig = plt.figure()
gs = GridSpec(nrows=1, ncols=2, width_ratios=[1, 1])


vmin=0
vmax=9
psdplot=np.transpose(psd[:,1:])
title='PSD [$10^{-4}$  ms$^{-1}$Hz$^{-1}$]'
titleeta='PSD [ mm$^{-1}$Hz$^{-1}$]'
xlab='Frequency [days$^{-1}$]'
ylab='Distance [km]'

#time= TIME.values.astype(int)
#time=dat.date2num(TIME)
markers=Line2D.filled_markers
markers=np.delete(markers,np.arange(2,5,1))
colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c']

ax = fig.add_subplot(gs[0, 0])
ax.set(xlabel=xlab, ylabel=ylab)
ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
ax.set_xticklabels(['0$^{-1}$', '1$^{-1}$', '2$^{-1}$', '3$^{-1}$','4$^{-1}$', '5$^{-1}$', '6$^{-1}$', '7$^{-1}$', '8$^{-1}$'])
ax.text(-0.08, 1.05, '(a)', fontweight='bold', color='k', 
        transform=ax.transAxes)
#val=VAL.values
cax = ax.pcolormesh(periodsplotfilt[1,:],dist,psdfilt[:,1:]*1e4,vmin=vmin,vmax=vmax, cmap=cmocean.cm.balance)

cb = plt.colorbar(cax)

cb.set_label(title)
vmineta=0
vmaxeta=20
ax1 = fig.add_subplot(gs[0, 1])
ax1.set(xlabel=xlab, ylabel=ylab)
ax1.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
ax1.set_xticklabels(['0$^{-1}$', '1$^{-1}$', '2$^{-1}$', '3$^{-1}$','4$^{-1}$', '5$^{-1}$', '6$^{-1}$', '7$^{-1}$', '8$^{-1}$'])
ax1.text(-0.08, 1.05, '(b)', fontweight='bold', color='k', 
        transform=ax1.transAxes)
#val=VAL.values
cax = ax1.pcolormesh(periodsplotfilt[1,:],disteta,psdfilteta[:,1:]*1e3,vmin=vmineta,vmax=vmaxeta, cmap=cmocean.cm.balance)

cb = plt.colorbar(cax)

cb.set_label(titleeta)

fig.tight_layout()

'''

times = (TIME)*60 #in s

t0 = 172800 # start is day 2 
dt1 = 600 # 20 min
dt2 = 1200 # 20 min 
freqlim = (1./dt)

nx = len(dist)
nt = int(len(times)/2)+2

psd = np.zeros((nx,nt))*np.nan
phase = np.zeros((nx,nt))*np.nan

psdfilt = np.zeros((nx,nt))*np.nan
phasefilt = np.zeros((nx,nt))*np.nan

freq =  np.zeros((nx,nt))*np.nan
freqfilt = np.zeros((nx,nt))*np.nan


ii=0
signalFFT1 = np.fft.rfft(Wdif[:432,ii])
signalFFTfilt1 = np.fft.rfft(Wfilt[:432,ii])

## Get Power Spectral Density
signalPSD1 = np.abs(signalFFT1)
signalPSDfilt1= np.abs(signalFFTfilt1)

## Get frequencies corresponding to signal 
fftFreq1 = np.fft.rfftfreq(len(Wdif[:432,ii]), dt1)
fftFreqfilt1 = np.fft.rfftfreq(len(Wfilt[:432,ii]), dt1)

signalFFT2 = np.fft.rfft(Wdif[432:,ii])
signalFFTfilt2 = np.fft.rfft(Wfilt[432:,ii])

## Get Power Spectral Density
signalPSD2 = np.abs(signalFFT2)
signalPSDfilt2 = np.abs(signalFFTfilt2)

## Get frequencies corresponding to signal 
fftFreq2 = np.fft.rfftfreq(len(Wdif[432:,ii]), dt2)
fftFreqfilt2 = np.fft.rfftfreq(len(Wfilt[432:,ii]), dt2)

psd[ii,:len(signalPSD1)] = signalPSD1[:]
psd[ii,len(signalPSD1):] = signalPSD2[:]
psdfilt[ii,:len(signalPSDfilt1)] = signalPSDfilt1[:]
psdfilt[ii,len(signalPSDfilt1):] = signalPSDfilt2[:]
freq[ii,:len(fftFreq1)] =  fftFreq1[:]
freq[ii,len(fftFreq1):] =  fftFreq2[:]
freqfilt[ii,:len(fftFreqfilt1)] = fftFreqfilt1[:]
freqfilt[ii,len(fftFreqfilt1):] = fftFreqfilt2[:]

