from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as pl
import matplotlib.pyplot as plt
import cmocean
from scipy.io import loadmat
import numpy as np
from SVBfunc import loadNetCDFs 


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
