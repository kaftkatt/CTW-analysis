sp=[]
shelf=np.zeros(len(dep))
shelfind=np.zeros(len(dep))	
for i in range(len(dep)):
	ind=np.where(dist[i]>=35)[0][0]	
	sp.append((dep[i][1:ind]-dep[i][:ind-1])/(dist[i][1:ind]-dist[i][:ind-1]))
	shelfind[i]=np.where(dep[i]>=300)[0][0]
	shelf[i]=dep[i][shelfind[i].astype(int)]
	
s=np.zeros(len(dep))

for i in range(len(dep)):
	depin=dep[i][:np.where(dist[i]>=40)[0][0]]
	distin=dist[i][:np.where(dist[i]>=40)[0][0]]	
	ho=np.mean(depin)
	lam=1/(distin[np.where(depin>=ho)[0][0]-1])
	hx=(depin[1:]-depin[:-1])/(distin[1:]-distin[:-1])
	
	s[i]=np.mean((hx/depin[:-1])*distin[:-1])
			
	fig=plt.figure()
	plt.plot(distin,-depin,color=colors[i],linewidth=2)
	plt.plot(distin[:-1],-((ho)*(lam*distin[:-1])**s[i]),'--')
	#plt.ylim([-3000,0])
	
plt.scatter(corrinds,s)
plt.show()

pathVEL='/home/athelandersson/NETCDFs/WVELAC.nc'
dsVEL= xr.open_dataset(pathVEL)

WVEL=dsVEL.ValfiltAll.values
distVEL=dsVEL.dist.values
TIMEVEL=dsVEL.time2.values
lat_acVEL=dsVEL.latAC.values
lon_acVEL=dsVEL.lonAC.values

hfa = np.ma.masked_values(hFacC[55, :, :], 0)
mask2 = np.ma.getmask(hfa)

fig = plt.figure()
gs = GridSpec(nrows=2, ncols=2)

ax = fig.add_subplot(gs[0, 0])
ax.set_facecolor('tan')
cax = ax.contourf(LON, LAT[lat_ac[:250]], np.ma.masked_array(depth[lat_ac[:250],:], mask=mask2[lat_ac[:250],:]), 50,vmin=0, vmax=5000, cmap=cmocean.cm.deep) 

 
cn = ax.contour(LON, LAT[lat_ac[:250]], depth[lat_ac[:250],:], colors=['0.2','0.4', '0.6', '0.8'],levels=[0,500, 1000, 2000])

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Depth [m]')
ax.text(-0.15, 1.2, (a), fontweight='bold', color='k', 
        transform=ax.transAxes)

ax.scatter(LON[lon_acVEL[hej2]],LAT[lat_acVEL[hej2]],s=0.5)
for i in range(len(hej2)):
	ax.axhline(y=LAT[lat_ac[hej2[i]]],color=colors[i],linewidth=2,alpha=0.4)

ax.set_xlabel('Lon [°]')
ax.set_ylabel('Lat [°]')

ax = fig.add_subplot(gs[0, 1])
vmin=-5
vmax=5
cbarall=1
SVBfunc.plot_HOVMOLLER(ax,distVEL,TIMEVEL,WVEL*1e6,'','Vertical velocity  [$10^{-6}$ ms$^{-1}$]',vmin,vmax,fig,lat_acVEL,lon_acVEL,1,cbarall,'(a)')
ax.set_ylim([0,distAC[250]])

for i in range(len(hej2)):#range(len(colors)):
	ax.axhline(y=distAC[hej2[i]],color=colors[i],linewidth=2,alpha=0.7)
	#ax.axhline(y=distAC[27*i],color=colors[i],linewidth=2,alpha=0.3)
'''
ax = fig.add_subplot(gs[0, 2])
ax.plot(s,distAC[:250])	
ax.set_ylim([0,distAC[250]])
ax.set(xlabel='s', ylabel='Distance from bay [km]')
'''

xlab='Distance from coast [km]'
ylab='Distance from bay [km]'

vmin=-200
vmax=200
ax = fig.add_subplot(gs[1, 0])
cax = ax.pcolormesh(dist[0][:ind-1],distAC[hej2],s,cmap=cmocean.cm.diff,vmin=vmin,vmax=vmax)
ax.scatter(dist[0][shelfind.astype(int)],distAC[hej2]) 


ax.set(xlabel=xlab, ylabel=ylab)

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Inclination of coast')
ax.text(-0.15, 1.2, '(c)', fontweight='bold', color='k', 
        transform=ax.transAxes)

for i in range(len(hej2)):
	ax.axhline(y=distAC[hej2[i]],color=colors[i],linewidth=2,alpha=0.7)

ax = fig.add_subplot(gs[1, 1])
cax = ax.pcolormesh(dist[0],distAC[:250],dep,cmap=cmocean.cm.diff,vmin=0,vmax=2000)
ax.scatter(dist[0][shelfind.astype(int)],distAC[:250]) 

divider = make_axes_locatable(ax)
axdiv = divider.new_vertical(size = '5%', pad = 0.5)
fig.add_axes(axdiv)
cbar_ax = plt.colorbar(cax, cax=axdiv,orientation='horizontal')
cbar_ax.ax.xaxis.set_label_position("top")
cbar_ax.set_label('Depth [m]')
ax.text(-0.15, 1.2, '(d)', fontweight='bold', color='k', 
        transform=ax.transAxes)


ax.set(xlabel=xlab, ylabel=ylab)
ax.set_xlim([0,dist[0][ind-1]])


colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999']

for i in range(250):
	ax.axhline(y=distAC[:250,color=colors[i],linewidth=2,alpha=0.7)

plt.savefig('/home/athelandersson/CTW-analysis/Figures/HovPerpInds.png')	



