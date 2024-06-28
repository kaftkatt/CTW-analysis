import SVBfunc
dsmoov=SVBfunc.sgolay2d ( depth.values, window_size=29, order=4, derivative=None)
dsmoov2=SVBfunc.sgolay2d ( depth.values, window_size=29, order=8, derivative=None)
dsmoov3=SVBfunc.sgolay2d ( depth.values, window_size=199, order=4, derivative=None)
     
fig,ax=plt.subplots()
ax.pcolormesh(dsmoov)
fig1,ax1=plt.subplots()
ax1.pcolormesh(dsmoov2)
fig2,ax2=plt.subplots()
ax2.pcolormesh(dsmoov3)
fig3,ax3=plt.subplots()
ax3.pcolormesh(depth)



fig,ax=plt.subplots()
ax.plot(dsmoov[4,:])
fig1,ax1=plt.subplots()
ax1.plot(dsmoov2[4,:])
fig2,ax2=plt.subplots()
ax2.plot(dsmoov3[4,:])
fig3,ax3=plt.subplots()
ax3.plot(depth[4,:])
 
