from xmitgcm import open_mdsdataset
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from SVBfunc import loadNetCDFs

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import pylab as pl
import cartopy.mpl.geoaxes

coast='original'
linesperpcoast=1

levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
          11, 12, 13, 14, 15, 16, 17,
          18, 19, 20, 21, 22, 23, 24, 25,
          26, 27, 28, 29, 30, 31,
          32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
          45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
          58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
          74, 79, 84, 89, 94, 99, ]
          
file = sio.loadmat('/home/athelandersson/CTW-analysis/Files/Older/N2_lin.mat')
N2 = file['N2']
z = np.arange(0, -3500, -10)
if coast == 'smooth':
	pathw = '/data/SO2/sio-kramosmusalem/exp11_512x612x100_smooth_SVB/01_febTS_1000x'
	pathn = '/data/SO2/sio-kramosmusalem/exp11_512x612x100_smooth/01_febTS_1000x'
	data_dirWITH = '/data/SO2/sio-kramosmusalem/exp11_512x612x100_smooth_SVB/01_febTS_1000x'
elif coast == 'original':
	pathw = '/data/SO2/sio-kramosmusalem/exp06_512x612x100_ORL_SVB/01_SVB_febTS_output/'
	pathn = '/data/SO2/sio-kramosmusalem/exp06_512x612x100_ORL/01_noSVB_febTS/'
	data_dirWITH = '/data/SO2/sio-kramosmusalem/exp06_512x612x100_ORL_SVB/01b_SVB_febTS_output/'


dsw = open_mdsdataset(pathw, pathw, prefix=['eta'], default_dtype='>f4', levels=levels)
dsn = open_mdsdataset(pathn, pathn, prefix=['eta'], default_dtype='>f4', levels=levels)

LAT = dsw.YC
LON = dsw.XC - 360
Z = dsw.Z
hFacC = dsw.hFacC
hfa = np.ma.masked_values(hFacC[0, :, :], 0)
mask = np.ma.getmask(hfa)


dswr = open_mdsdataset(data_dirWITH, data_dirWITH, prefix=['rhoRef'], default_dtype='>f4', levels=levels)

depth = dsw.Depth
depthno = dsn.Depth
rho = dswr.rhoRef

ind_lon = [-115.11813068276555, -115.939167, -116.605833, -117.1625, -118.24368, -119.714167, -120.471439,
           -120.7586085906775]
ind_lat = [27.850440699318973, 30.556389, 31.857778, 32.715, 34.05223, 34.425833, 34.448113, 35.17364705813524]

matfile=sio.loadmat('/home/athelandersson/CTW-analysis/Files/' + str(coast) + '/BT_P_res30.mat')
lon,lat=matfile['lon'][0],matfile['lat'][0]


params = {'font.size': 16,
            'figure.figsize': (14, 8),
            'font.family': 'sans'}
pl.rcParams.update(params)
  
fig = plt.figure()
gs = GridSpec(nrows=2, ncols=2, width_ratios=[2.5, 0.5])


ax = fig.add_subplot(gs[0:, 0])
ax.set_facecolor('tan')

pc = ax.contourf(LON, LAT, np.ma.masked_array(depth, mask=mask), 50,
                 vmin=0, vmax=5000, cmap=cmo.cm.deep, zorder=1) 

cb = plt.colorbar(pc)  
cn = ax.contour(LON, LAT, depth, colors=['0.2', '0.4', '0.6', '0.8'],
                levels=[200, 500, 1000, 2000], zorder=2)
cb.set_label('Depth [m]')
ax.contour(LON, LAT, depthno, levels=[0], colors='brown', linestyles=':', linewidths=2.5, zorder=3)

axins = inset_axes(ax, width="28%", height="28%", loc='upper right',
                   axes_class=cartopy.mpl.geoaxes.GeoAxes,
                   axes_kwargs=dict(map_projection=cartopy.crs.Orthographic(central_latitude=32,
                                                                            central_longitude=-118)))
axins.add_feature(cartopy.feature.OCEAN, zorder=0)
axins.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
axins.gridlines()
axins.stock_img()

axins.plot(LON, np.ones_like(LON) * 35.3, 's', color='r', markersize=14, fillstyle='none',
           transform=cartopy.crs.Orthographic(-118, 30.7))

if linesperpcoast==1:
	colors=[ '#4daf4a', '#a65628', '#984ea3',
                   '#e41a1c', '#dede00','#377eb8'
       ,'#ff7f00','#f781bf','#999999','tab:blue']

	for i in range(len(lon)):
		ax.scatter(lon[i][0],lat[i][0],color=colors[i],linewidth=1, zorder=4)


for kk, ll, lab in zip(ind_lon, ind_lat,
                       ['Punta \n Eugenia', 'San Quintín', 'Ensenada', 'San Diego', 'Los Angeles', 'Santa Barbara',
                        'Point Conception', 'Port San Luis']):
	ax.plot(kk, ll, 'o',markersize=14, color='r', markeredgecolor='k')
	if lab == 'Point Conception':
		ax.text(kk - 0.06, ll + 0.25, lab, fontsize=18)
	elif lab == 'Santa Barbara':
		ax.text(kk + 0.2, ll - 0.05, lab, fontsize=18)
	elif lab == 'Port San Luis':
		ax.text(kk + 0.2, ll - 0.1, lab, fontsize=18)
	elif lab == 'Punta \n Eugenia':
		ax.text(kk + 0.6, ll - 0.1, lab, fontsize=18, horizontalalignment='center')
	else:
		ax.text(kk + 0.16, ll - 0.05, lab, fontsize=18)


ax.text(0.93, 0.17, 'SVB', fontsize=24, horizontalalignment='center', fontweight='bold',
        transform=ax.transAxes)
        
ax.plot(LON[45], LAT[40], '*', color='gold', markersize=20, markeredgecolor='w')

ax.text(0.08, 0.75, 'Santa \n Barbara \n Channel', color='w', fontsize=15, transform=ax.transAxes,
        fontweight='demibold', horizontalalignment='center')
ax.text(0.18, 0.62, 'Santa \n Rosa \n Ridge', color='w', fontsize=15, transform=ax.transAxes,
        fontweight='bold', horizontalalignment='center')
hej = [58, 85, 205, 227]

markers = Line2D.filled_markers
markers = np.delete(markers, np.arange(2, 5, 1))
colors = ['b', 'g', 'r', 'k', 'm', 'c', 'y', 'brown', 'lime', 'fuchsia', 'beige']

ax.set_xlabel('Lon [°]')
ax.set_ylabel('Lat [°]')
ax.set_xlim(238 - 360, 246 - 360)
ax.set_ylim(27, 35.3)
ax.set_aspect(1)
ax.text(-0.08, 1.05, '(a)', fontweight='bold', color='k',
        transform=ax.transAxes)

ax1 = fig.add_subplot(gs[0:, 1])
ax1.plot(rho[Z >= -1800], Z[Z >= -1800], 'tab:blue', linewidth='5')
ax1.set_xlabel('Density [kg/m³]', color='tab:blue', labelpad=10)
ax1.tick_params(axis='x', labelcolor='tab:blue', pad=0)
ax1.set_ylabel('Depth [m]')
ax1.set_ylim(-1500, 0)
ax1.yaxis.tick_right()
ax1.yaxis.set_label_position("right")

ax2 = ax1.twiny()
ax2.plot(N2[z >= -1500] * 10 ** 4, z[z >= -1500], ':', color='tab:red', linewidth='5')
ax2.tick_params(axis='x', labelcolor='tab:red', pad=0)
ax2.set_xlabel('N2 [$10^{-4} s^{-1}$]', color='tab:red', labelpad=10)
ax2.text(-0.3, 1.05, '(b)', fontweight='bold',
         color='k', transform=ax2.transAxes)
fig.tight_layout()

plt.savefig('/home/athelandersson/CTW-analysis/Figures/' + str(coast) + 'map' + str(linesperpcoast) + '.png')


