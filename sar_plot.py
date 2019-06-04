import numpy as np
import netCDF4
import datetime as dt
import matplotlib.pyplot as plt
from cartopy import feature
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from cartopy.io.shapereader import Reader
from my_color import my_colors
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.image as image
import scipy.io as sio

#import cartopy.io.shapereader as shpreader
#f=shpreader.natural_earth(resolution='10m',
#        category='cultural',
#        name='populated_places',
#        )
#print(f)
#load logo
def draw_info(ax,title, txt_box,\
    left=0.70,bottom=0.00,w=0.30,h=0.13,\
    h_title=0.04, xpad=0.01, ypad=0.05, tpad=0.03):
    ''' Draw info box
    ax - axis
    title - title text
    txt_box - body text
    '''
    top=bottom+h
    right=left+w

    l_rec = mpatches.Rectangle(
            (left,bottom), w, h,
            fill=False, 
            transform=ax.transAxes,
            zorder=17,
                        )
    ax.add_patch(l_rec)
    l_line= mlines.Line2D([left,right],[top-h_title],
            color='black',
            transform=ax.transAxes,
            zorder=17,
            )
    ax.add_line(l_line)
    #props = dict(
    #        facecolor='white',
    #        alpha=1.0,
    #        )

    ax.text(left+xpad,top-tpad,
            title,
            transform=ax.transAxes,
            fontsize=14,
            fontfamily='serif',#'monospace',#'sans-serif',
            #bbox=props,
            )
    ax.text(left+xpad,top-ypad,
            txt_box,
            transform=ax.transAxes,
            verticalalignment='top',
            ha='left',
            #bbox=props,
            #color='white',
            ) 

def getdata(url):
    '''get common data from netCDF file
    return data dict
    '''
    data={}
    with netCDF4.Dataset(url) as nc_data:
        data['lon']=nc_data['lon'][:]-360
        data['lat']=nc_data['lat'][:]
        #triangle structure
        data['nv']=nc_data['nv'][:].T-1
        #triangle center
        data['lonc']=nc_data['lonc'][:]-360
        data['latc']=nc_data['latc'][:]
        #data
        #data=nc_data[var][:,0,:]
        #vector data
        data['u']=nc_data['u'][:,0,:]
        data['v']=nc_data['v'][:,0,:]
        data['magnitude'] = (data['u'] ** 2 + data['v'] ** 2) ** 0.5
        #time
        times=[]
        for i in nc_data['Times'][:]:
            times.append(i.tobytes().decode('utf-8'))
        data['time']= times
        #print(data['time'], )
        return data

def plot_base(data, itime, var_levels, oname):
    #land,coastline,states,reader=features
#Size
#WXGA size
#plt.figure(figsize=(12.80,7.20),dpi=100)
    t1= dt.datetime.now()
#FHD size
    fig=plt.figure(figsize=(19.20,10.80),dpi=100)

#Adjusting margins figure
    fig.subplots_adjust(left=0.05,right=0.95,top=0.97,bottom=0.03)
    #ax=fig.subplot(1,1,1,
    ax=fig.add_subplot(1,1,1,
            projection=proj,
            )
#load features
#ax.add_feature(ocean,
#        facecolor=feature.COLORS['water'],
#        )
    ax.add_feature(land,
            zorder=10,
            )
    ax.add_feature(feature.BORDERS,
            zorder=11,
            )
    ax.add_feature(coastline, 
        edgecolor='black',
        zorder=12,
        )
    ax.add_feature(states, 
        edgecolor='gray',
        zorder=12,
        )
    #print('features:',dt.datetime.now()-t1, end=' ')
#url = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
#layer = 'VIIRS_CityLights_2012'
#ax.add_wmts(url, layer)

    #city names
    #for i,record in enumerate(reader.records()):
    #    name = record.attributes['name_es']
    #    geometry = record.geometry
    #    #if i==0:
    #        #print(record.attributes.keys())
    #    x,y = geometry.centroid.x, geometry.centroid.y
    #    ax.text(x, y, name,
    #            #fontsize=2,
    #            zorder=17,
    #            clip_on=True,
    #            ha='center',
    #            va='top',
    #            transform=ccrs.PlateCarree(),
    #            color='gray',
    #            )
    #print('names:',dt.datetime.now()-t1, end=' ')

    #logo plot
    im = image.imread('logo_atmosfera.png')
    ax.imshow(im,
        aspect='equal',
        extent=(lon_min+0.1,lon_min+1.1, lat_max,lat_max-0.27),
        zorder=15,
        )
    #print('logo:',dt.datetime.now()-t1, end=' ')
    #info text
    txt_box='\n'.join((
            '>Distribución del sargazo + Temperatura ',
            #' Corriente en superficie',
            '>Modelo: FVCOM',
            '>Grupo Interacción Océano-Atmósfera',
            )
            )
#
    #bbox_props = dict(boxstyle="pad=0.8", fc="white", ec="b", lw=2)
    #info box plot
    draw_info(ax,title,txt_box)
    #print('box::',dt.datetime.now()-t1, end=' ')
#wind plot
#ind=np.arange(0,len(u),20)
#ax.quiver(
        #lonc,
        #latc,
        #u,
        #v,
        #scale=30,
        ##transform=ccrs.PlateCarree(),
        #regrid_shape=20,
        #)
#wplot=ax.streamplot(
#        lonc,
#        latc,
#        u,
#        v,
#        #linewidth=1,
#        density=2,
#        color=magnitude,
#        cmap=plt.cm.tab10,
#        #levels=np.arange(0,2,0.2),
#        zorder=5,
#        )
    #city points plot
    ax.scatter(xp,
               yp,
               100,
               marker='*',
               color='white',
               transform=ccrs.PlateCarree(),
               zorder=15,
               )
    #print('points:',dt.datetime.now()-t1, end=' ')
    #apply limits
    plt.axis(ax_lim)
    ax_leg=ax.gridlines(
            draw_labels=True,
            linestyle='--',
            )
    ax_leg.xlabels_top=False
    ax_leg.ylabels_right=False
    ax_leg.xlocator=mticker.FixedLocator(range(int(lon_min),int(lon_max)+1,1))
    ax_leg.xformatter=LONGITUDE_FORMATTER
    ax_leg.ylocator=mticker.FixedLocator(range(int(lat_min),int(lat_max)+2,1))
    ax_leg.yformatter=LATITUDE_FORMATTER

    #plot
    l_map=ax.tricontour(
            data['lon'],
            data['lat'],
            data['nv'],
            data['var'][itime],
            colors='gray',
            levels=var_levels,
        )
        #ax.clabel(z_map,
        #        )
    #print('lines:',dt.datetime.now()-t1, end=' ')
    var_map=plt.tricontourf(
        data['lon'],
        data['lat'],
        data['nv'],
        data['var'][itime],
        levels=var_levels,
        #transform=ccrs.Mercator(),
        #cmap=plt.cm.coolwarm,
        cmap=my_cbar,
        )
    #print('fill:',dt.datetime.now()-t1, end=' ')
    #print(data['time'][itime])
    plt.title(data['time'][itime])

#color bar
    cbar=plt.colorbar(
            var_map,
            #ax=ax,
            #ticks=t_ticks,
            ticks=var_levels,
            #shrink=0.75,
            #orientation='horizontal',
            aspect=50,
            pad=0.02,
            fraction=0.03,
            )
    cbar.ax.set_ylabel('Elevación de la superficie [m]')
    #mat plot
    ax.plot(
            data['par_lon'],
            data['par_lat'],
            'r,',
            #color='red',
            #marker='*',
            #markersize=12,
            transform=ccrs.PlateCarree(),
            zorder=1,
            )
    #print('mat:',dt.datetime.now()-t1, )
    plt.savefig('{}_{:03}.png'.format(oname,itime))
    return 0

import sys
if __name__=='__main__':
    #plot title (info box)
    var=sys.argv[1]
    #var='temp'
    title=sys.argv[2]
    #title='ELEVACIÓN DEL MAR'
    url=sys.argv[3]
    #url='lustre_test/Sargazo01_0001.nc'
    ofile=sys.argv[4]
    #ofilename='testfig'
    n_cpu=int(sys.argv[5])
#load color bar
    my_cbar=ListedColormap(my_colors)
#config vars
#projection
    proj=ccrs.PlateCarree()
    lat_min=17.0
    lat_max=22.8
    lon_min=-89
    lon_max=-80.5
#map limits
    ax_lim= [lon_min, lon_max, lat_min, lat_max]
#var levels
    var_levels={}

#zeta
    lev_a=np.arange(-1,0,0.1)
    lev_b=np.arange(-0.05,0.1,0.05)
    lev_c=np.arange(0.2,1,0.1)
#zlevels=np.concatenate((lev_a,lev_b,lev_c))
    var_levels['zeta']=np.concatenate((lev_a,lev_b,lev_c))

#temperature
#tlevels=np.arange(23,30,0.5)
    var_levels['temp']=np.arange(23,30,0.5)
    t_ticks=np.arange(20,30,0.5)

#load shapes
    states = NaturalEarthFeature(category="cultural", scale="10m",
            facecolor="none",
            name='admin_1_states_provinces_lines',
            )
    coastline = NaturalEarthFeature(category="physical", scale="10m",
            facecolor="none",
            name='coastline',
            )
    land= NaturalEarthFeature(category="physical", scale="10m",
            #facecolor=feature.COLORS['land'],
            facecolor='lightgray',
            name='land',
            )
    fname='files/ne_10m_populated_places.shp'
    reader=Reader(fname)
#city points
    points=[list(point.coords) for point in reader.geometries()]
    xp=[point[0][0] for point in points]
    yp=[point[0][1] for point in points]
    #features=[land,coastline,states,reader]
    
#ocean= feature.ShapelyFeature(
#        Reader('ne_10m_ocean.shp').geometries(),
#        ccrs.PlateCarree(),
#        )

    stime=dt.datetime.now()
    #load netcdf data
    data=getdata(url)
    with netCDF4.Dataset(url) as nc_data:
        if var == 'temp':
            #only 0 layer
            data['var']=nc_data[var][:,0,:]
        if var == 'zeta':
            #zeta:elevación de la superficie del mar
            data['var']=nc_data[var][:,:]

    #load mat file
    par_filename='1992-01-01_h_23.mat'
    data_stru=sio.loadmat(par_filename)['Particles_S']
    data['par_lat']=data_stru['Lat'][0,0][0]
    data['par_lon']=data_stru['Lon'][0,0][0]+10
    #print(par_lon.shape)
    #print(par_lon[0])
    #print(par_lon)
    t1=dt.datetime.now()-stime
    stime=dt.datetime.now()
    itime=0
    import multiprocessing as mp
    pars=[(data,itime,var_levels[var],ofile) for itime in range(1,121) ]
    #print(type(pars),len(pars),type(pars[3][4]))
    
    with mp.Pool(n_cpu) as pool:
        pool.starmap(plot_base, pars)
    t2=dt.datetime.now()-stime
    print('time:', t1, t2)

