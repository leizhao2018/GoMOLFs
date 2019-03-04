#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 15:37:42 2019
get the data from gomofs

@author: leizhao
"""
import netCDF4
import datetime
import zlconversions as zl

import os
import conda
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import numpy as np


def get_gomofs_url(date):
    """
    the format of date is:datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    input date and return the url of data
    """
    date_ymdh=date.strftime('%Y%m%d%H')
    ym=date_ymdh[:6]
    ymd=date_ymdh[:8]
    hour=date_ymdh[8:10]
    url='http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/GOMOFS/MODELS/'\
    +ym+'/nos.gomofs.stations.forecast.'+ymd+'.t'+hour+'z.nc'
    return url

def get_gomofs(time,lat,lon,depth):
    """
    the format time is: datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    lat and lon use decimal degrees,minutes
    if the depth is under the water, please must add the marker of '-'
    input time,lat,lon,depth return the temperature of specify location (or return temperature,nc,rho_index,ocean_time_index)
    
    return the temperature of specify location
    """
    url=get_gomofs_url(time)

    try:
        nc=netCDF4.Dataset(str(url))
    except:
        print('please check the website or internet!')
    gomofs_lons=nc.variables['lon_rho'][:]
    gomofs_lats=nc.variables['lat_rho'][:]
    gomofs_temp=nc.variables['temp'][:]
    gomofs_time=nc.variables['ocean_time'][:]
    gomofs_h=nc.variables['h'][:]
    gomofs_rho=nc.variables['s_rho'][:]
#    time=datetime.datetime.strptime(time,'%Y-%m-%d %H:%M:%S')
    diff_time=datetime.timedelta(days=10) #specify the initial value of diff_time
    # find the nearest time of input time
    for i in range(len(gomofs_time)):
        ocean_time=datetime.datetime.strptime('2016-01-01 00:00:00','%Y-%m-%d %H:%M:%S')+datetime.timedelta(seconds=gomofs_time[i])
        if abs(diff_time)>abs(ocean_time-time):
            ocean_time_index=i
        else:
            ocean_time_index=i-1
        if ocean_time-time>datetime.timedelta(seconds=0):
            break
        else:
            diff_time=ocean_time-time
            
    #caculate the index of the nearest four points        
    location_index=0  #specify the initial index of location
    location_index1=1
    location_index2=2
    location_index3=3
    distance0,distance1,distance2,distance3=5000,5000,5000,5000 #specify the initial value of distance
    for j in range(len(gomofs_lats)):
        distance=zl.dist(lat1=lat,lon1=lon,lat2=gomofs_lats[j],lon2=gomofs_lons[j])
        if distance<=distance0:
            distance3,distance2,distance1,distance0=distance2,distance1,distance0,distance
            location_index3,location_index2,location_index1=location_index2,location_index1,location_index
            location_index=j
        elif distance<=distance1:
            distance3,distance2,distance1=distance2,distance1,distance
            location_index3,location_index2=location_index2,location_index1
            location_index1=j
        elif distance<=distance2:
            distance3,distance2=distance2,distance
            location_index3=location_index2
            location_index2=j
        elif distance<=distance3:
            distance3=distance
            location_index3=j
        else:
            continue
    # estimate the bottom depth of point location       
    points_h=[[gomofs_lats[location_index],gomofs_lons[location_index],gomofs_h[location_index]],
             [gomofs_lats[location_index1],gomofs_lons[location_index1],gomofs_h[location_index1]],
             [gomofs_lats[location_index2],gomofs_lons[location_index2],gomofs_h[location_index2]],
             [gomofs_lats[location_index3],gomofs_lons[location_index3],gomofs_h[location_index3]]]
    point_h=zl.fitting(points_h,lat,lon) 
    # caculate the rho index
    if depth=='bottom':
        rho_index=0
    else:
        distance_h=gomofs_rho[0]*point_h-depth      
        for k in range(len(gomofs_rho)):
            if abs(distance_h)>=abs(gomofs_rho[k]*point_h-depth):
                distance_h=gomofs_rho[k]*point_h-depth
                rho_index=k
    #estimate the temperature of point location
    points_temp=[[gomofs_lats[location_index],gomofs_lons[location_index],gomofs_temp[ocean_time_index][location_index][rho_index]],
             [gomofs_lats[location_index1],gomofs_lons[location_index1],gomofs_temp[ocean_time_index][location_index1][rho_index]],
             [gomofs_lats[location_index2],gomofs_lons[location_index2],gomofs_temp[ocean_time_index][location_index2][rho_index]],
             [gomofs_lats[location_index3],gomofs_lons[location_index3],gomofs_temp[ocean_time_index][location_index3][rho_index]]]
    temperature=zl.fitting(points_temp,lat,lon)
    # if input depth os out of the bottom, print the prompt message
    if depth!='bottom':
        if abs(point_h)<abs(depth):
            print ("the depth is out of the bottom:"+str(point_h))
            temperature=9999
    return temperature
#    return temperature,nc,rho_index,ocean_time_index


#
#input_date_time='2018-11-11 12:00:00'
#input_lat=41.784712
#input_lon=-69.231081
def countour_depth_temp_gomfs(output_path,date_time,lat=41.784712,lon=-69.231081,depth='bottom',addlon=3,addlat=3,mod_points='yes',depth_contours_interval=[20, 50,100,150,200,500]):
    """Draw contours and isothermal layers on the map
    feb 11:Geographic box with full map in “insert” use “addlon=0.3”
    Allow users hardcode at top of code to defive depth_contours_to_plot=[20, 50,100]
    Allow users option of posting model node points
    notice:
    addlon,addlat: edges around point to include in the zoomed in plot
    depth: if depth=='bottom' print the bottom data.
    the format time is: datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    mod_points:do you want to post model grid nodes,if mod_points='yes', print model grid nodes;if other string, skip 
    
    
    """
    temperature,nc,rho_index,ocean_time_index=get_gomofs(date_time,lat,lon,depth)
    gomofs_lons=nc.variables['lon_rho'][:]
    gomofs_lats=nc.variables['lat_rho'][:]
    gomofs_temp=nc.variables['temp'][:]
#    gomofs_time=nc.variables['ocean_time'][:]
    gomofs_h=nc.variables['h'][:]
#    gomofs_rho=nc.variables['s_rho'][:]
    
    
#    xData = list(gomofs_lats)
#    x = list(set(xData))
#    x.sort(key=xData.index)
# 
#    # y
#    yData = list(gomofs_lons)
#    y = list(set(yData))
#    y.sort(key=yData.index)
# 
##    # to X and Y
##    X, Y = np.meshgrid(y, x)
# 
    
    lats,lons,Z,TEMP=[],[],[],[]
    for i in range(len(gomofs_lons)):
        if gomofs_lons[i]>100000 or gomofs_lats[i]>100000 or gomofs_temp[ocean_time_index][i][rho_index]=='--':
            continue
        else:
            lats.append(gomofs_lats[i])
            lons.append(gomofs_lons[i])
            Z.append(gomofs_h[i])
            TEMP.append(gomofs_temp[ocean_time_index][i][rho_index])
    lats=np.array(lats)
    lons=np.array(lons)
    Z=np.array(list(Z))
    TEMP=np.array(list(TEMP))

    
    #creat map
    #Create a blank canvas       
    fig=plt.figure(figsize = (20, 20))
    fig.suptitle('GoMOFs model bottom temp(deg C) and depth(meter)',fontsize=35, fontweight='bold')

    #Draw contour lines and temperature maps in detail
    ax1=fig.add_axes([0.07,0.03,0.85,0.95])
    ax1.set_title(date_time.strftime('%Y-%m-%d %H:%M:%S'), loc='center')
    ax1.axes.title.set_size(24)
    
    service = 'Ocean_Basemap'
    map=Basemap(llcrnrlat=lat-addlat,urcrnrlat=lat+addlat,llcrnrlon=lon-addlon,urcrnrlon=lon+addlon,\
            resolution='f',projection='tmerc',lat_0=lat,lon_0=lon,epsg = 4269)
    map.arcgisimage(service=service, xpixels = 5000, verbose= False)

    #label the latitude and longitude
    parallels = np.arange(0.,90,1)
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=20,linewidth=0.0)
    # draw meridians
    meridians = np.arange(180.,360.,1)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=20,linewidth=0.0)
    m_lon,m_lat=map(lons,lats)

    dept_clevs=depth_contours_interval
    dept_cs=map.contour(m_lon,m_lat,[Z],dept_clevs,colors='black')
    plt.clabel(dept_cs, inline = True, fontsize =20,fmt="%1.0f")
    
    temp_cs=map.contourf(m_lon,m_lat,TEMP,7)
    temp_cbar=map.colorbar(temp_cs,location='right',size="5%",pad="1%")
    temp_cbar.set_label('deg C',size=25)
    temp_cbar.ax.set_yticklabels(temp_cbar.ax.get_yticklabels(), fontsize=20)
    if temperature==9999:
        citys=['the depth is out of the bottom depth']
    else:
        citys=[str(round(temperature,2))]
    lat_point=[lat]
    lon_point=[lon]
    x,y=map(lon_point,lat_point)
    if mod_points=='yes':
        x1,y1=map(gomofs_lons,gomofs_lats)
        ax1.plot(x1,y1,'yo')
    ax1.plot(x,y,'ro')
    ax1.text(x[0]+0.02,y[0]-0.01,citys[0],bbox=dict(facecolor='yellow',alpha=0.5),fontsize =30)
    #indert a map that have mmore screen 
    ax2=fig.add_axes([0.09,0.68,0.2,0.2])
    #Build a map background
    map1=Basemap(llcrnrlat=int(lat)-5,urcrnrlat=int(lat)+5,llcrnrlon=int(lon)-5,urcrnrlon=int(lon)+5,\
            resolution='f',projection='tmerc',lat_0=int(lat),lon_0=int(lon),epsg = 4269)
    map1.arcgisimage(service=service, xpixels = 5000, verbose= False)


    parallels = np.arange(0.,90.,3)
    map1.drawparallels(parallels,labels=[0,1,0,0],fontsize=10,linewidth=0.0)
    # draw meridians
    meridians = np.arange(180.,360.,3)
    map1.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10,linewidth=0.0)
    x,y=map1(lon_point,lat_point)
    #Draw contour lines and temperature maps
    ax2.plot(x,y,'ro')
#    plt.savefig(output_path+'contour_depth_tem_doppio2.png',dpi=300)
    plt.show()
#    return 1


def transform(x,y,z):
    xData = list(x)
    x = list(set(xData))
    x.sort(key=xData.index)
 
    # y
    yData = list(y)
    y = list(set(yData))
    y.sort(key=yData.index)
 
    # to X and Y
    X, Y = np.meshgrid(y, x)
 
    # z
    z = z
    Z = np.array(z)
    Z = Z.reshape((len(x), len(y)))
    




#output_path='/home/jmanning/Desktop/testout/doppio/'
#date_time=datetime.datetime.strptime('20180331 060000','%Y%m%d %H%M%S')
interval=[20,50,100,150,200,500]
#a=countour_depth_temp_gomfs(output_path,date_time,lat=41.784712,lon=-69.231081,depth='bottom',addlon=3,addlat=3,mod_points='yes',depth_contours_interval=interval)



#hardcodes
time=datetime.datetime.strptime('20180331 060000','%Y%m%d %H%M%S')
lat=41
lon=-71
depth=-900
###############
##
##
#
temp,nc,rho_index,ocean_time_index=get_gomofs(time,lat,lon,depth)
print(temp)
gomofs_lons=nc.variables['lon_rho'][:]
gomofs_lats=nc.variables['lat_rho'][:]
gomofs_temp=nc.variables['temp'][:]
gomofs_time=nc.variables['ocean_time'][:]
gomofs_h=nc.variables['h']
gomofs_rho=nc.variables['s_rho'][:]
#
service = 'Ocean_Basemap'
xpixels = 5000 
fig=plt.figure(figsize=(8,8.5))
map=Basemap(projection='mill',llcrnrlat=38,urcrnrlat=43,llcrnrlon=-73,urcrnrlon=-68,\
                resolution='f',lat_0=42,lon_0=-67,epsg = 4269)
map.arcgisimage(service=service, xpixels = xpixels, verbose= False)
step=1.0
parallels = np.arange(0.,90.0,step)
map.drawparallels(parallels,labels=[0,1,0,0],fontsize=10,linewidth=0.2)
# draw meridians
meridians = np.arange(180.,360.,step)
map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10,linewidth=0.2)
tele_x,tele_y=map(gomofs_lons,gomofs_lats)
plt.plot(tele_x,tele_y,'b*',markersize=8,alpha=0.5,label='telemetry')
dept_clevs=interval
dept_cs=map.contourf(tele_x,tele_y,gomofs_h,colors='black')
plt.show()
plt.savefig('/home/jmanning/Desktop/gomofdataplot.png',dpi=500)

#http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/GOMOFS/MODELS/201803/nos.gomofs.stations.forecast.20180331.t06z.nc
#http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/GOMOFS/MODELS/201803/nos.gomofs.stations.nowcast.'+ymd+'.t'+hour+'z.nc



