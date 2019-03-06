#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 15:37:42 2019
get the data from GoMOFs

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
import sys
import numpy as np
import math
import time


def get_gomofs_url(date):
    """
    the format of date is:datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    input date and return the url of data
    """
    
    print('start calculate the url!') 
    #below is a method to calculate this part of string in url('/201902/nos.gomofs.fields.n006.20190227.t12z.nc')
    date=date+datetime.timedelta(hours=4.5) #4.5 is a parameter use to get the correct time string("n003","t12") 
    date_ymdh=date.strftime('%Y%m%d%H%M%S')  #timestring,include year,month, day, hour, minetes, seconds
    ym=date_ymdh[:6]   #timestring,include year,month
    ymd=date_ymdh[:8] #timestring,include year,month,day,hour
    hours=int(date_ymdh[8:10])+int(date_ymdh[10:12])/60.+int(date_ymdh[12:14])/3600.  #hours,mintes,seconds change to hours  for examole: 1 h 30 min 0 sec=1.5h
    t=math.floor((hours)/6.0)*6
    if len(str(t))==1:
        tstr='t0'+str(t)+'z'
    else:
        tstr='t'+str(t)+'z'
    if round((hours)/3.0-1.5,0)==t/3:
        nstr='n006'
    else:
        nstr='n003'
    url='http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/GOMOFS/MODELS/'\
    +ym+'/nos.gomofs.fields.'+nstr+'.'+ymd+'.'+tstr+'.nc'
    return url

def get_gomofs(date_time,lat,lon,depth,mindistance=20):
    """
    the format time is: datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    lat and lon use decimal degrees
    if the depth is under the water, please must add the marker of '-'
    input time,lat,lon,depth return the temperature of specify location (or return temperature,nc,rho_index,ocean_time_index)
    
    return the temperature of specify location
    """
    #the data start time is '2018-07-01 00:00:00'
    if date_time<datetime.datetime.strptime('2018-07-01 00:00:00','%Y-%m-%d %H:%M:%S'):
        print('Time out of range')
        sys.exit()
        
        
    #start upload data
    check,count=0,1
    while(check==0):  #upload the data, if upload failed, re_upload several times
        try:
            url=get_gomofs_url(date_time)
            print('calculate the url finished!')
            nc=netCDF4.Dataset(str(url))
            gomofs_lons=nc.variables['lon_rho'][:]
            gomofs_lats=nc.variables['lat_rho'][:]
            gomofs_temp=nc.variables['temp'][:]
            gomofs_h=nc.variables['h'][:]
            gomofs_rho=nc.variables['s_rho'][:]
            check=1    #if upload succeed, end loop 
        except:
            time.sleep(30)
            count=count+1
            print('the '+str(int(count))+' times to upload data.')
        if count==5:   #loop five times will end the loop 
            print('end the loop to upload data, please check the website or file is exist!')
            sys.exit()
   

    """
    #upload the data
    try:
        url=get_gomofs_url(time)
        print('calculate the url finished!')
        nc=netCDF4.Dataset(str(url))
    except:
        print('please check the website or internet!')
    gomofs_lons=nc.variables['lon_rho'][:]
    gomofs_lats=nc.variables['lat_rho'][:]
    gomofs_temp=nc.variables['temp'][:]
    gomofs_h=nc.variables['h'][:]
    gomofs_rho=nc.variables['s_rho'][:]
    """
    
    #caculate the index of the nearest four points    
    print('start caculate the nearest four points!')
    eta_rho,xi_rho=0,0  #specify the initial index of location
    eta_rho1,xi_rho1=0,1
    eta_rho2,xi_rho2=1,0
    eta_rho3,xi_rho3=1,1
    distance0,distance1,distance2,distance3=5000,5000,5000,5000 #specify the initial value of distance
    for i in range(777):
        for j in range(1173):
            distance=zl.dist(lat1=lat,lon1=lon,lat2=gomofs_lats[i][j],lon2=gomofs_lons[i][j])
            if distance<=distance0:
                distance3,distance2,distance1,distance0=distance2,distance1,distance0,distance
                eta_rho3,xi_rho3,eta_rho2,xi_rho2,eta_rho1,xi_rho1=eta_rho2,xi_rho2,eta_rho1,xi_rho1,eta_rho,xi_rho
                eta_rho,xi_rho=i,j
            elif distance<=distance1:
                distance3,distance2,distance1=distance2,distance1,distance
                eta_rho3,xi_rho3,eta_rho2,xi_rho2=eta_rho2,xi_rho2,eta_rho1,xi_rho1
                eta_rho2,xi_rho2=i,j
            elif distance<=distance2:
                distance3,distance2=distance2,distance
                eta_rho3,xi_rho3=eta_rho2,xi_rho2
                eta_rho2,xi_rho2=i,j
            elif distance<=distance3:
                distance3=distance
                eta_rho3,xi_rho3=i,j
            else:
                continue
    if distance0>mindistance:
        print('THE location is out of range')
        sys.exit()
    # estimate the bottom depth of point location  
    print('start caculate the bottom depth of point location!') 
    points_h=[[gomofs_lats[eta_rho][xi_rho],gomofs_lons[eta_rho][xi_rho],gomofs_h[eta_rho][xi_rho]],
             [gomofs_lats[eta_rho1][xi_rho1],gomofs_lons[eta_rho1][xi_rho1],gomofs_h[eta_rho1][xi_rho1]],
             [gomofs_lats[eta_rho2][xi_rho2],gomofs_lons[eta_rho2][xi_rho2],gomofs_h[eta_rho2][xi_rho2]],
             [gomofs_lats[eta_rho3][xi_rho3],gomofs_lons[eta_rho3][xi_rho3],gomofs_h[eta_rho3][xi_rho3]]]
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
    print('start caculate the temperature of point location!')
    points_temp=[[gomofs_lats[eta_rho][xi_rho],gomofs_lons[eta_rho][xi_rho],gomofs_temp[0][rho_index][eta_rho][xi_rho]],
             [gomofs_lats[eta_rho1][xi_rho1],gomofs_lons[eta_rho1][xi_rho1],gomofs_temp[0][rho_index][eta_rho1][xi_rho1]],
             [gomofs_lats[eta_rho2][xi_rho2],gomofs_lons[eta_rho2][xi_rho2],gomofs_temp[0][rho_index][eta_rho2][xi_rho2]],
             [gomofs_lats[eta_rho3][xi_rho3],gomofs_lons[eta_rho3][xi_rho3],gomofs_temp[0][rho_index][eta_rho3][xi_rho3]]]
    temperature=zl.fitting(points_temp,lat,lon)
    # if input depth out of the bottom, print the prompt message
    if depth!='bottom':
        if abs(point_h)<abs(depth):
            print ("the depth is out of the bottom:"+str(point_h))
            sys.exit()
    return temperature,nc,rho_index,eta_rho,xi_rho
 #   return temperature

def countour_depth_temp_gomfs(output_path,date_time,lat=41.784712,lon=-69.231081,depth='bottom',addlon=.3,addlat=.3,mod_points='yes',depth_contours_interval=[20, 50,100,150,200,500]):
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
    '''
    #prepare the data
    temperature,nc,rho_index,eta_rho,xi_rho=get_gomofs(date_time,lat,lon,depth)
    gomofs_lons=nc.variables['lon_rho'][:]
    gomofs_lats=nc.variables['lat_rho'][:]
    gomofs_temp=nc.variables['temp'][:]
    gomofs_h=nc.variables['h'][:]
    '''
    
    #prepare the data
    temperature,nc,rho_index,eta_rho,xi_rho=get_gomofs(date_time,lat,lon,depth)
    check,count=0,0
    while(check==0):  #upload the data, if upload failed, re_upload several times
        try:
            gomofs_lons=nc.variables['lon_rho'][:]
            gomofs_lats=nc.variables['lat_rho'][:]
            gomofs_temp=nc.variables['temp'][:]
            gomofs_h=nc.variables['h'][:]
            check=1
        except:
            time.sleep(30)
            count=count+1 
            print('the '+str(int(count))+' times to upload data.')
    
    #creat map   
    print('start draw map!')   
    #Create a blank canvas
    fig=plt.figure(figsize = (20, 20))
    fig.suptitle('GoMOFs model bottom temp(deg C) and depth(meter)',fontsize=35, fontweight='bold')
    #Draw contour lines and temperature maps in detail
    ax1=fig.add_axes([0.07,0.03,0.85,0.95])
    ax1.set_title(date_time.strftime('%Y-%m-%d %H:%M:%S'), loc='center')
    ax1.axes.title.set_size(24)
    while(not zl.isConnected()):
        time.sleep(120)
    check,count=0,0
    service = 'Ocean_Basemap'
    while(check==0):
        try:
            map=Basemap(llcrnrlat=lat-addlat,urcrnrlat=lat+addlat,llcrnrlon=lon-addlon,urcrnrlon=lon+addlon,\
                        resolution='f',projection='tmerc',lat_0=lat,lon_0=lon,epsg = 4269)
            map.arcgisimage(service=service, xpixels = 5000, verbose= False)
            check=1
        except:
            check=0
            count=count+1
            print('start '+str(count)+' times add arcgisimage!')
    if mod_points=='yes':
        x1,y1=map(gomofs_lons,gomofs_lats)
        ax1.plot(x1,y1,'yo',markersize=0.1)
    #label the latitude and longitude
    parallels = np.arange(0.,90,1)
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=20,linewidth=0.0)
    # draw meridians
    meridians = np.arange(180.,360.,1)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=20,linewidth=0.0)
    m_lon,m_lat=map(gomofs_lons,gomofs_lats)
    print('start contour depth!')
    dept_clevs=depth_contours_interval
    dept_cs=map.contour(m_lon,m_lat,gomofs_h,dept_clevs,colors='black')
    plt.clabel(dept_cs, inline = True, fontsize =20,fmt="%1.0f")
    print('start contour temperature!')
    temp_cs=map.contourf(m_lon,m_lat,gomofs_temp[0][rho_index],7)
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
    ax1.plot(x,y,'ro')
    ax1.text(x[0]+0.02,y[0]-0.01,citys[0],bbox=dict(facecolor='yellow',alpha=0.5),fontsize =30)
    #indert a map that have mmore screen 
    ax2=fig.add_axes([0.09,0.68,0.2,0.2])
    #Build a map background
    print('start draw insert map!')
    while(not zl.isConnected()):
        time.sleep(120)
    check,count=0,0
    while(check==0):
        try:
            map2=Basemap(llcrnrlat=int(lat)-5,urcrnrlat=int(lat)+5,llcrnrlon=int(lon)-5,urcrnrlon=int(lon)+5,\
                         resolution='f',projection='tmerc',lat_0=int(lat),lon_0=int(lon),epsg = 4269)
            map2.arcgisimage(service=service, xpixels = 5000, verbose= False)
            check=1
        except:
            check=0
            count=count+1
            print('start '+str(count)+' times add arcgisimage!')

    parallels = np.arange(0.,90.,3)
    map2.drawparallels(parallels,labels=[0,1,0,0],fontsize=10,linewidth=0.0)
    # draw meridians
    meridians = np.arange(180.,360.,3)
    map2.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10,linewidth=0.0)
    x2,y2=map2(lon_point,lat_point)
    ax2.plot(x2,y2,'ro')
#    plt.savefig(output_path+'contour_depth_tem_GoMOFs.png',dpi=300)
    plt.show()
