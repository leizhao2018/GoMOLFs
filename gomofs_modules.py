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
import sys
import numpy as np


def get_gomofs_url(date,data_str):
    """
    the format of date is:datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    input date and return the url of data
    """
    date_ymdh=date.strftime('%Y%m%d%H')
    ym=date_ymdh[:6]
    ymd=date_ymdh[:8]
    hour=str(int(int(date_ymdh[8:10])/6)*6)
    if len(hour)==1:
        hour='0'+hour
    url='http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/GOMOFS/MODELS/'\
    +ym+'/nos.gomofs.stations.'+data_str+'.'+ymd+'.t'+hour+'z.nc'
    return url

def get_gomofs(time,lat,lon,depth,mindistance=20):
    """
    the time start in 2016-11-09 00:00:00
    the format time is: datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    lat and lon use decimal degrees
    if the depth is under the water, please must add the marker of '-'
    input time,lat,lon,depth return the temperature of specify location (or return temperature,nc,rho_index,ocean_time_index)
    return the temperature of specify location
    """
    if time<datetime.datetime.strptime('20161109 000000','%Y%m%d %H%M%S'):
        print('Time out of range!')
        sys.exit()
    try:
        try:
            url=get_gomofs_url(time,'nowcast')
            nc=netCDF4.Dataset(str(url))
        except:
            url=get_gomofs_url(time,'forecast')
            nc=netCDF4.Dataset(str(url))
    except:
        print('please check the website or internet!')
    gomofs_lons=nc.variables['lon_rho'][:]
    gomofs_lats=nc.variables['lat_rho'][:]
    gomofs_temp=nc.variables['temp'][:]
    gomofs_time=nc.variables['ocean_time'][:]
    gomofs_h=nc.variables['h'][:]
    gomofs_rho=nc.variables['s_rho'][:]
    
    # find the nearest time of input time
    diff_time=datetime.timedelta(days=10) #specify the initial value of diff_time
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
    if distance0>mindistance:
        print('The location is out of range!')
        sys.exit()
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
