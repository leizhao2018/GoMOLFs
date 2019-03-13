#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 15:37:42 2019
get the data from GoMOFs
update function get_gomofs in download data(correct the part name "start upload data" to donload data)
add a function(get_gomofs_url_forcast(date,forcastdate=1))
march 13:update function get_gomofs in download data
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
import math
import time


def get_gomofs_url(date):
    """
    the format of date is:datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    input date and return the url of data
    """
#    print('start calculate the url!') 
    date=date+datetime.timedelta(hours=4.5)
    date_str=date.strftime('%Y%m%d%H%M%S')
    hours=int(date_str[8:10])+int(date_str[10:12])/60.+int(date_str[12:14])/3600.
    tn=int(math.floor((hours)/6.0)*6)  ## for examole: t12z the number is 12
    if len(str(tn))==1:
        tstr='t0'+str(tn)+'z'   # tstr in url represent hour string :t00z
    else:
        tstr='t'+str(tn)+'z'
    if round((hours)/3.0-1.5,0)==tn/3:
        nstr='n006'       # nstr in url represent nowcast string: n003 or n006
    else:
        nstr='n003'
    url='http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/GOMOFS/MODELS/'\
    +date_str[:6]+'/nos.gomofs.fields.'+nstr+'.'+date_str[:8]+'.'+tstr+'.nc'
    return url

def get_gomofs_url_forcast(date,forcastdate=True):
    """
    the date use to choose file
    the format of date is GMT TIME:datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    forcastdate like date or True
    input date and return the url of data
    """
    if forcastdate==True:  #if forcastdate is True: default the forcast date equal to the time of choose file.
        forcastdate=date
    date=date-datetime.timedelta(hours=1.5)  #the parameter of calculate txx(eg:t00,t06 and so on)
    tn=int(math.floor(date.hour/6.0)*6)  #the numer of hours in time index: eg: t12, the number is 12
    ymdh=date.strftime('%Y%m%d%H%M%S')[:10]  #for example:2019011112(YYYYmmddHH)
    if len(str(tn))==1:
        tstr='t0'+str(tn)+'z'  #tstr: for example: t12
    else:
        tstr='t'+str(tn)+'z'
    fnstr=str(3+3*math.floor((forcastdate-datetime.timedelta(hours=1.5+tn)-datetime.datetime.strptime(ymdh[:8],'%Y%m%d')).seconds/3600./3.))#fnstr:the number in forcast index, for example f006 the number is 6
    if len(fnstr)==1:   
        fstr='f00'+fnstr  #fstr: forcast index:for example: f006
    else:
        fstr='f0'+fnstr
    url='http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/GOMOFS/MODELS/'\
    +ymdh[:6]+'/nos.gomofs.fields.'+fstr+'.'+ymdh[:8]+'.'+tstr+'.nc'
    return url

def get_gomofs(date_time,lat,lon,depth,mindistance=20,autocheck=True,forcontours=False):
    """
    the format time(GMT) is: datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    lat and lon use decimal degrees
    if the depth is under the water, please must add the marker of '-'
    input time,lat,lon,depth return the temperature of specify location (or return temperature,nc,rho_index,ocean_time_index)
    the unit is mile of distance
    return the temperature of specify location
    """
    if date_time<datetime.datetime.strptime('2018-07-01 00:00:00','%Y-%m-%d %H:%M:%S'):
        print('Time out of range, time start :2018-07-01 00:00:00z')
        return np.nan
    if date_time>datetime.datetime.now()+datetime.timedelta(days=3): #forecast time under 3 days
        print('forecast time under 3 days')
        return np.nan
   
    #start download data
    forecastdate=date_time  #forecast time equal input date_time
    changefile,filecheck=1,1  #changefile means whether we need to change the file to get data, filecheck means check the file exist or not.

    while(changefile==1):  
        count=1
        while(filecheck==1):  #download the data
            try:
                if forecastdate==date_time:   #the forcastdate is input date_time, if the date_time changed yet,we will use the forecast data
                    url=get_gomofs_url(date_time)
                    nc=netCDF4.Dataset(str(url))
                    print('download nowcast data.')
                else:
                    url=get_gomofs_url_forcast(date_time,forecastdate)
                    nc=netCDF4.Dataset(str(url))
                    print('download forecast data.')
                filecheck,readcheck=0,1      # if the file is there, filecheck=0,readcheck use to check the file whether read successfully               
            except OSError:
                try:
                    url=get_gomofs_url_forcast(date_time,forecastdate)
                    nc=netCDF4.Dataset(str(url))
                    print('download forecast data.')
                    filecheck,readcheck=0,1  
                except OSError:
                    date_time=date_time-datetime.timedelta(hours=6)
                    if (forecastdate-date_time)>datetime.timedelta(days=3):  #every file only have 3 days data.
                        print('please check the website or file is exist!')
                        return np.nan
                except:
                    return np.nan
            except:
                return np.nan
        print('start read data.')
        while(readcheck==1):  #read data,  if readcheck==1 start loop
            try:
                gomofs_lons=nc.variables['lon_rho'][:]
                gomofs_lats=nc.variables['lat_rho'][:]
                gomofs_temp=nc.variables['temp'][:]
                gomofs_h=nc.variables['h'][:]
                gomofs_rho=nc.variables['s_rho'][:]
                readcheck,changefile=0,0   #if read data successfully, we do not need to loop
                print('end read data.')
            except RuntimeError: 
                count=count+1
                if count>5:
                    if autocheck==True:
                        return np.nan
                    while True:
                        print('it will return nan, if you do not need read again.')
                        cmd = input("whether need read again(y/n)?：")
                        if cmd.lower() == "y":
                            count=1
                            break
                        elif cmd.lower() == "n":
                            cmd2 = input("whether need change file(y/n)?：")
                            if cmd2.lower()=="y":
                                date_time=date_time-datetime.timedelta(hours=6)
                                readcheck,filecheck=0,1
                                break
                            else:
                                print('interrupt read data.')
                                return np.nan
                        else:
                            break
                time.sleep(20)   #every time to reread data need stop 20s
                print('the '+str(int(count))+' times to read data.')
            except:
                return np.nan

    #caculate the index of the nearest four points    
    print('start caculate the nearest four points!')
    target_distance=2*zl.dist(lat1=gomofs_lats[0][0],lon1=gomofs_lons[0][0],lat2=gomofs_lats[0][1],lon2=gomofs_lons[0][1])
    eta_rho,xi_rho=zl.find_nd(target=target_distance,lat=lat,lon=lon,lats=gomofs_lats,lons=gomofs_lons)
    
    if zl.dist(lat1=lat,lon1=lon,lat2=gomofs_lats[eta_rho][xi_rho],lon2=gomofs_lons[eta_rho][xi_rho])>mindistance:
        print('THE location is out of range')
        return np.nan
    
    # estimate the bottom depth of point location 
    if eta_rho==0:
        eta_rho=1
    if eta_rho==len(gomofs_lats)-1:
        eta_rho=len(gomofs_lats)-2
    if xi_rho==len(gomofs_lats[0])-1:
        eta_rho=len(gomofs_lats[0])-2
    print('start caculate the bottom depth of point location!') 
    points_h=[[gomofs_lats[eta_rho][xi_rho],gomofs_lons[eta_rho][xi_rho],gomofs_h[eta_rho][xi_rho]],
             [gomofs_lats[eta_rho][xi_rho-1],gomofs_lons[eta_rho][xi_rho-1],gomofs_h[eta_rho][xi_rho-1]],
             [gomofs_lats[eta_rho][xi_rho+1],gomofs_lons[eta_rho][xi_rho+1],gomofs_h[eta_rho][xi_rho+1]],
             [gomofs_lats[eta_rho-1][xi_rho],gomofs_lons[eta_rho-1][xi_rho],gomofs_h[eta_rho-1][xi_rho]],
             [gomofs_lats[eta_rho+1][xi_rho],gomofs_lons[eta_rho+1][xi_rho],gomofs_h[eta_rho+1][xi_rho]]]
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
             [gomofs_lats[eta_rho][xi_rho-1],gomofs_lons[eta_rho][xi_rho-1],gomofs_temp[0][rho_index][eta_rho][xi_rho-1]],
             [gomofs_lats[eta_rho][xi_rho+1],gomofs_lons[eta_rho][xi_rho+1],gomofs_temp[0][rho_index][eta_rho][xi_rho+1]],
             [gomofs_lats[eta_rho-1][xi_rho],gomofs_lons[eta_rho-1][xi_rho],gomofs_temp[0][rho_index][eta_rho-1][xi_rho]],
             [gomofs_lats[eta_rho-1][xi_rho],gomofs_lons[eta_rho-1][xi_rho],gomofs_temp[0][rho_index][eta_rho-1][xi_rho]]]
    temperature=zl.fitting(points_temp,lat,lon)
#    temperature=nc.variables['temp'][0][rho_index][eta_rho][xi_rho]
    # if input depth out of the bottom, print the prompt message
    if depth!='bottom':
        if abs(point_h)<abs(depth):
            print ("the depth is out of the bottom:"+str(point_h))
            return np.nan
    if forcontours==False:
        return temperature
    else:
        return temperature,rho_index,gomofs_temp,gomofs_h,gomofs_lats,gomofs_lons
    

def contours_depth_temp_gomfs(output_path,date_time,lat=41.784712,lon=-69.231081,depth='bottom',addlon=.3,addlat=.3,mod_points='yes',depth_contours_interval=[20, 50,100,150,200,500]):
    """Draw contours and isothermal layers on the map
    notice:
    addlon,addlat: edges around point to include in the zoomed in plot
    depth: if depth=='bottom' print the bottom data.
    the format time is: datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    mod_points:do you want to post model grid nodes,if mod_points='yes', print model grid nodes;if other string, skip 
    """
    #prepare the data
    try:
        temperature,rho_index,gomofs_temp,gomofs_h,gomofs_lats,gomofs_lons=get_gomofs(date_time,lat,lon,depth,forcontours=True)
    except:
        return 0
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
    temp_cs=map.contourf(m_lon,m_lat,gomofs_temp[:][0][rho_index],7)
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
    plt.savefig(output_path+'contour_depth_tem_GoMOFs.png',dpi=300)
#    plt.show()
    return 1
