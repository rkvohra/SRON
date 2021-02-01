#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 10:45:20 2020

@author: rushilv
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
from sron_toolbox import * 
import math
import json
from Regions import *

def globarea(im=360,jm=180,silent=True):
    deg2rad = np.pi/180.
    radius  = 6.371e6
    dxx=360.0/im*deg2rad
    dyy=180.0/jm*deg2rad
    lat=np.arange(-90*deg2rad,90*deg2rad,dyy)
    dxy=dxx*(np.sin(lat+dyy)-np.sin(lat))*radius**2
    area=np.resize(np.repeat(dxy,im,axis=0) ,[jm,im])
    if not silent:
        print('total area of field = ',sum(area.flat))
        print('total earth area    = ',4*pi*radius**2)
    print 'shape of area ', np.shape(area)
    return area

def sageseries(y, m, region, area):
    rowmin = (90+region.get('latmin'))*4
    rowmax = (90+region.get('latmax'))*4
    colmin = (180+region.get('longmin'))*4
    colmax = (180+region.get('longmax'))*4
    day_sec = float(60.*60.*24.)
    if m < 10:
        data = nc_read_vars('Andreae2/'+str(y)+'/SAGE-IGP+GFEDv4s_'+str(y)+'_0'+str(m)+'_CO.nc', ['lat','lon','CO'])

    else:
        data = nc_read_vars('Andreae2/'+str(y)+'/SAGE-IGP+GFEDv4s_'+str(y)+'_'+str(m)+'_CO.nc', ['lat','lon','CO'])
    print 'shape of CO ', np.shape(data['CO'])
    for i in range(np.shape(data['CO'])[0]):
        data['CO'][i] = day_sec*area*data['CO'][i]
    data['CO'] = data['CO'][:,rowmin:rowmax,colmin:colmax].tolist()
    data['lat'] = data['lat'][rowmin:rowmax].tolist()
    data['lon'] = data['lon'][colmin:colmax].tolist()
    print 'lengths of lon, lat: ', len(data['lon']), len(data['lon'])

    day0, day1 = 1, len(data['CO'])
    mjd0, mjd1 = int(utc2mjd(y, m, day0, 0, 0, 0)), int(utc2mjd(y, m, day1, 0, 0, 0))
    mjdlist = range(mjd0, mjd1+1)
    return mjdlist, data['lat'], data['lon'], data['CO']

start_mjd = 6575
end_mjd = 6605
area = globarea(1440,720)
region = {'latmax':india2['latmax'], 'latmin':india2['latmin'], 'longmax':india2['longmax'], 'longmin':india2['longmin']}

m = 9
y = 2019
while m < 13:
    dates, lats, lons, co  = sageseries(y, m, region, area)
    saged = {'mjdlist':dates, 'lats':lats, 'lons':lons, 'data':co}
    with open('sage_india2/sage_'+str(y)+str(m)+'.json', 'w') as f:
        json.dump(saged, f)
    f.close()
    print len(dates), type(co)
    del dates, lats, lons, co, saged
    m+=1

