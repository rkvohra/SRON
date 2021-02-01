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
    "function to calculate area of patches on surface of the globe based on lat, lon"
    
    deg2rad = np.pi/180.
    radius  = 6.371e6
    dxx=360.0/im*deg2rad
    dyy=180.0/jm*deg2rad
    lat=np.arange(-90*deg2rad,90*deg2rad,dyy)
    dxy=dxx*(sin(lat+dyy)-sin(lat))*radius**2
    area=np.resize(np.repeat(dxy,im,axis=0) ,[jm,im])
    if not silent:
        print 'total area of field = ',sum(area.flat)
        print 'total earth area    = ',4*pi*radius**2
    return area

def gfasseries(y, m, region, area):
    "function to convert GFAS data to monthly json files. returns list of dates as mjd (modified julian date) and list of GFAS emissions over given region on those dates"
    rowmin = (90-region['latmax'])*10
    rowmax = (90-region['latmin'])*10
    colmin = region['longmin']*10
    colmax = region['longmax']*10
    day_sec = float(60.*60.*24.)
    if m < 10:
        data = (day_sec*((area*nc_read_vars('/deos/ivarvdv/GFAS_data/gfas_'+str(y)+'0'+str(m)+'.nc', ['cofire'])['cofire'])[:,rowmin:rowmax,colmin:colmax])).tolist()
    else:
        data = (day_sec*((area*nc_read_vars('/deos/ivarvdv/GFAS_data/gfas_'+str(y)+str(m)+'.nc', ['cofire'])['cofire'])[:,rowmin:rowmax,colmin:colmax])).tolist()
    day0, day1 = 1, len(data)
    mjd0, mjd1 = int(utc2mjd(y,m,day0,0,0,0)), int(utc2mjd(y,m,day1,0,0,0))
    mjdlist = range(mjd0, mjd1+1)
    return mjdlist, data

start_mjd = 6575
end_mjd = 6605
area = globarea(3600,1800)
region = {'latmax':india2['latmax'], 'latmin':india2['latmin'], 'longmax':india2['longmax'],  'longmin':india2['longmin']}

#beginning month and year
m = 9
y = 2017

#loop through all remaining months in year, save json files
while m < 13:
    dates, data = gfasseries(y, m, region, area)
    gfasd = {'dates':dates, 'data':data}
    with open('gfas_india2_test/gfas_'+str(y)+str(m)+'.json', 'w') as f:
        json.dump(gfasd, f)
    f.close()
    del dates, data, gfasd
    m+=1
    
