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

def gfedseries(y, m, region):
    rowmin = (90-region.get('latmax'))*4-1
    rowmax = (90-region.get('latmin'))*4-1
    colmin = (180+region.get('longmin'))*4-1
    colmax = (180+region.get('longmax'))*4-1

    if m == 1 or m == 3 or m == 5 or m == 7 or m == 8:
        data = pickle.load(open("/deos/ivarvdv/GFED_data/ba6_"+str(y)+'0'+str(m)+".p", "rb"))[0][:,rowmin:rowmax,colmin:colmax].tolist()
    elif m == 2:
        data = pickle.load(open("/deos/ivarvdv/GFED_data/ba6_"+str(y)+'0'+str(m)+".p", "rb"))[0][0:-3][:,rowmin:rowmax,colmin:colmax].tolist()
    elif m == 4 or m == 6 or m == 9:
        data = pickle.load(open("/deos/ivarvdv/GFED_data/ba6_"+str(y)+'0'+str(m)+".p", "rb"))[0][0:-1][:,rowmin:rowmax,colmin:colmax].tolist()
    elif m == 11:
        data = pickle.load(open("/deos/ivarvdv/GFED_data/ba6_"+str(y)+str(m)+".p", "rb"))[0][0:-1][:,rowmin:rowmax,colmin:colmax].tolist()
    elif m == 10 or m == 12:
        data = pickle.load(open("/deos/ivarvdv/GFED_data/ba6_"+str(y)+str(m)+".p", "rb"))[0][:,rowmin:rowmax,colmin:colmax].tolist()

    day0, day1 = 1, len(data)
    mjd0, mjd1 = int(utc2mjd(y, m, day0, 0, 0, 0)), int(utc2mjd(y, m, day1, 0, 0, 0))
    mjdlist = range(mjd0, mjd1+1)
    return mjdlist, data

region = {'latmax':india2['latmax'], 'latmin':india2['latmin'], 'longmax':india2['longmax'],  'longmin':india2['longmin']}

m = 9
y = 2017
while m < 13:
    dates, data = gfedseries(y, m, region)
    gfedd = {'dates':dates, 'data':data}
    with open('gfed_india2/gfed__'+str(y)+str(m)+'.json', 'w') as f:
        json.dump(gfedd, f)
    f.close()
    print len(dates), type(data), np.shape(np.array(data))
    del dates, data, gfedd
    m+=1
    
