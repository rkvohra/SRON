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

def viirsseries(region):
    with open('/deos/ivarvdv/VIIRS/fire_archive_V1_108597.json', 'r') as f:
        viirs0 = json.load(f)
    with open('/deos/ivarvdv/VIIRS/fire_nrt_V1_108597.json', 'r') as f:
        viirs0 = viirs0 + json.load(f)
    coords = [(region['longmin'], region['latmin']), (region['longmin'], region['latmax']), (region['longmax'], region['latmax']), (region['longmax'], region['latmin'])]
    viirs = []
    for i in viirs0:
        if point_in_poly(i['longitude'], i['latitude'], coords):
            viirs.append(i)
    frpdaily = []
    coorddaily = []
    mjdlist = []
    for i in viirs:
        mjd = utc2mjd(int(i['acq_date'][0:4]), int(i['acq_date'][5:7]), int(i['acq_date'][8:10]),0,0,0)
        i['mjd'] = mjd
        try:
            if mjdlist[-1] != mjd:
                mjdlist.append(mjd)
            else:
                pass
        except:
            mjdlist.append(mjd)
    for mjd in mjdlist:
        frp_mjd = []
        coord_mjd = []
        for i in viirs:
            if i['mjd'] == mjd:
                frp_mjd.append(i['frp'])
                coord_mjd.append((i['longitude'], i['latitude']))
        coorddaily.append(coord_mjd)
        frpdaily.append((frp_mjd))
        del frp_mjd
    
    return mjdlist, frpdaily, coorddaily

#start_mjd = 6575
#end_mjd = 6605
#area = globarea(3600,1800)
#region = { 'latmax':34, 'latmin':21, 'longmax':92, 'longmin':70    }
region = {'latmax':india2['latmax'], 'latmin':india2['latmin'], 'longmax':india2['longmax'],  'longmin':india2['longmin']}

#while start_mjd < 7277:
#    dates, data = gfasseries(start_mjd, end_mjd, area)
#    gfasd = {'dates':dates, 'data':data}
#    with open('gfas_'+str(start_mjd)+'.json', 'w') as f:
#        json.dump(gfasd, f)
#    f.close()
#    del dates, data, gfasd
#    start_mjd, end_mjd = start_mjd+31, end_mjd+31

mjdlist, frpdaily, coorddaily = viirsseries(region)
viirsd = {'mjd':mjdlist, 'frp':frpdaily, 'coords':coorddaily}
with open('viirs_india2.json', 'w') as f:
    json.dump(viirsd, f)
    f.close()
