# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:30:17 2020

@author: rushi
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
from sron_toolbox import * 
import math
import json
from Regions import *
import glob
import netCDF4 as nc

def wrfseries(y, m, region):
    files = glob.glob('/deos/ivarvdv/WRF_OUTPUT_201809-12_withbckgrnd/wrfout_d01_%s-%02d*_s'%(y,m))
    files.sort()
    files = files[:]
    hr = 0
    Clist = []
    C1list = []
    C2list = []
    uu = []
    vv = []
    ulist = []
    vlist = []
    for f in files:
        ncf=nc.Dataset(f)
        U = ncf.variables['U'][:]
        V = ncf.variables['V'][:]
        C = ncf.variables['tracer_ens'][:]
        C1 = ncf.variables['tracer_11'][:]
        C2 = ncf.variables['tracer_12'][:]
        PSFC = ncf.variables['PSFC'][:]
        P = ncf.variables['P'][:]
        Ptest = ncf.variables['PB'][:]
        
        wrf_p = np.zeros((33,99,122))
        wrf_p[0,:,:] = PSFC[0,:,:]
        wrf_p[1:,:,:] = P[0,:,:,:]+Ptest[0,:,:,:]
        delta_p = np.diff(wrf_p,axis=0)
        C = -delta_p[None,:,:,:]*C/(0.028964*9.80665)/1.e6
        C1 = -delta_p[None,:,:,:]*C1/(0.028964*9.80665)/1.e6
        C2 = -delta_p[None,:,:,:]*C2/(0.028964*9.80665)/1.e6
        U = U[:,:,:,:-1]
        V = V[:,:,:-1,:]
        #U = U[:,:,::10,::10]
        #V = V[:,:,::10,::10]
        u_norm = U/np.sqrt(U**2.0+V**2.0)
        v_norm = V/np.sqrt(U**2.0+V**2.0)
        uu.append(u_norm)
        vv.append(v_norm)

        xco2 = C[hr,:,:,:].sum(0)
        C = np.ma.masked_where(xco2[:,:]!=xco2[:,:],C[hr,:,:,:].sum(0))
        C1 = np.ma.masked_where(xco2[:,:]!=xco2[:,:],C1[hr,:,:,:].sum(0))
        C2 = np.ma.masked_where(xco2[:,:]!=xco2[:,:],C2[hr,:,:,:].sum(0))
        C = ((C*0.028964*9.80665)/PSFC[0,:,:])*1.e9
        C1 = ((C1*0.028964*9.80665)/PSFC[0,:,:])*1.e9
        C2 = ((C2*0.028964*9.80665)/PSFC[0,:,:])*1.e9
        C = C.tolist()
        C1 = C1.tolist()
        C2 = C2.tolist()
        XLAT = ncf.variables['XLAT'][:]
        XLONG = ncf.variables['XLONG'][:]
        lon2 = XLONG.tolist()
        lat2 = XLAT.tolist()
        lat = np.ma.masked_where(xco2[:,:]!=xco2[:,:],XLAT[:,:])[:,0].tolist()
        lon = np.ma.masked_where(xco2[:,:]!=xco2[:,:],XLONG[:,:])[0].tolist()
        ncf.close()
        Clist.append(C)
        C1list.append(C1)
        C2list.append(C2)
#        print np.shape(u_norm[hr,:,:,:])
        u_norm = np.mean(u_norm[hr,0:8,:,:],axis=0)
        v_norm = np.mean(v_norm[hr,0:8,:,:],axis=0)
#        u_norm = u_norm[hr,5,:,:]
#        v_norm = v_norm[hr,5,:,:]
        ulist.append(u_norm.tolist())
        vlist.append(v_norm.tolist())
#    print Clist[0][0][0]
#    print lat
#    print np.array(vlist[-1]).shape
    d0, d1 = 1, len(files)
    mjd0, mjd1 = int(utc2mjd(y, m, d0, 0, 0, 0)), int(utc2mjd(y, m, d1, 0, 0, 0))
    mjdlist = range(mjd0, mjd1+1)
    return mjdlist, lat, lon, Clist, C1list, C2list, ulist, vlist, lat2, lon2

y = 2018
m = 9
region = {'latmax':india2['latmax'], 'latmin':india2['latmin'], 'longmax':india2['longmax'],  'longmin':india2['longmin']}
while m < 13:
    dates, lat, lon, C, C1, C2, u, v, wlat, wlon = wrfseries(y, m, region)
    wrfd = {'mjdlist':dates, 'lats':lat, 'lons':lon, 'edgargfed':C, 'gfed':C1, 'edgar':C2, 'u':u, 'v':v, 'wlats':wlat, 'wlons':wlon}
    with open('wrf_india2_gfed/wrf_'+str(y)+str(m)+'.json', 'w') as f:
        json.dump(wrfd, f)
    f.close()
    del dates, lat, lon, C, C1, C2, u, v, wlat, wlon, wrfd
    m+=1
