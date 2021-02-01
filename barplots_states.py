import pickle
import numpy as np
import matplotlib.pyplot as plt
from sron_toolbox import * 
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import rgb2hex, Normalize
import datetime as dt
import math
from pylab import *
import json
import collections as col
from scipy.optimize import minimize
from Regions import *
import ast
import os
from matplotlib import gridspec
import matplotlib.cm as cm

#from pklload2 import *

####### Functions

def readgfas(start_mjd, end_mjd, region, polygon):
    "function to read json files containing GFAS data and filter it for the given province or igp3."

    y0, m0, d0 = mjd2date(start_mjd).year, mjd2date(start_mjd).month, mjd2date(start_mjd).day
    y1, m1, d1 = mjd2date(end_mjd-1).year, mjd2date(end_mjd-1).month, mjd2date(end_mjd-1).day

    mjd = start_mjd
    m=m0
    y=y0
    gfasmonthly = []
    months = []
    monthlystd = []
    while y<y1+1:
        while m<m1+1:
            with open('gfas_india2/gfas_'+str(y)+str(m)+'.json', 'r') as f:
                gfasd = json.load(f)
                gfasmonthly.append(np.array(gfasd['data']))
            f.close()
            months.append((y,m))
            m+=1
        y+=1
        m=1

    gfasmonthly[0], gfasmonthly[-1] = gfasmonthly[0][d0-1:], gfasmonthly[-1][:d1]

    lats = np.flipud(np.arange(region.get('latmin'), region.get('latmax'), 0.1))
    lons = np.arange(region.get('longmin'), region.get('longmax'), 0.1)
    mask = np.zeros(np.shape(gfasmonthly[0][0]))
    gfasmonthlymask = []
    gfasmonthlysum = []
    gfasdaily = []
    for i in range(len(lats)):
        for j in range(len(lons)):
            if point_in_poly(lons[j], lats[i], polygon['coordinates']) == False:
                mask[i][j] = 1
    for month in range(len(gfasmonthly)):
        monthly = []
        for day in range(gfasmonthly[month].shape[0]):
            gfasmask = np.ma.masked_array(gfasmonthly[month][day], mask=mask)
            gfasmonthlymask.append(gfasmask)
            gfasdaily.append(np.sum(gfasmask)/10**9)
            monthly.append(np.sum(gfasmask)/10**9)
        monthlystd.append(np.sqrt(np.var(monthly)*len(monthly)))
        gfasmonthlysum.append(sum(monthly))
    mjdlist = range(start_mjd, end_mjd)
    gfasdates = []
    for i in mjdlist:
        gfasdates.append(dt.date(mjd2utc(i)['year'], mjd2utc(i)['month'], mjd2utc(i)['day']))
    return mjdlist, gfasdates, gfasdaily, months, gfasmonthly, monthlystd

def readgfed(start_mjd, end_mjd, region, polygon):
    y0, m0, d0 = mjd2date(start_mjd).year, mjd2date(start_mjd).month, mjd2date(start_mjd).day
    y1, m1, d1 = mjd2date(end_mjd-1).year, mjd2date(end_mjd-1).month, mjd2date(end_mjd-1).day
    mjd = start_mjd
    gfedmonthly = []
    monthlystd = []
    months = []
    m=m0
    y=y0
    while y<y1+1:
        while m<m1+1:
            with open('gfed_india2/gfed__'+str(y)+str(m)+'.json', 'r') as f:
                gfedd = json.load(f)
                gfedmonthly.append(np.array(gfedd['data']))
            f.close()
            months.append((y,m))
            m+=1
        y+=1
        m=1
    gfedmonthly[0], gfedmonthly[-1] = gfedmonthly[0][d0-1:], gfedmonthly[-1][:d1]
    lats = np.flipud(np.arange(region.get('latmin'), region.get('latmax'), 0.25))
    lons = np.arange(region.get('longmin'), region.get('longmax'), 0.25)
    mask = np.zeros(np.shape(gfedmonthly[0][0]))
    gfedmonthlymask = []
    gfedmonthlysum = []
    gfeddaily = []
    for i in range(len(lats)):
        for j in range(len(lons)):
            if point_in_poly(lons[j], lats[i], polygon['coordinates']) == False:
                mask[i][j] = 1
#    imshow(mask)
#    show()
    for i in range(len(gfedmonthly)):
        monthly = []
        for j in range(gfedmonthly[i].shape[0]):
            gfedmask = np.ma.masked_array(gfedmonthly[i][j], mask=mask)
            gfedmonthlymask.append(gfedmask)
            monthly.append(np.sum(gfedmask)/10**12)
            gfeddaily.append(np.sum(gfedmask)/10**12)#*region.get('area'))
        monthlystd.append(np.sqrt(np.var(monthly)*len(monthly)))
        gfedmonthlysum.append(sum(monthly))
    
    mjdlist = range(start_mjd, end_mjd)
    gfeddates = []
    for i in mjdlist:
        gfeddates.append(dt.date(mjd2utc(i).get('year'), mjd2utc(i).get('month'), mjd2utc(i).get('day')))
    return mjdlist, gfeddates, gfeddaily, months, gfedmonthlysum, monthlystd

def readviirs(start_mjd, end_mjd, polygon):
    with open('viirs.json', 'r') as f:
        viirsd = json.load(f)
    try:
        i0, i1 = viirsd['mjd'].index(start_mjd), viirsd['mjd'].index(end_mjd)
    except:
        i0, i1 = viirsd['mjd'].index(min(viirsd['mjd'])), viirsd['mjd'].index(end_mjd)
    mjdlist = viirsd['mjd'][i0:i1]
    frpdaily = np.array(viirsd['frp'][i0:i1])
    coords = np.array(viirsd['coords'][i0:i1])
#    for j in range(len(coords)):
#        for i in range(len(coords[j])):
#            plt.scatter(coords[j][i][0], coords[j][i][1])
#    plt.show()
    frpmean = []
#    mask = np.copy(frpdaily)
    viirsregion = []
#    print frpdaily
    for i in range(len(frpdaily)):
        frpday = []
        for j in range(len(frpdaily[i])):
            if point_in_poly(coords[i][j][0], coords[i][j][1], polygon['coordinates']) == True:
                frpday.append(frpdaily[i][j])
#                plt.scatter(coords[i][j][0], coords[i][j][1])
        frpmean.append(sum(frpday)/1000000.)
        del frpday
#    plt.show()
    viirsdates = []
    for i in mjdlist:
        viirsdates.append(dt.date(mjd2utc(i).get('year'), mjd2utc(i).get('month'), mjd2utc(i).get('day')))

    return mjdlist, viirsdates, frpmean

def readsage(start_mjd, end_mjd, region, polygon):
    y0, m0, d0 = mjd2date(start_mjd).year, mjd2date(start_mjd).month, mjd2date(start_mjd).day
    y1, m1, d1 = mjd2date(end_mjd-1).year, mjd2date(end_mjd-1).month, mjd2date(end_mjd-1).day
    mjd = start_mjd
    sagemonthly = []
    sagemonthlysum = []
    monthlystd = []
    mjdlist = []
    months = []
    m=m0
    y=y0
    while y<y1+1:
        while m<m1+1:
            try:
                with open('sage_india2/sage_'+str(y)+str(m)+'.json', 'r') as f:
                    saged = json.load(f)
#                   saged = ast.literal_eval(json.dumps(saged))
                    sagemonthly.append((saged))
                f.close()
                mjdlist+=sagemonthly[-1]['mjdlist']
                months.append((y,m))
                m+=1
#                y, m = mjd2date(mjd).year, mjd2date(mjd).month
            except:
                m+=1
        y+=1
        m=1
#               y, m = mjd2date(mjd).year, mjd2date(mjd).month
#       print 'SAGE ', mjdlist == range(start_mjd, end_mjd)
    lats = np.array(sagemonthly[0]['lats'])
    lons = np.array(sagemonthly[0]['lons'])
    for i in range(len(sagemonthly)):
        sagemonthly[i] = np.array(sagemonthly[i]['data'])
    sagemonthly[0], sagemonthly[-1] = sagemonthly[0][d0-1:], sagemonthly[-1][:d1]

    mask = np.zeros(np.shape(sagemonthly[0][0]))
    sagemonthlymask = []
    sagedaily = []
    
    for i in range(len(lats)):
        for j in range(len(lons)):
            if point_in_poly(lons[j], lats[i], polygon['coordinates']) == False:
                mask[i][j] = 1

    for month in range(len(sagemonthly)):
        monthly = []
        for day in range(sagemonthly[month].shape[0]):
            sagemask = np.ma.masked_array(sagemonthly[month][day], mask=mask)
            sagemonthlymask.append(sagemask)
            sagedaily.append(np.sum(sagemask)/10**9)
            monthly.append(np.sum(sagemask)/10**9)
        monthlystd.append(np.sqrt(np.var(monthly)*len(monthly)))
        sagemonthlysum.append(sum(monthly))

    sagedates = []
    for i in mjdlist:
        sagedates.append(dt.date(mjd2utc(i)['year'], mjd2utc(i)['month'], mjd2utc(i)['day']))
    print 'SAGE ',len(sagedates), len(sagedaily)
    return mjdlist, sagedates, sagedaily, months, sagemonthlysum, monthlystd

def mjd2date(mjd):
    "converts mjd to datetime object"
    return dt.date(mjd2utc(mjd).get('year'), mjd2utc(mjd).get('month'), mjd2utc(mjd).get('day'))
    
def normalize(array):
    mi = np.min(array)
    array = array-np.min(array)
    array = array/np.max(array)
    return array

def runningmean(x, N):
    avg = np.convolve(x, np.ones((N,))/N, mode='valid')
    return avg 

#function to sum total emissions in a given season
def total(season, data, mjdlist):
    start_mjd, end_mjd = season['start_mjd'], season['end_mjd']
    try:
        i0, i1 = mjdlist.index(start_mjd), mjdlist.index(end_mjd)
    except:
        mjddiff = [abs(mjdlist[i]-start_mjd) for i in range(len(mjdlist))]
        i0 = mjdlist[mjddiff.index(min(mjddiff))]  
        mjddiff = [abs(mjdlist[i]-end_mjd) for i in range(len(mjdlist))]
        i1 = mjdlist[mjddiff.index(min(mjddiff))]
    try:
        totalco = sum(data[i0:i1])
        var = np.var(data[i0:i1])*(i1-i0)
        std = np.sqrt(var)
    except:
        totalco = 0
        std = 0
    return totalco, std

def bargraph(seasons, ylabels, ylist, stdlist):
"function which makes barplot of total emissions in different seasons.
Requires as input:
    seasons = list of seasons dictionaries (defined below)
    ylabels = labels of emissions product names 
    ylist = sums of emissions
    stdlist = standard deviations on sums"

    n = float(len(ylist))
    x = np.arange(len(seasons))  # the label locations
    width = 0.3  # the width of the bars
    fig, ax = plt.subplots()
    bars = []
    ax.bar(x-width*1.5, ylist[0], width, yerr=stdlist[0], color='r', label=ylabels[0])
    ax.bar(x-width*0.5, ylist[1], width, yerr=stdlist[1], color='b', label=ylabels[1])
    ax.bar(x+width*0.5, ylist[2], width, yerr=stdlist[2], color='k', label=ylabels[2])
    ax.set_ylabel('Fires Emission (GFAS, GFED, SAGE), Tg of CO')
    ax.set_title('Emissions Products, '+str(polygon['name']))
    ax.set_xticks(x)
    ax.set_xticklabels(['2017 SOND','2018 MAMJ','2018 SOND','2019 MAMJ','2019 SOND'])
#    ax2 = ax.twinx()
#    ax2.bar(x+width, ylist[3], width, yerr=stdlist[3], color='m', label=ylabels[3])
#    ax2.set_ylabel('Fires Emission (VIIRS), FRP in TW')
    ax.legend(loc='upper right')
#    ax2.legend(loc='upper right')
#    for bar in bars:
#        autolabel(bar)
    fig.tight_layout()
    plt.show()
#    plt.savefig('barplot_'+str(polygon), dpi=150)
#    plt.close()

def bar_states(polygonlist, seasons, statenames):
    gfassumlist = []
    gfasstdlist = []
    gfedsumlist = []
    gfedstdlist = []
    sagesumlist = []
    sagestdlist = []
    for k in range(len(polygonlist)):
        gfasmjd, gfasdates, gfas, gfasmonths, gfassum, gfasmonthlystd = readgfas(start_mjd, end_mjd, region, polygonlist[k])
        gfedmjd, gfeddates, gfed, gfedmonths, gfedsum, gfedmonthlystd = readgfed(start_mjd, end_mjd, region, polygonlist[k])
        sagemjd, sagedates, sage, sagemonths, sagesum, sagemonthlystd = readsage(start_mjd, end_mjd, region, polygonlist[k])
        gfassum = []
        gfasstd = []
        gfedsum = []
        gfedstd = []
        sagesum = []
        sagestd = []
        for j in range(len(seasons)):
            gfassum.append(total(seasons[j], gfas, gfasmjd)[0])
            gfasstd.append(total(seasons[j], gfas, gfasmjd)[1])
            gfedsum.append(total(seasons[j], gfed, gfedmjd)[0])
            gfedstd.append(total(seasons[j], gfed, gfedmjd)[1])
            sagesum.append(total(seasons[j], sage, sagemjd)[0])
            sagestd.append(total(seasons[j], sage, sagemjd)[1])
        gfassumlist.append(gfassum)
        gfasstdlist.append(gfasstd)
        gfedsumlist.append(gfedsum)
        gfedstdlist.append(gfedstd)
        sagesumlist.append(sagesum)
        sagestdlist.append(sagestd)
    gfassumlist = np.array(gfassumlist)
    gfasstdlist = np.array(gfasstdlist)
    gfedsumlist = np.array(gfedsumlist)
    gfedstdlist = np.array(gfedstdlist)
    sagesumlist = np.array(sagesumlist)
    sagestdlist = np.array(sagestdlist)
    x = np.arange(len(polygonlist))
    width = 0.1  # the width of the bars
    fig, (ax1,ax2,ax3) = plt.subplots(3,1)
    bars = []
    seasonnames = ['2017 SOND', '2018 MAMJ', '2018 SOND', '2019 MAMJ', '2019 SOND']
    colours = ['b', 'g', 'y', 'r', 'm']
    for i in range(5):
        r, g, b = 0.1, 0.1, 0.1
        ax1.bar(x-0.25+width*i, gfassumlist[:,i], width, yerr=gfasstdlist[:,i],ecolor='k', color=(0.1,0.3,0.8,0.15+float(i)/7.), label=seasonnames[i])
        ax2.bar(x-0.25+width*i, gfedsumlist[:,i], width, yerr=gfedstdlist[:,i],ecolor='k',color=(0.3,0.8,0.4,0.15+float(i)/7.), label=seasonnames[i])
        ax3.bar(x-0.25+width*i, sagesumlist[:,i], width, yerr=sagestdlist[:,i],ecolor='k',color=(0.9,0.3,0.2,0.15+float(i)/7.), label=seasonnames[i])

    ax1.set_ylabel('GFAS, Tg of CO')
    ax2.set_ylabel('GFED, Tg of CO')
    ax3.set_ylabel('SAGE, Tg of CO')
    ax1.set_ylim(0,0.7)
    ax2.set_ylim(0,0.7)
    ax3.set_ylim(0,0.7)
    ax1.set_title('Emissions Products')
    ax1.set_xticks(x)
    ax2.set_xticks(x)
    ax3.set_xticks(x)
#    plt.yticks(np.arange(0, 0.7, 0.05))
    ax1.set_xticklabels(statenames)
    ax2.set_xticklabels(statenames)
    ax3.set_xticklabels(statenames)
    ax1.legend(loc='upper right')
    ax2.legend(loc='upper right')
    ax3.legend(loc='upper right')
    plt.show()
    
def bar_months(metric):
    fig, (ax1,ax2) = plt.subplots(2,1)
    sumlist, sumlist2 = [], []
    stdlist, stdlist2 = [], []
    for k in range(len(polygonlist)):
        gfedmjd, gfeddates, gfed, gfedmonths, gfedsum, gfedmonthlystd = readgfed(start_mjd, end_mjd, region, polygonlist[k])
        gfed_dict = {'mjd': gfedmjd, 'dates': gfeddates, 'daily':gfed, 'months':gfedmonths, 'monthly':gfedsum, 'monthlystd':gfedmonthlystd}
        sumlist.append(gfed_dict['monthly'][6:10])
        sumlist2.append(gfed_dict['monthly'][18:22])
        stdlist.append(gfed_dict['monthlystd'][6:10])
        stdlist2.append(gfed_dict['monthlystd'][18:22])

    #months = metric['months']
    #totals = metric['monthly']
    #stdlist = metric['monthlystd']
    x = np.arange(len(sumlist[0]))
    print len(polygon), len(sumlist[0])
    sumlist = np.array(sumlist)
    stdlist = np.array(stdlist)
    labels = ['Punjab(Pak)', 'Punjab(Ind)', 'Haryana', 'Uttar Pradesh', 'Bihar']
    monthlist = ['March','April','May','June']
    bars = []
    colours = ['b', 'g', 'y', 'r', 'm']
    width = 0.1
    for i in range(len(sumlist)):
        ax1.bar(x-0.25+width*i, sumlist[i,:], width, yerr=stdlist[i,:], ecolor='k',color=colours[i], label=labels[i])
    #ax.bar(x-width*0.5, ylist[1], width, yerr=stdlist[1], color='b', label=ylabels[1])
    #print [len(ylist[i]) for i in range(len(ylist))]
    #ax.bar(x+width*0.5, ylist[2], width, yerr=stdlist[2], color='k', label=ylabels[2])
    #ax.bar(x+width, ylist[3], width, yerr=stdlist[3], color='m', label=ylabels[3])
#    for i in range(int(n)):
#        bar = ax.bar(x-width/n, ylist[i], width, label=ylabels[i])
#        bars.append(bar)
    ax1.set_ylabel('Fires Emission (GFED), Tg of CO')
    ax1.set_title('Emission 2018 Pre-Monsoon')
    ax1.set_xticks(x)
#    plt.yticks(np.arange(0, 0.7, 0.05))
    ax1.set_xticklabels(monthlist)
#    ax2 = ax.twinx()
#    ax2.bar(x+width, ylist[3], width, yerr=stdlist[3], color='m', label=ylabels[3])
#    ax2.set_ylabel('Fires Emission (VIIRS), FRP in TW')
    ax1.legend(loc='upper right')

    #months = metric['months']
    #totals = metric['monthly']
    #stdlist = metric['monthlystd']
    x = np.arange(len(sumlist2[0]))
    sumlist = np.array(sumlist2)
    stdlist = np.array(stdlist2)
    bars = []
    for i in range(len(sumlist)):
        ax2.bar(x-0.25+width*i, sumlist[i,:], width, yerr=stdlist[i,:], ecolor='k',color=colours[i], label=labels[i])
    ax2.set_ylabel('Fires Emission (GFED), Tg of CO')
    ax2.set_title('Emission 2019 Pre-Monsoon')
    ax2.set_xticks(x)
#    plt.yticks(np.arange(0, 0.7, 0.05))
    ax2.set_xticklabels(monthlist)
#    ax2.legend(loc='upper right')
#    for bar in bars:
#        autolabel(bar)
    
    fig.tight_layout()
    plt.show()
#    plt.savefig('barplot_'+str(polygon), dpi=150)
#    plt.close()
    
def seasonfilter(season, data, dates):
    start_mjd = season['start_mjd']
    end_mjd = season['end_mjd']+1
    start, end = mjd2date(start_mjd), mjd2date(end_mjd)
    while start not in dates:
        start_mjd+=1
        start = mjd2date(start_mjd)
    while end not in dates:
        end_mjd+=1
        end = mjd2date(end_mjd)
#    print start, end
    i0, i1 = dates.index(start), dates.index(end)
    data = data[i0:i1]
    dates = dates[i0:i1]
    return data, dates

############ Main
#Set start and end dates, region and polygon(province)
y0, m0, d0 = 2017, 9, 1
y1, m1, d1 = 2019, 12, 31
region = india2
polygon = igp3 #igp3, punjab_pak, punjab_ind, haryana, up, or bihar
dirname = 'barplots01/'

start_mjd = int(utc2mjd(y0, m0, d0, 0, 0, 0))
end_mjd = int(utc2mjd(y1, m1, d1, 0, 0, 0))+1

#dictionaries of seasons
sond2017 = {'name':'2017 SOND', 'start_mjd':int(utc2mjd(2017,9,1,0,0,0)), 'end_mjd':int(utc2mjd(2017,12,31,0,0,0))}
mamj2018 = {'name':'2018 MAMJ', 'start_mjd':int(utc2mjd(2018,3,1,0,0,0)), 'end_mjd':int(utc2mjd(2018,6,30,0,0,0))}
sond2018 = {'name':'2018 SOND', 'start_mjd':int(utc2mjd(2018,9,1,0,0,0)), 'end_mjd':int(utc2mjd(2018,12,31,0,0,0))}
mamj2019 = {'name':'2019 MAMJ', 'start_mjd':int(utc2mjd(2019,3,1,0,0,0)), 'end_mjd':int(utc2mjd(2019,6,30,0,0,0))}
sond2019 = {'name':'2019 SOND', 'start_mjd':int(utc2mjd(2019,9,1,0,0,0)), 'end_mjd':int(utc2mjd(2019,12,31,0,0,0))}

seasons = [sond2017,mamj2018,sond2018,mamj2019,sond2019]
polygonlist = [punjab_ind, haryana, up, bihar] #imported from Regions.py
statenames = ['Punjab (Ind)', 'Haryana', 'Uttar Pradesh', 'Bihar']

gfasmjd, gfasdates, gfas, gfasmonths, gfassum, gfasmonthlystd = readgfas(start_mjd, end_mjd, region, polygon)
gfedmjd, gfeddates, gfed, gfedmonths, gfedsum, gfedmonthlystd = readgfed(start_mjd, end_mjd, region, polygon)
sagemjd, sagedates, sage, sagemonths, sagesum, sagemonthlystd = readsage(start_mjd, end_mjd, region, polygon)

####################Fires bar graphs

#lists of sums of seasonal emissions and lists of standard deviations of each
gfassum, gfedsum, sagesum, viirssum = [], [], [], []
gfasstd, gfedstd, sagestd, viirsstd = [], [], [], []
for j in range(len(seasons)):
    gfassum.append(total(seasons[j], gfas, gfasmjd)[0])
    gfasstd.append(total(seasons[j], gfas, gfasmjd)[1])
    gfedsum.append(total(seasons[j], gfed, gfedmjd)[0])
    gfedstd.append(total(seasons[j], gfed, gfedmjd)[1])
    sagesum.append(total(seasons[j], sage, sagemjd)[0])
    sagestd.append(total(seasons[j], sage, sagemjd)[1])
#    viirssum.append(total(seasons.keys()[j], viirs, viirsmjd)[0])
#    viirsstd.append(total(seasons.keys()[j], viirs, viirsmjd)[1])
ylist = [gfassum, gfedsum, sagesum, viirssum]
stdlist = [gfasstd, gfedstd, sagestd, viirsstd]

bargraph(seasons, ['GFAS', 'GFED', 'SAGE'], ylist, stdlist)

bar_states(polygonlist, seasons, statenames)


##################PHASE SHIFT
#gfasphase = []
#gfedphase = []
#viirsphase = []
#for season in seasons.keys():
#    gfasr = []
#    gfedr = []
#    viirsr = []
#    for i in range(15):
#        gfasopt = correlation(i, seasons[season], gfas, gfasdates, 'GFAS', codaily, datelist, 3)
#        gfedopt = correlation(i, seasons[season], gfed, gfeddates, 'GFED', codaily, datelist, 3)
#        viirsopt = correlation(i, seasons[season], viirs, viirsdates, 'VIIRS', codaily, datelist, 3)      
#        gfasr.append(gfasopt[0])
#        gfedr.append(gfedopt[0])
#        viirsr.append(viirsopt[0])
#    gfasphase.append([season, max(gfasr), gfasr.index(max(gfasr))])
#    gfedphase.append([season, max(gfedr), gfedr.index(max(gfedr))])
#    viirsphase.append([season, max(viirsr), viirsr.index(max(viirsr))])
#    del gfasr, gfedr, viirsr
#        
#print gfasphase, gfedphase, viirsphase
#delays = [0,12,13,3]
#gfasdatai = []
#gfeddatai = []
#viirsdatai = []
#codatai = []
#for i in range(len(delays)):
#    gfasi = correlation(delays[i], seasons[seasons.keys()[i]], gfas, gfasdates, 'GFAS', codaily, datelist, 3)
#    gfedi = correlation(delays[i], seasons[seasons.keys()[i]], gfed, gfeddates, 'GFED', codaily, datelist, 3)
#    viirsi = correlation(delays[i], seasons[seasons.keys()[i]], viirs, viirsdates, 'VIIRS', codaily, datelist, 3)
#    gfasdatai.append(gfasi[1])
#    gfeddatai.append(gfedi[1])
#    viirsdatai.append(viirsi[1])
#    codatai.append(gfasi[2])
#    print seasons.keys()[i]
#    plt.plot(gfasi[1], label='GFAS')
#    plt.plot(gfedi[1])
#    plt.plot(viirsi[1])
#    plt.plot(gfasi[2], label='CO Column')
#    plt.title(str(seasons.keys()[i])+', delay = '+str(delays[i])+' days')
#    plt.ylabel('Normalized Emission and CO Concentrations')
#    plt.legend(loc='upper right')
#    plt.show()


#plt.plot(1,1)
#plt.show()



