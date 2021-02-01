import matplotlib
matplotlib.use('Agg')
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.basemap import Basemap,maskoceans
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import *
from pylab import normpdf
from matplotlib.lines import Line2D
import pickle
import datetime as dt
from distributed_data_structure2 import *
from sron_toolbox import *
from Regions import *
import os

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    import numpy as n
    import scipy.interpolate
    import scipy.ndimage
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)

    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print("[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions.")
        return None
    newdims = n.asarray( newdims, dtype=int )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = n.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = n.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = n.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [n.arange(i, dtype = n.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )
#        print(ndims)
        trorder = [ndims - 1] + list(n.arange( ndims - 1 ))
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = range(n.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (n.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print ("Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported.")
        return None

def date2mjd(day, month, year, hour= 0, minute= 1, second= 1, millisecond=0):
    """docstring for utc2mjd
       Convert date to modified Julian date 2000
       This function is vector-safe.
       Parameters:
           year: the year
           month: the month
           day: the day
           hour: the hour
           minute: the minute
           second: the second
           millisecond: the millisecond (default: 0)
       Return value:
           result: the time/date converted to modified Julian date 2000
    """
    import numpy as np
    import calendar

    t2000 = calendar.timegm((2000, 1, 1, 0, 0, 0)) # seconds in epoch 1970 corresponding to 1 Jan 2000
    if isinstance(year, np.ndarray): # if input is vectors
        is_scalar = False
        if isinstance(millisecond,int): # if millisecond is an int, which means that the optional argument has presumably not been set
            millisecond = np.zeros(year.shape, dtype=np.int) # set millisecond to an array of matching size
        elif any(millisecond != 0) and any(second%1 != 0): # if both millisecond and fractional second are given
            print("Warning: both milliseconds and fractional seconds given! Ignoring fractional seconds.")
            second = np.floor(second).astype(np.int)
    else: # input is scalars
        is_scalar = True
        if millisecond != 0 and second%1 != 0: # if both millisecond and fractional second are given
            print("Warning: both milliseconds and fractional seconds given! Ignoring fractional seconds.")
            second = np.int(np.floor(second)) # cut off fractional seconds
        year = np.array([year]) # convert to arrays
        month = np.array([month])
        day = np.array([day])
        hour = np.array([hour])
        minute = np.array([minute])
        second = np.array([second])
        millisecond = np.array([millisecond])
    result = np.zeros(year.shape, dtype=np.float64) # initialise field for result
    for ind in range(len(year)): # loop over entries of vectors (one entry for scalars)
        t = calendar.timegm((year[ind], month[ind], day[ind], hour[ind], minute[ind], second[ind])) # convert to seconds in epoch 1970
        result[ind] = ( np.float64(t - t2000) + np.float64(millisecond[ind])/1000.0 ) / 86400.0 # convert to (fractional) days since 1 Jan 2000
    if is_scalar: # if input has been scalars
        result = result.item() # convert result back to scalar
    return result


#######################PLOT OPERATIONAL
def rebin( a, newshape ):
    '''Rebin an array to a new shape.
    '''
    assert len(a.shape) == len(newshape)

    slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    coordinates = mgrid[slices]
    indices = coordinates.astype('i')   #choose the biggest smaller integer index
    return a[tuple(indices)]

def removeNANarray(yy= []):
    y = []
    kk=isnan(yy)
    for ii in arange(len(yy)):
        if not kk[ii] :
            y.append(yy[ii])
    return array(y)

def ArrayParam(ar, nam= '' , head= True ):
    dd= sort(removeNANarray(ar.flatten()))
    return ''

def filter_nan(data):
    idx= (np.isnan(data["lat"])==False)

    return idx

def filter_func_tropomi_co(data):
    idx = (data["qa_value"]>0.5) & (np.isnan(data["co_column"])==False)
    #idx=(data["processing_quality_flags"] < 100) & (np.isnan(data["co_column"])==False) & (((data["viewing_azimuth_angle"] > 0) & (data["viewing_zenith_angle"]     > 65)) == False)
    #idx_cloudy=(data["aerosol_height"]<5000) & (data["aerosol_optical_thickness"]>0.5)
    #idx_clear=(data["aerosol_height"]<500) & (data["aerosol_optical_thickness"]<0.5) & (data["surface_altitude"]>0)

    """
    idx_cloud = (data["cloud_fraction"][:,0]>-1.) & (data["cloud_fraction"][:,1]>-1.) & (data["cloud_fraction"][:,2]>-1.) & (data["cloud_fraction"][:,3]>-1.)
    idx_cloud = idx_cloud & (data["cloud_fraction"][:,0]<0.001) & (data["cloud_fraction"][:,1]<0.001) & (data["cloud_fraction"][:,2]<0.001) & (data["cloud_fraction"][:,3]<0.001) 
    idx_noscat = (data["weak_ch4_column"]>0.) & (data["strong_ch4_column"]>0.) &(data["weak_ch4_column"]/data["strong_ch4_column"]>0.94) & \
                 (data["weak_ch4_column"]/data["strong_ch4_column"]<1.06)  
    idx_noscat = idx_noscat & (data["weak_h2o_column"]>0.) & (data["strong_h2o_column"]>0.) &(data["weak_h2o_column"]/data["strong_h2o_column"]>0.92) &   (data["weak_h2o_column"]/data["strong_h2o_column"]<1.08) 
    #    idx_noscat = idx_noscat &  (data["stdv_ch4_ratio"]<0.05) #&  (data["stdv_h2o_ratio"]<0.05)     
    if data["aerosol_optical_thickness"].ndim==1 :


        idx_quality = (np.isnan(data["xch4"])==False) & (data["aerosol_optical_thickness"]<0.06) & (data["chi2"]<100.) & (data["sza"]<70.) & (data["vza"]<60.) & (data["xch4_precision"]<10.0) & (data["surface_altitude_stdv"]<100.)
    if data["aerosol_optical_thickness"].ndim==2:
        print shape(data["xch4"]),  shape(data["aerosol_optical_thickness"]) , shape(data['chi2']) , shape(data['vza']), shape(data['xch4_precision']) , shape(data['surface_altitude_stdv']) 
        
#        (data["chi2"]<100.0) & (data["sza"]<70.) & (data["vza"]<60.) & (data["xch4_precision"]<10.0) & (data["surface_altitude_stdv"]<100.)
        idx_quality = (np.isnan(data["xch4"])==False) & (data["aerosol_optical_thickness"][:,0]<0.3) & (data["chi2"]<100.0) & (data["sza"]<70.) & (data["vza"]<60.) & (data["xch4_precision"]<10.0) & (data["surface_altitude_stdv"]<100.)
    idx_land = (data["glintflag"][:]==0)  

    idx = idx_quality & idx_cloud & idx_noscat & idx_land 
    """

    #    idx = idx_land 
    return  idx #& idx_cloudy

def basemap_plot_mat(fig,lat,lon,v,m,GfedSymbol='s',contourfire=False,LUmap=False,Gfed=False,lognorm=False,vrange="",smoothing=False,colorbar=True, mask_ocean=False, resolution="c", clabel = ' ' ,cmap=sron_cmap("rainbow_WhBr"), alpha = 1.0):
    """docstring for basemap_plot_mat
    
    This function allows to plot data in matrix structure on a basemap map projection

    input:

            lat: latitude grid points (vector with length n)
            lon: longitude grid points (vector with length m)
            v: data to be plotted (matrix mxn)
            m: basemap instance (e.g. configured with the function config_basemap)


            vrange= define the range of the data to be plotted 
            colorbar=False plot no colorbar
            mask_ocean=True mask data points over the ocean
            resolution =  resolution of the ocean contours (c,l,h,...)
            cmap= chose a color map
	    alpha = transparancy of the data on the plot (between 0 and 1)

    output: plot instance
    
    """
    import numpy.ma as ma
    from matplotlib.colors import rgb2hex, Normalize, LogNorm, PowerNorm
    from matplotlib.widgets import Slider,Button
    from squeezed_norm import SqueezedNorm
    from matplotlib.colorbar import ColorbarBase,make_axes

    # transform lon lat -> x y
    lon, lat = np.meshgrid(lon,lat)
    x, y = m(lon,lat)
    
    # create raster plot	    alpha = transparancy of the data on the plot (between 0 and 1)
    Zm = ma.array(v,mask=(np.isnan(v)))
    Array = np.zeros((Zm.shape[0],Zm.shape[1]),)
    lon2 = lon[1:,1:].copy()
    lat2= lat[1:,1:].copy()

    if mask_ocean ==True:
        # mask ocean
        Zm = maskoceans(x[:,:], y[:,:], Zm[:,:].T, inlands=True, resolution=resolution, grid=1.25).T

        pass
                    
    if vrange != "":
        test=np.nan_to_num(Zm)
        #norm = Normalize(vmin=vrange[0], vmax=vrange[1])
        #norm = LogNorm(vmin=vrange[0], vmax=vrange[1])
        #norm = PowerNorm(0.4,vmin=vrange[0], vmax=vrange[1])
        #norm=SqueezedNorm(vmin=vrange[0], vmax=vrange[1], mid=130, s1=1, s2=1)
        if lognorm:
            norm = LogNorm(vmin=vrange[0], vmax=vrange[1])
        else:    
            norm = Normalize(vmin=vrange[0], vmax=vrange[1])
        if LUmap==True:
            Zm2=Zm[0,:,:].copy()
            for i in range(6):
                Zm2[Zm[i,:,:]>0.5] =i+1
     
            Zm2 = np.ma.masked_where(Zm2==0,Zm2)
            try:
                #cs1=m.contourf(x[:-1,:-1],y[:-1,:-1],Zm[i,:,:].T.T,color='red',linewidth=1,zorder=2) # norm=norm,cmap=cmap,zorder=2,linewidths=1,levels=[0,1,3,4,5,6])
                cs1=m.pcolormesh(x[:-1,:-1],y[:-1,:-1],Zm2[:,:],vmin=1,vmax=6,edgecolors='None',shading='gouraud',zorder=2,norm=norm,rasterized=True,cmap=cmap,alpha=alpha)
            except:
                pass
        if Gfed==False: # or LUmap==False:
            cs1=m.pcolormesh(x[:,:],y[:,:], Zm.T,vmin=vrange[0],vmax=vrange[1],edgecolors='None',zorder=2,norm=norm,rasterized=True,cmap=cmap,alpha=alpha)

    if colorbar==True and Gfed==False:
        divider = make_axes_locatable(ax5)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cs1,ax=ax5,cax=cax)
        cb = plt.colorbar(cs1,ax=ax5,cax=cax)
    if colorbar==True and Gfed==True:
        cbaxes1=fig.add_axes([0.1, 0.1, 0.02, 0.8])
        cb = ColorbarBase(cbaxes1,cmap=cmap,norm=norm,orientation="vertical",extend="both" ,label=clabel)   
        
    fig.subplots_adjust(bottom=0.15)
    return cs1, cb

def main(co_gfed4s_days_mean,co_gfed4s_days_mean2,plot='africa', proj="",Step=0.25, Stepx=0.26, conconly=False,PlotData=True,contourfire=False,gfedonly=False,RegridData=False,mask_ocean=False,Species='co',dpi=100,smoothing=False,GfedSymbol='s',Gfed=False,outputdir = "",origpath="",timeperiod = 'daily',st_time=[12,12,2017] ,ed_time= [15,12,2017]  , reReadData= True,PolygonPlot=False,wrf=False,kernel=False):

    # calculate output grid 
    # the grid is chosen bigger to prevent plotting issues
    ulc_lat,ulc_lon,lrc_lat, lrc_lon,dlat, dlon,  step_a= giveRegionBoxCord(plot)    

    ulc_lat, lrc_lat, ulc_lon, lrc_lon = 37.22, 13.59, 62, 94
    #lrc_lat=37
    #ulc_lat=13
    #lrc_lon=94
    #ulc_lon=62
    prj = 'cyl' #if plot=='global' else 'cyl'
    
    step_a= Step  # User defined resolution
    step_b= Stepx
    lat_grid_a=np.arange(lrc_lat,ulc_lat+2*step_a,step_a) 
    lon_grid_a=np.arange(ulc_lon,lrc_lon+2*step_b,step_b)
    
    tit_add= ' '    
    h=12 # h x h database grid
    d=1 # time grid, [days]
    cs=300 # cache size
    path = origpath #'/deos/ivarvdv/cluster/test/Output/output_dds_201807_noprefilter_destripe_latloncorners/'
    this=distributed_data_structure2(h,d,cs,path,file_format="json")

    mjd_start= date2mjd (st_time.day,st_time.month ,st_time.year) 
    mjd_stop= date2mjd (ed_time.day,ed_time.month,ed_time.year)
    
    resolution="i"

    data = {} 

    if kernel==True:
        varrs= ["altitude_levels","co_column_averaging_kernel","co_column","surface_pressure","dry_air_column","lat","lon","mjd"]
    else:
        if Species=='co':
            varrs= ["co_column","surface_pressure","pressure_levels","lat","lon","lat_con","lon_con","mjd"]

    if not  os.path.exists('%s'%(proj)) :        os.mkdir('%s'%(proj))

    if not  os.path.exists('%s/%s'%(proj,outputdir)) :        os.mkdir('%s/%s'%(proj,outputdir))

    if not  os.path.exists('%s/%s/%s'%(proj,outputdir,plot)) :        os.mkdir('%s/%s/%s'%(proj,outputdir,plot))

    if not  os.path.exists('%s/%s/%s/%s'%(proj,outputdir,plot,timeperiod)) :        os.mkdir('%s/%s/%s/%s'%(proj,outputdir,plot,timeperiod))    
 
    def readfromJSONFiles():
        if Species=='co':
            temp=this.get_average(lat_grid_a, lon_grid_a, mjd_start, mjd_stop,varrs,par_sel_err=[""],operator="mean",filter_func=filter_func_tropomi_co)

        print  'reading files between (julian dates)' , mjd_start , mjd_stop    
        for ii  in arange(len(varrs)): 
            data[varrs[ii]]= temp[ii] 
        data['count']= temp[-1] 
   
        print 'sasa',len(data['count']) 
        if len(data['count']) == 0 :
            print '\nno *****data found %s = %s, \n try another dates'%(st_time, ed_time)      
            return 
            print  '\n  !!!!!!!! saving %s/%4.4i-%4.4i.pkl  !!!!!!!\n'%(plot, mjd_start , mjd_stop) 
            pickle.dump(data, open('%s/%s/%s/%s/%4.4i-%4.4i_%s.pkl'%(proj,outputdir,plot,timeperiod,mjd_start , mjd_stop ,Species), 'w'))

        else:
            pickle.dump(data, open('%s/%s/%s/%s/%4.4i-%4.4i_%s.pkl'%(proj,outputdir,plot,timeperiod,mjd_start , mjd_stop ,Species), 'w'))

    if reReadData:     
        readfromJSONFiles()
    
    try:         
        data = pickle.load( open('%s/%s/%s/%s/%4.4i-%4.4i_%s.pkl'%(proj,outputdir,plot,timeperiod,mjd_start , mjd_stop,Species), 'r'))
    except:
        print 'file not found for %s-%s. May there is no data'%(st_time,ed_time)
        return
    if wrf == False:
        if kernel==True:
            ak=data['co_column_averaging_kernel']
        xco=data['co_column']*6.022141e19
        if Species=='co':
            xco = col2vmr(xco,data['surface_pressure'],is_profile=False)*1.e9
            xco = np.ma.masked_where(xco<0.,xco)
        print 'MEAN ', np.ma.mean(xco[np.where(xco>0.)])
        print np.shape(xco)
    try:
        lat2 = data['lat_con']
        lon2 = data['lon_con']
    except:
        print 'error lat_con'
        pass

    tit= tit_add

    day_str= '%2.2i,%2.2i to %2.2i,%2.2i (%4.4i)  '%(st_time.day, st_time.month, ed_time.day, ed_time.month, ed_time.year)

    ArrayParam(xco, 'xco', True)

    day_str= '%2.2i - %2.2i - %4.4i  '%(st_time.month, st_time.day, st_time.year)
    data_a = xco.copy()
    #lat_grid_a3 = lat_grid_a#.copy()
    #lon_grid_a3 = lon_grid_a#.copy()
    #lat_grid_a4 = lat_grid_a4
    #lon_grid_a4 = lon_grid_a4
    #subplot(152)

    #ra,fig=plotSubMap(data_a,'XCO')

    #if plot=="west_USA":
        #fig = plt.figure(figsize=(15, 8))
    #else:
        #fig = plt.figure(figsize=(5, 3))
    #ax = fig.add_subplot(1, 1, 1)
    #ax5 = fig.add_axes([0.6, 0.1, 0.2, 0.8])
    #cbaxes1 = fig.add_axes([0.9, 0.1, 0.03, 0.8]) 
    #cbaxes2 = fig.add_axes([0.1, 0.1, 0.03, 0.8]) 
    cbaxes1=0.
    mea= mean(ma.masked_invalid(data_a)) ; sd= std(ma.masked_invalid(data_a)) * 3
    
    if 'alt' in tit.lower() and 'co' not in tit.lower() : 
        clabel = ('meters')
        range= [0, mea +sd]
    else:
        if Species=='co': 
            clabel ='ppb'
            SpeciesLabel='CO [ppb]'
            gfed_label='Tg CO/grid/day'
            vrange=[0.00001,0.003]
            vrange=[0.00001,0.01]
            range = [-sd , sd] if (mea-sd)* (mea+sd) < 0 else [mea-sd, mea+sd]
            range=[0,200]#60 200
    m5=basemap_config(ax=ax5,prj=prj,show_x=False,show_y=False,latlonbox=[lrc_lat,ulc_lon,ulc_lat,lrc_lon],lw_coast=1.0,resolution=resolution,meridians=np.arange(0.,360.0,dlon), parallels=np.arange(-90,90,dlat), center_latlon=[0,0], fs_latlon= 6, drawstates=True)
    m5.drawmapboundary(fill_color='lightgrey')         
    if not PolygonPlot:
        if gfedonly==True:
            pass
        else: 
            lat_grid_a=lat_grid_a[:-1]
            lon_grid_a=lon_grid_a[:-1]
            #print 'LATS ', lon_grid_a
        co_max = np.max(data_a.transpose())
        basemap_plot_mat(fig,lat_grid_a,lon_grid_a,data_a.transpose(),m5,colorbar=True,clabel=SpeciesLabel,smoothing=smoothing,vrange=range,cmap=sron_cmap("rainbow_WhBr",nan_transparent=True),mask_ocean=mask_ocean,resolution="h",alpha=1)
    
#    m5.drawcountries(linewidth= 0.25)         
    m5.drawcoastlines() 
        #m.drawmapboundary(fill_color='aqua')

    title( tit_add )
    #return range , fig

    suptitle(day_str, fontsize= 15 )
        
    filename_out='%s'%(proj)
    if not  os.path.exists('%s'%(proj)) :        os.mkdir('%s'%(proj))

    filename_out='%s/plots_%s'%(proj,outputdir)
    if not  os.path.exists('%s/plots_%s'%(proj,outputdir)) :        os.mkdir('%s/plots_%s/'%(proj,outputdir))

    filename_out='%s/plots_%s/%s'%(proj,outputdir,plot)
    if not  os.path.exists('%s/plots_%s/%s'%(proj,outputdir,plot)) :        os.mkdir('%s/plots_%s/%s'%(proj,outputdir,plot))

    filename_out='%s/plots_%s/%s/%s'%(proj,outputdir,plot,timeperiod)
    if not  os.path.exists('%s/plots_%s/%s/%s'%(proj,outputdir,plot,timeperiod)) :        os.mkdir('%s/plots_%s/%s/%s'%(proj,outputdir,plot,timeperiod))

    filename_out= '%s/plots_%s/%s/%s/%s_%s.png'%(proj,outputdir,plot,timeperiod,st_time.strftime('%Y-%m-%d'),Species)
  
    if gfedonly==True:
        filename_out= '%s/plots_%s/%s/%s/%s_%s_gfedonly.png'%(proj,outputdir,plot,timeperiod,st_time.strftime('%Y-%m-%d'),Species)

    if contourfire==True:
        filename_out= '%s/plots_%s/%s/%s/%s_%s_contourgfed.png'%(proj,outputdir,plot,timeperiod,st_time.strftime('%Y-%m-%d'),Species)            

    if conconly==True:
        filename_out= '%s/plots_%s/%s/%s/%s_%s_conconly.png'%(proj,outputdir,plot,timeperiod,st_time.strftime('%Y-%m-%d'),Species)

    print  '\n \n ****** ' , filename_out
    return m5, lat_grid_a, lon_grid_a, data_a.transpose()
#######################PLOT OPERATIONAL


atm5 = np.array([0.000000,	22.835938,	424.414063,	1387.546875,	3057.265625,	5564.382813,	8163.375000,	11901.339844,	14898.453125,	17471.839844,	19290.226563,	20361.816406,	20337.863281,	19859.390625,	19031.289063,	18308.433594,	17008.789063,	15508.256836,	13881.331055,	12766.873047,	11116.662109,	9562.682617,	8608.525391,	7311.869141,	6156.074219,	4490.817383,	3381.743652,	2265.431641,	1423.770142,	823.967834,	427.592499,	191.338562,	69.520576,	18.608931])

btm5 = np.array([1.000000,	0.991984,	0.969513,	0.931881,	0.873929,	0.790717,	0.704669,	0.576692,	0.466003,	0.358254,	0.263242,	0.168910,	0.111505,	0.077958,	0.051773,	0.038026,	0.022355,	0.011806,	0.005378,	0.002857,	0.000890,	0.000199,	0.000059,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000])

yr=2018
month=11

#directory to get wrf output from
files = glob.glob('/deos/ivarvdv/WRF_OUTPUT_201809-12_withbckgrnd_sage/wrfout_d01_%s-%02d*_s'%(yr,month))
files.sort()
files=files[:]
print 'FILES ', len(files)
colordefor='red' #'cadetblue'
colorsav= 'blue' #'mediumslateblue'
colorbackground='purple'
fs=16

DCO=[]
DCObg=[]
DCObg_perc=[]
DNO2=[]
DNO2bg=[]
DNO2bg_perc=[]
DCOsav=[]
DCOsavbg=[]
DCOsavbg_perc=[]
DNO2sav=[]
DNO2savbg=[]
DNO2savbg_perc=[]
uu=[]
vv=[]
C4=[]
NO24=[]
trop_copres=[]
trop_coak=[]
trop_no2pres=[]
trop_no2ak=[]
i=0
st_time=dt.datetime(yr,month,1)
ed_time=dt.datetime(yr,month,2)

#######PLOT OPERATIONAL INPUTS
step = 0.243
Stepx = 0.265
tp = 'yearly' #create output folder daily, weekly or monthly
tp2= 1 #average data over number of days
##year=2019
##month=12
##day_start=18
##day_end=31
region = 'india' #'northern_hemisphere' #regions specified in plotting.py
species = 'co' #type species

#OUTPUT DIRECTORY
project="/deos/rushilv/test_india2"

rereaddata =True  #read JSON and write region-specific data to binary pkl files. Set to False if you only want to read the existing pkl files
regriddata= False #test option
plotdata = True
polygon = False
##step=0.25
dpi=150
wrf=False
kernel = False
gfed=False
gfedonly=False
conconly=False
smoothing =False
gfedsymbol='+'
contourfire = False
mask_ocean = False
output = 'operational_gridded_%sx%s_test'%(step,step)
#Set here the directory to get tropomi data from
origpath = '/deos/ivarvdv/cluster/test/Output/output_dds_201801-201812_with_prefilter/'
dlon, dlat = 1,1

for file in files:
    #if '09-30' in file:
    #    continue
    print(file)
    offset=0
    ncf=nc.Dataset(file)
    U=ncf.variables["U"][:]
    V=ncf.variables["V"][:]
    P = ncf.variables["P"][:]
    Ptest = ncf.variables["PB"][:]
    PB=ncf.variables["PH"][:]
    PHB=ncf.variables["PHB"][:]
    PSFC = ncf.variables['PSFC'][:] 
    C = ncf.variables["tracer_ens"][:]
    print 'CDIM = ', np.shape(C)
    C1 = ncf.variables["tracer_11"][:]
    C2 = ncf.variables["tracer_12"][:]
    XLAT=ncf.variables['XLAT'][:]
    XLONG=ncf.variables['XLONG'][:]
    #print ncf.variables['PSFC']
    ncf.close()

    mjd_start= date2mjd (st_time.day,st_time.month ,st_time.year)
    mjd_stop= date2mjd (ed_time.day,ed_time.month,ed_time.year)
    daymonth=[31,28,31,30,31,30,31,31,30,31,30,31]    
    if i+1<daymonth[month-1]:
        st_time=dt.datetime(yr,month,1+i)
        ed_time=dt.datetime(yr,month,2+i)
    if i+1==daymonth[month-1]:
        st_time=dt.datetime(yr,month,1+i)
        ed_time=dt.datetime(yr,month+1,1)
    mjd_start= date2mjd (st_time.day,st_time.month ,st_time.year)
    mjd_stop= date2mjd (ed_time.day,ed_time.month,ed_time.year)

    dx_grid=30000.
    dy_grid=30000.
    g_we=129 #254
    g_sn=97 #249.
    res='i'
    thres=1000
    lat_0=26
    lon_0=78

    tlat1= 13 # 25 #13# -54 
    tlat2= 37 # 34 #37 #2
    tlon1= 62 # 70 #62 #111
    tlon2= 94 # 83 #94#179

    lrc_lat=13 #-50
    ulc_lat=37 #-5
    lrc_lon=62#180
    ulc_lon=94 #105

#    i0,i1,i2,i3 = 50, 90, 30, 80 #for isolating punjab

    wrf_p=np.zeros((33,99,122))
    wrf_p[0,:,:] = PSFC[0,:,:]
    wrf_p[1:,:,:] = P[0,:,:,:]+Ptest[0,:,:,:]
    #print np.shape(Ptest[0,:,:,:])
    delta_p = np.diff(wrf_p, axis=0)
    C = -delta_p[None,:,:,:]*C/(0.028964*9.80665)/1.e6
    C1 = -delta_p[None,:,:,:]*C1/(0.028964*9.80665)/1.e6
    C2 = -delta_p[None,:,:,:]*C2/(0.028964*9.80665)/1.e6
    for hr in [0]:

        U = U[:,:,:,:-1]
        V = V[:,:,:-1,:]
#        U = U[:,:,i0:i1,i2:i3] #[:,:,::2,::2] #for isolating punjab
#        V = V[:,:,i0:i1,i2:i3] #[:,:,::2,::2]
        U = U[:,:,::4,::4]
        V = V[:,:,::4,::4]

        u_norm = U / np.sqrt(U ** 2.0 + V ** 2.0)
        v_norm = V / np.sqrt(U ** 2.0 + V ** 2.0)

        uu.append(u_norm)
        vv.append(v_norm)
#        lon2 = XLONG[i0:i1,i2:i3] #[::2,::2] #for isolating punjab
#        lat2 = XLAT[i0:i1,i2:i3] #[::2,::2]
        lon2 = XLONG[::4,::4]
        lat2 = XLAT[::4,::4]
        wind_mask = np.zeros(np.shape(u_norm))
        for y in range(len(lat2[:,0])):
            for z in range(len(lon2[0])):
                if not point_in_poly(lon2[0][z], lat2[:,0][y], punjab_ind['coordinates']):
#                    u_norm[:,5,y,z]=0
#                    v_norm[:,5,y,z]=0
                    wind_mask[:,5,y,z]=1
#        u_norm = np.ma.masked_array(u_norm,mask=wind_mask)
#        v_norm = np.ma.masked_array(v_norm,mask=wind_mask)
#        print lat2, lon2
        ##xco2 = congrid(xco[21:130,200:366],(99,122),method='nearest',centre=False, minusone=False)
        ##lontropomi2 = congrid(lontropomi[21:130,200:366],(99,122),method='nearest',centre=False, minusone=False)       
        ##lattropomi2= congrid(lattropomi[21:130,200:366],(99,122),method='nearest',centre=False, minusone=False)
#        print 'CDIM = ', np.shape(C)
        xco2=C[hr,:,:,:].sum(0)
        C = np.ma.masked_where(xco2[:,:]!=xco2[:,:],C[hr,:,:,:].sum(0))
        #print 'C = ', C[0]
        C1 = np.ma.masked_where(xco2[:,:]!=xco2[:,:],C1[hr,:,:,:].sum(0))
        C2 = np.ma.masked_where(xco2[:,:]!=xco2[:,:],C2[hr,:,:,:].sum(0))
        #C = np.ma.masked_where(xco2[:,:]!=xco2[:,:],C[hr,:,:,:].mean(0)*1000.)
        #C1 = np.ma.masked_where(xco2[:,:]!=xco2[:,:],C1[hr,:,:,:].mean(0)*1000.)
        #C2 = np.ma.masked_where(xco2[:,:]!=xco2[:,:],C2[hr,:,:,:].mean(0)*1000.)
        C  =  ((C* 0.028964 * 9.81) / PSFC[0,:,:])*1.e9
        C1 = ((C1* 0.028964 * 9.81) / PSFC[0,:,:])*1.e9
        C2 = ((C2* 0.028964 * 9.81) / PSFC[0,:,:])*1.e9
        lon = np.ma.masked_where(xco2[:,:]!=xco2[:,:],XLONG[:,:])
        lat = np.ma.masked_where(xco2[:,:]!=xco2[:,:],XLAT[:,:])

#        lattest = lat[:,0]
#        lontest = lon[0]
#        print np.shape(C), len(lattest), len(lontest)
#        polygon = punjab_ind
#        for r in range(len(lattest)):
#            for m in range(len(lontest)):
#                if not point_in_poly(lontest[m], lattest[r], polygon['coordinates']):
#                    C[r][m] = 0
        #count = 
        #print 'Cavg = ', np.mean(C[np.where(C>0)])
            
        #xco2 =  np.ma.masked_where(xco2[:,:]!=xco2[:,:],xco2)
        ##fig, ((ax1, ax2,ax3,ax4),(ax5, ax6,ax7,ax8),(ax9, ax10,ax11,ax12),(ax13, ax14,ax15,ax16)) = plt.subplots(4,4,figsize=(20,20))
        fig, ((ax1, ax2, ax3), (ax14, ax5, ax6)) = plt.subplots(2,3,figsize=(20,10))
#        fig, ax14 = plt.subplots()
        ##m3 = Basemap(projection='merc',llcrnrlat=tlat1,urcrnrlat=tlat2,\
        ##        llcrnrlon=tlon1,urcrnrlon=tlon2,lat_ts=0,resolution=res,area_thresh=thres,ax=ax3)  

        ##m2 = Basemap(projection='merc',llcrnrlat=tlat1,urcrnrlat=tlat2,\
        ##        llcrnrlon=tlon1,urcrnrlon=tlon2,lat_ts=0,resolution=res,area_thresh=thres,ax=ax2)
 
        m1 = Basemap(projection='cyl',llcrnrlat=tlat1,urcrnrlat=tlat2,\
                llcrnrlon=tlon1,urcrnrlon=tlon2,lat_ts=0,resolution=res,area_thresh=thres,ax=ax1)

        m2 = Basemap(projection='cyl',llcrnrlat=tlat1,urcrnrlat=tlat2,\
                llcrnrlon=tlon1,urcrnrlon=tlon2,lat_ts=0,resolution=res,area_thresh=thres,ax=ax2)

        m3 = Basemap(projection='cyl',llcrnrlat=tlat1,urcrnrlat=tlat2,\
                llcrnrlon=tlon1,urcrnrlon=tlon2,lat_ts=0,resolution=res,area_thresh=thres,ax=ax3)

        m14 = Basemap(projection='cyl',llcrnrlat=tlat1,urcrnrlat=tlat2,\
                llcrnrlon=tlon1,urcrnrlon=tlon2,lat_ts=0,resolution=res,area_thresh=thres,ax=ax14)
        
        m6 = Basemap(projection='cyl',llcrnrlat=tlat1,urcrnrlat=tlat2,\
                llcrnrlon=tlon1,urcrnrlon=tlon2,lat_ts=0,resolution=res,area_thresh=thres,ax=ax6)
        m6.drawmapboundary(fill_color='lightgrey')

        x, y = m1(lon,lat)
        x2, y2 = m14(lon2,lat2)
        part_1=[[232, 236, 251], [221, 216, 239], [209, 193, 225], [195, 168, 209], [181, 143, 194], [167, 120, 180], [155, 98, 167], [140, 78, 153]]
        part_2=[[111, 76, 155], [96, 89, 169], [85, 104, 184], [78, 121, 197], [77, 138, 198], [78, 150, 188], [84, 158, 179], [89, 165, 169], [96, 171, 158], [105, 177, 144], [119, 183, 125], [140, 188, 104], [166, 190, 84], [190, 188, 72], [209, 181, 65], [221, 170, 60], [228, 156, 57], [231, 140, 53], [230, 121, 50], [228, 99, 45], [223, 72, 40], [218, 34, 34]]
        part_3=[[184, 34, 30], [149, 33, 27], [114, 30, 23], [82, 26, 19]]
        cmap_def =part_1 + part_2 + part_3
        cmap_def2 = [[0, 0, 255], [15, 15, 255], [30, 30, 255], [45, 45, 255], [60, 60, 255], [75, 75, 255], [90, 90, 255], [105, 105, 255], [120, 120, 255], [135, 135, 255], [150, 150, 255], [165, 165, 255], [180, 180, 255], [195, 195, 255], [210, 210, 255], [225, 225, 255], [240, 240, 255], [255, 255, 255], [255, 240, 240], [255, 225, 225], [255, 210, 210], [255, 195, 195], [255, 180, 180], [255, 165, 165], [255, 150, 150], [255, 135, 135], [255, 120, 120], [255, 105, 105], [255, 90, 90], [255, 75, 75], [255, 60, 60], [255, 45, 45], [255, 30, 30], [255, 15, 15], [255, 0, 0]]
        result = matplotlib.colors.LinearSegmentedColormap.from_list("rainbow_WhBr", np.array(cmap_def)/256.0)
        result2 = matplotlib.colors.LinearSegmentedColormap.from_list("rainbow_WhBr", np.array(cmap_def2)/256.0)

        aa1=m1.pcolormesh(x,y,C[:,:],vmin=0,vmax=200,cmap=result,ax=ax1)
        aa2=m2.pcolormesh(x,y,C1[:,:],vmin=0,vmax=200,cmap=result,ax=ax2)
        aa3=m3.pcolormesh(x,y,C2[:,:],vmin=0,vmax=200,cmap=result,ax=ax3)
        wrf_lat = XLAT[:,0]
        wrf_lon = XLONG[0]
        #print 'WRF LATS ', XLONG[0]
        #print 'WRF LATS ', wrf_lat
        #for i in range(len(wrf_lon)-1):
        #    print wrf_lon[i+1]-wrf_lon[i]
        #for i in range(len(XLONG[0,:])-1):
        #    print XLONG[0,i]-XLONG[0,i+1]
        lat_grid_a=np.arange(lrc_lat,ulc_lat+2*0.2,0.2)[:-1]
        lon_grid_a=np.arange(ulc_lon,lrc_lon+2*0.2,0.2)[:-1]
        lon, lat = np.meshgrid(lon_grid_a,lat_grid_a,sparse=False)
        print np.shape(u_norm[hr,:,:,:])
        aa14=m14.quiver(x2, y2, np.mean(u_norm[hr,0:8,:,:],axis=0), np.mean(v_norm[hr,0:8,:,:],axis=0), pivot='middle',ax=ax14)
#        print np.mean(u_norm[hr,5,:,:]), np.mean(v_norm[hr,5,:,:])
        counter=1

######################### PLOT OPERATIONAL
        co_gfed4s_days_mean2=np.zeros((720,1440),)
        co_gfed4s_days_mean4=np.zeros((720,1440),)
        co_gfed4s_days_mean2 = 0.

        m5, trop_lat, trop_lon, trop = main(co_gfed4s_days_mean2,co_gfed4s_days_mean4,conconly=conconly,PlotData=plotdata,plot=region, \
             contourfire=contourfire,gfedonly=gfedonly,mask_ocean=mask_ocean,proj=project,RegridData=regriddata, \
             smoothing=smoothing,dpi=dpi, Step=step, Stepx=Stepx,Species= species,GfedSymbol=gfedsymbol,Gfed=gfed,outputdir=output, \
             origpath=origpath,timeperiod = tp,st_time= st_time  ,ed_time= ed_time  ,reReadData= rereaddata, \
             PolygonPlot=polygon,wrf=wrf, kernel=kernel)
######################### PLOT OPERATIONAL

        C3 = np.subtract(C[:,:], trop.transpose())
        C3 = ma.array(C3,mask=(np.isnan(C3)))
        print 'MAX = ', np.max(C3), ' MIN = ', np.min(C3)
        aa6=m6.pcolormesh(x,y,C3[:,:],vmin=-200,vmax=200,cmap=result2,ax=ax6)

#        aa5=m5.pcolormesh(x,y,C3[:,:],vmin=0,vmax=200,cmap=result,ax=ax5)
        for m in [m1,m2,m3,m14,m5,m6]:
#        for m in [m14]:
            m.drawcoastlines()
#            m.drawcountries()
            parallels = np.arange(-90.,90,10.)
#            m.drawparallels(parallels,labels=[False,False,False,False])
            meridians = np.arange(-180.,180.,20.)
            m.drawmeridians(meridians,labels=[False,False,False,False])
            for state in ['Punjab_india', 'Haryana_india', 'UP', 'Bihar']:
                m.readshapefile('state_outlines/'+state+'/'+state+'-polygon',state)
        for ax,aa in zip([ax1,ax2,ax3,ax6],[aa1,aa2,aa3,aa6]):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(aa,ax=ax,cax=cax)
#            if ax==ax1 or ax==ax2 or ax==ax3 or ax==ax6:
                ##if ax==ax1:
                ##    ax.set_title('TROPOMI Carbon monoxide [mol/m2]')        
            if ax==ax1:
                    #ax.set_title('GFED + EDGAR \n CO [mol/m2]')
                ax.set_title('GFED + EDGAR CO [ppb]')
            if ax==ax2:
                   #ax.set_title('GFED CO [mol/m2]')
                ax.set_title('GFED CO [ppb]')
            if ax==ax3:
                    #ax.set_title('EDGAR CO [mol/m2]')
                ax.set_title('EDGAR CO [ppb]')
            if ax==ax6:
                ax.set_title('GFED + EDGAR - TROPOMI CO [ppb]')
                ##if ax==ax3:
                ##    ax.set_title('Tracer_Ens Carbon monoxide [mol/m2]')
        ax14.set_title('Wind field')
        ax5.set_title('TROPOMI CO [ppb]')
        if month==4:
            plt.suptitle("April %02d, %02d, UTC %02d"%(i+1,yr,8), fontsize=18)
        if month==5:
            plt.suptitle("May %02d, %02d, UTC %02d"%(i+1,yr,8), fontsize=18)
        if month==6:
            plt.suptitle("June %02d, %02d, UTC %02d"%(i+1,yr,8), fontsize=18)
        if month==7:
            plt.suptitle("July %02d, %02d, UTC %02d"%(i+1,yr,8), fontsize=18)
        if month==8:
            plt.suptitle("August %02d, %02d, UTC %02d"%(i+1,yr,8), fontsize=18)
        if month==9:
            plt.suptitle("September %02d, %02d, UTC %02d"%(i+1,yr,8), fontsize=18)
        if month==10:
            plt.suptitle("October %02d, %02d, UTC %02d"%(i+1,yr,8), fontsize=18)
        if month==11:
            plt.suptitle("November %02d, %02d, UTC %02d"%(i+1,yr,8), fontsize=18)
        if month==12:
            plt.suptitle("December %02d, %02d, UTC %02d"%(i+1,yr,8), fontsize=18)

        #fig.tight_layout()
#SET OUTPUT DIRECTORY HERE
        dirname = '16June01/'
        if not os.path.isdir('/deos/rushilv/'+dirname):
            os.mkdir('/deos/rushilv/'+str(dirname))
        plt.savefig(str(dirname)+'/comparison_%s_%02d_%02d_%02d_windfield5.png'%(yr,month,i+1+offset,hr),dpi=150)
        plt.show()
        plt.close()
        #print wrf_lat
        #print np.subtract(wrf_lat, trop_lat)
        #print trop.transpose()[0], C3[0]
        i = i+1

