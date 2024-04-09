
import cftime
import xarray as xr
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import time
import numpy.ma as ma

def basic_composite(tos,regbox,ninobox,times,strt,fin,strtstr,endstr,yrs):
    '''
    shay magahey UCSB
    03/06/2024
    identifies el nino years and averages over those times, calculates composites based on htat
    tos: SST variable
    regbox: total region to consider
    ninobox: region for compositing, usually n34
    times: string with times of interest
    strt,fin: start and end times as integers (ie 2070,2100)
    strtstr,endstr: start and end times as strings with month only
    yrs: array with year strings
    returns ninocomp,ninacomp: el nino and la nina composites over regbox for each member (mem x lat x lon)
    '''
    
    ens_sst=tos.sel(lat=slice(regbox[0],regbox[1]),lon=slice(regbox[2],regbox[3]),time=slice(times[0],times[1]))
    ens_n34=tos.sel(lat=slice(ninobox[0],ninobox[1]),lon=slice(ninobox[2],ninobox[3]),time=slice(times[0],times[1]))
        
    # Calculate NINO3.4 regional average
    ens_n34=ens_n34.mean(dim="lat").mean(dim="lon").squeeze()
                
    # Detrend NINO3.4 by removing the ensemble mean
    ens_n34=ens_n34-ens_n34.mean('member')
        
    # Remove climatology
    clim=ens_n34.groupby('time.month').mean()
    ens_n34=ens_n34.groupby('time.month')-clim  

    clim=ens_sst.groupby('time.month').mean()
    ens_sst=ens_sst.groupby('time.month')-clim  

    # Calculate DJF averages
    tmpdjf=ens_n34 #.resample(time='QS-DEC').mean(dim="time")
    tmpdjf=tmpdjf.where(tmpdjf.time.dt.month == 1,drop="True") # n34 sst
    ens_n34_djf=tmpdjf
    
    tmpdjf1=ens_sst #.resample(time='QS-DEC').mean(dim="time")
    tmpdjf1=tmpdjf1.where(tmpdjf1.time.dt.month == 1,drop="True") # all sst
    ens_sst_djf=tmpdjf1   
    
    # Identify El Nino, La Nina years
    tmpn34std=ens_n34_djf.std('time').mean('member').squeeze() 
    msk_en=ens_n34_djf > tmpn34std
    en_evt=ens_sst_djf.where(msk_en, drop=True).squeeze()

    msk_ln=ens_n34_djf < -1*tmpn34std
    ln_evt=ens_sst_djf.where(msk_ln, drop=True).squeeze()

    yrrange=np.arange(strt,fin,1)



    enyrs_per=[]

    for i in range(0,16):
        A=msk_en[i,:]
        #print(A)
        fut=ma.masked_array(yrrange,mask=~A[0:30],fill_value=-99999)
        #m=A[0:30]
        #fut=xr.where(futyrs!=np.nan,futyrs,np.nan)
        #print(fut)
        fut=fut.filled(-99999999)
        enyrs_per.append(fut)
     # this is giving the nino years!!  
    ninoyrs_mem=np.stack(enyrs_per,axis=0)
    ninoyrs_mem=np.where(ninoyrs_mem!=-99999999,ninoyrs_mem,np.nan)

    lnyrs_per=[]
    for i in range(0,16):
        B=msk_ln[i,:]
        #print(A)
        futln=ma.masked_array(yrrange,mask=~B[0:30],fill_value=-99999)
        #m=A[0:30]
        #fut=xr.where(futyrs!=np.nan,futyrs,np.nan)
        #print(fut)
        futln.astype(int)
        futln=futln.filled(-99999999)
        lnyrs_per.append(futln)
     # this is giving the nino years!!  
    ninayrs_mem=np.stack(lnyrs_per,axis=0)
    ninayrs_mem=np.where(ninayrs_mem!=-99999999,ninayrs_mem,np.nan)

    # convert things to datetime because cftime is moderately useless
    # next want to just get SSTs for years that are el nino/la nina

    
    time = pd.date_range(start=strtstr, periods=30, freq='AS')
    time_coord_np = np.array(time)
    yrs=yrs.rename({'dim_0':'time'})
    # Assign the time coordinate to the data array
    yrs['time'] = ('time', time_coord_np)
    yrs['time'] = np.array(yrs['time'].values, dtype='datetime64')
    dts=pd.date_range(start=strtstr,end=endstr,freq='MS')
    #sst['time']=dts

    mems=np.array([0,1,2,3,4,5,6,7,8,10,11,12,13,14,15])
    enlist=[]
    for m in mems:
        print('member=',m)
        thismem=xr.DataArray(ninoyrs_mem[m,:]).rename({'dim_0':'time'}) # set of yrs for this member
        thisyrsen=yrs.where(~thismem.isnull(),drop=True) # el nino yrs
        sstEN=ens_sst[m,:,:,:] # 361x10x160
        #fill_value = -999
        #arr_filled=thisyrs.where(~thisyrs.isnull(), fill_value)
        sstEN['time']=dts
        sstEN=sstEN.convert_calendar("proleptic_gregorian",use_cftime=False)
        #print(thisyrsen)
        #print(sstEN.time)
        sstENyrs=sstEN.sel(time=thisyrsen).mean('time')
        enlist.append(sstENyrs)
    ninocomp=xr.concat(enlist,dim='member')
    # ERROR: not all index values found in time: print time: probably what you were doung yesterday w format
    
    lnlist=[]
    for m in mems:
        print('member=',m)
        thismem=xr.DataArray(ninayrs_mem[m,:]).rename({'dim_0':'time'}) # set of yrs for this member
        thisyrsln=yrs.where(~thismem.isnull(),drop=True) # el nino yrs
        sstLN=ens_sst[m,:,:,:] # 361x10x160
        #fill_value = -999
        #arr_filled=thisyrs.where(~thisyrs.isnull(), fill_value)
        sstLN['time']=dts
        sstLN=sstLN.convert_calendar("proleptic_gregorian",use_cftime=False)
        sstLNyrs=sstLN.sel(time=thisyrsln).mean('time')
        lnlist.append(sstLNyrs)
    ninacomp=xr.concat(lnlist,dim='member')
    return ninocomp,ninacomp


def composite_anyvar(tos,var,regbox,ninobox,times,strt,fin,strtstr,endstr,yrs):
    '''
    shay magahey UCSB
    03/06/2024
    identifies el nino years and averages over those times, calculates composites based on htat
    tos: SST variable
    var: variable to be composited based on nino yrs. should alrady be detrended
    regbox: total region to consider
    ninobox: region for compositing, usually n34
    times: string with times of interest
    strt,fin: start and end times as integers (ie 2070,2100)
    strtstr,endstr: start and end times as strings with month only
    yrs: array with year strings
    returns ninocomp,ninacomp: el nino and la nina composites over regbox for each member (mem x lat x lon)
    '''
    
    ens_sst=tos.sel(lat=slice(regbox[0],regbox[1]),lon=slice(regbox[2],regbox[3]),time=slice(times[0],times[1]))
    ens_n34=tos.sel(lat=slice(ninobox[0],ninobox[1]),lon=slice(ninobox[2],ninobox[3]),time=slice(times[0],times[1]))

    anyvar=anyvar.sel(lat=slice(regbox[0],regbox[1]),lon=slice(regbox[2],regbox[3]),time=slice(times[0],times[1]))
    
    # Calculate NINO3.4 regional average
    ens_n34=ens_n34.mean(dim="lat").mean(dim="lon").squeeze()
                
    # Detrend NINO3.4 by removing the ensemble mean
    ens_n34=ens_n34-ens_n34.mean('member')
        
    # Remove climatology
    clim=ens_n34.groupby('time.month').mean()
    ens_n34=ens_n34.groupby('time.month')-clim  

    clim=ens_sst.groupby('time.month').mean()
    ens_sst=ens_sst.groupby('time.month')-clim  

    # Calculate DJF averages
    tmpdjf=ens_n34 #.resample(time='QS-DEC').mean(dim="time")
    tmpdjf=tmpdjf.where(tmpdjf.time.dt.month == 1,drop="True") # n34 sst
    ens_n34_djf=tmpdjf
    
    tmpdjf1=ens_sst #.resample(time='QS-DEC').mean(dim="time")
    tmpdjf1=tmpdjf1.where(tmpdjf1.time.dt.month == 1,drop="True") # all sst
    ens_sst_djf=tmpdjf1   
    
    # Identify El Nino, La Nina years
    tmpn34std=ens_n34_djf.std('time').mean('member').squeeze() 
    msk_en=ens_n34_djf > tmpn34std
    en_evt=ens_sst_djf.where(msk_en, drop=True).squeeze()

    msk_ln=ens_n34_djf < -1*tmpn34std
    ln_evt=ens_sst_djf.where(msk_ln, drop=True).squeeze()

    yrrange=np.arange(strt,fin,1)

    enyrs_per=[]

    for i in range(0,16):
        A=msk_en[i,:]
        #print(A)
        fut=ma.masked_array(yrrange,mask=~A[0:30],fill_value=-99999)
        #m=A[0:30]
        #fut=xr.where(futyrs!=np.nan,futyrs,np.nan)
        #print(fut)
        fut=fut.filled(-99999999)
        enyrs_per.append(fut)
     # this is giving the nino years!!  
    ninoyrs_mem=np.stack(enyrs_per,axis=0)
    ninoyrs_mem=np.where(ninoyrs_mem!=-99999999,ninoyrs_mem,np.nan)

    lnyrs_per=[]
    for i in range(0,16):
        B=msk_ln[i,:]
        #print(A)
        futln=ma.masked_array(yrrange,mask=~B[0:30],fill_value=-99999)
        #m=A[0:30]
        #fut=xr.where(futyrs!=np.nan,futyrs,np.nan)
        #print(fut)
        futln.astype(int)
        futln=futln.filled(-99999999)
        lnyrs_per.append(futln)
     # this is giving the nino years!!  
    ninayrs_mem=np.stack(lnyrs_per,axis=0)
    ninayrs_mem=np.where(ninayrs_mem!=-99999999,ninayrs_mem,np.nan)

    # convert things to datetime because cftime is moderately useless
    # next want to just get SSTs for years that are el nino/la nina

    
    time = pd.date_range(start=strtstr, periods=30, freq='AS')
    time_coord_np = np.array(time)
    yrs=yrs.rename({'dim_0':'time'})
    # Assign the time coordinate to the data array
    yrs['time'] = ('time', time_coord_np)
    yrs['time'] = np.array(yrs['time'].values, dtype='datetime64')
    dts=pd.date_range(start=strtstr,end=endstr,freq='MS')
    #sst['time']=dts

    mems=np.array([0,1,2,3,4,5,6,7,8,10,11,12,13,14,15])
    enlist=[]
    for m in mems:
        print('member=',m)
        thismem=xr.DataArray(ninoyrs_mem[m,:]).rename({'dim_0':'time'}) # set of yrs for this member
        thisyrsen=yrs.where(~thismem.isnull(),drop=True) # el nino yrs
        sstEN=anyvar[m,:,:,:] # 361x10x160
        #fill_value = -999
        #arr_filled=thisyrs.where(~thisyrs.isnull(), fill_value)
        sstEN['time']=dts
        sstEN=sstEN.convert_calendar("proleptic_gregorian",use_cftime=False)
        #print(thisyrsen)
        #print(sstEN.time)
        sstENyrs=sstEN.sel(time=thisyrsen).mean('time')
        enlist.append(sstENyrs)
    ninoany=xr.concat(enlist,dim='member')
    # ERROR: not all index values found in time: print time: probably what you were doung yesterday w format
    
    lnlist=[]
    for m in mems:
        print('member=',m)
        thismem=xr.DataArray(ninayrs_mem[m,:]).rename({'dim_0':'time'}) # set of yrs for this member
        thisyrsln=yrs.where(~thismem.isnull(),drop=True) # el nino yrs
        sstLN=anyvar[m,:,:,:] # 361x10x160
        
        sstLN['time']=dts
        sstLN=sstLN.convert_calendar("proleptic_gregorian",use_cftime=False)
        sstLNyrs=sstLN.sel(time=thisyrsln).mean('time')
        lnlist.append(sstLNyrs)
    ninaany=xr.concat(lnlist,dim='member')
    return ninoany,ninaany


def evolution_composites(tos,regbox,ninobox,times,strt,fin,strtstr,endstr,yrs):
    '''
    shay magahey ucsb
    03/06/2024
    function to get the hovmoller style data: event evolution composites
    tos: sst data
    regbox: overall AOI
    ninobox: nino region for compositing (n34 probs)
    times: times of interest
    strt,fin: start and finish years as integers
    strtstr,endstr: strings of start and end dates
    yrs: array of years as strings
    returns composites with monthly evolution
    
    '''

    ens_sst=tos.sel(lat=slice(regbox[0],regbox[1]),lon=slice(regbox[2],regbox[3]),time=slice(times[0],times[1]))
    ens_n34=tos.sel(lat=slice(ninobox[0],ninobox[1]),lon=slice(ninobox[2],ninobox[3]),time=slice(times[0],times[1]))
        
    # Calculate NINO3.4 regional average
    ens_n34=ens_n34.mean(dim="lat").mean(dim="lon").squeeze()
                
    # Detrend NINO3.4 by removing the ensemble mean
    ens_n34=ens_n34-ens_n34.mean('member')
        
    # Remove climatology
    clim=ens_n34.groupby('time.month').mean()
    ens_n34=ens_n34.groupby('time.month')-clim  

    clim=ens_sst.groupby('time.month').mean()
    ens_sst=ens_sst.groupby('time.month')-clim  

    # Calculate DJF averages
    tmpdjf=ens_n34 #.resample(time='QS-DEC').mean(dim="time")
    tmpdjf=tmpdjf.where(tmpdjf.time.dt.month == 1,drop="True") # n34 sst
    ens_n34_djf=tmpdjf
    
    tmpdjf1=ens_sst #.resample(time='QS-DEC').mean(dim="time")
    tmpdjf1=tmpdjf1.where(tmpdjf1.time.dt.month == 1,drop="True") # all sst
    ens_sst_djf=tmpdjf1   
    
    # Identify El Nino, La Nina years
    tmpn34std=ens_n34_djf.std('time').mean('member').squeeze() 
    msk_en=ens_n34_djf > tmpn34std
    en_evt=ens_sst_djf.where(msk_en, drop=True).squeeze()

    msk_ln=ens_n34_djf < -1*tmpn34std
    ln_evt=ens_sst_djf.where(msk_ln, drop=True).squeeze()

    yrrange=np.arange(strt,fin,1)

    enyrs_per=[]

    for i in range(0,16):
        A=msk_en[i,:]
        #print(A)
        fut=ma.masked_array(yrrange,mask=~A[0:30],fill_value=-99999)
        #m=A[0:30]
        #fut=xr.where(futyrs!=np.nan,futyrs,np.nan)
        #print(fut)
        fut=fut.filled(-99999999)
        enyrs_per.append(fut)
     # this is giving the nino years!!  
    ninoyrs_mem=np.stack(enyrs_per,axis=0)
    ninoyrs_mem=np.where(ninoyrs_mem!=-99999999,ninoyrs_mem,np.nan)

    lnyrs_per=[]
    for i in range(0,16):
        B=msk_ln[i,:]
        #print(A)
        futln=ma.masked_array(yrrange,mask=~B[0:30],fill_value=-99999)
        #m=A[0:30]
        #fut=xr.where(futyrs!=np.nan,futyrs,np.nan)
        #print(fut)
        futln.astype(int)
        futln=futln.filled(-99999999)
        lnyrs_per.append(futln)
     # this is giving the nino years!!  
    ninayrs_mem=np.stack(lnyrs_per,axis=0)
    ninayrs_mem=np.where(ninayrs_mem!=-99999999,ninayrs_mem,np.nan)

    # convert things to datetime because cftime is moderately useless
    #yrs=xr.DataArray(['2070','2071','2072','2073','2074','2075','2076','2077','2078','2079','2080','2081','2082','2083','2084','2085','2086','2087','2088','2089','2090','2091','2092','2093','2094','2095','2096','2097','2098','2099'])
    time = pd.date_range(start=strtstr, periods=30, freq='AS')
    time_coord_np = np.array(time)
    yrs=yrs.rename({'dim_0':'time'})
    # Assign the time coordinate to the data array
    yrs['time'] = ('time', time_coord_np)
    yrs['time'] = np.array(yrs['time'].values, dtype='datetime64')
    dts=pd.date_range(start=strtstr,end=endstr,freq='MS')
    
    memEN=[]
    mems=np.array([0,2,3,4,5,6,7,8,10,11,12,13,14,15])
    for m in mems:
        print('member=',m)
        thismem=xr.DataArray(ninoyrs_mem[m,:]).rename({'dim_0':'time'}) # set of yrs for this member
        thisyrsen=yrs.where(~thismem.isnull(),drop=True)
        sstEN=ens_sst[m,:,:,:]
        #fill_value = -999
        #arr_filled=thisyrs.where(~thisyrs.isnull(), fill_value)
        sstEN['time']=dts
        sstEN=sstEN.convert_calendar("proleptic_gregorian",use_cftime=False)
        monthsEN=[]
        #rangeyrs=np.arange(1,len(thisyrs),1)
    
        for s in range(1,len(thisyrsen)-1):   
            print(s)
            start=np.array(thisyrsen[s].values,dtype='datetime64[M]')-np.timedelta64(12,'M')
            end=np.array(thisyrsen[s].values,dtype='datetime64[M]')+np.timedelta64(12,'M')
            print("start=",start)
            print("end=",end)
            sst1=sstEN.sel(time=slice(start,end),lat=slice(-5,5),lon=slice(120,280))
            #print("sst=",sst1) 
            monthsEN.append(sst1)
        #events=xr.concat(months,dim='mo') 
        eventsEN=np.stack(np.array(monthsEN),axis=0)
        meanmonthsEN=np.nanmean(eventsEN,axis=0) # mean of empty slice?
        #print("meanmonths shape=",meanmonths.shape)
        memEN.append(meanmonthsEN)

    monthsmemsEN=np.stack(memEN,axis=0)

    memLN=[]
    for m in mems:
        print('member=',m)
        thismem=xr.DataArray(ninayrs_mem[m,:]).rename({'dim_0':'time'}) # set of yrs for this member
        thisyrsln=yrs.where(~thismem.isnull(),drop=True)
        sstLN=ens_sst[m,:,:,:]
        #fill_value = -999
        #arr_filled=thisyrs.where(~thisyrs.isnull(), fill_value)
        sstLN['time']=dts
        sstLN=sstLN.convert_calendar("proleptic_gregorian",use_cftime=False)
        monthsLN=[]
        #rangeyrs=np.arange(1,len(thisyrs),1)
    
        for s in range(1,len(thisyrsln)-1):   
            print(s)
            start=np.array(thisyrsln[s].values,dtype='datetime64[M]')-np.timedelta64(12,'M')
            end=np.array(thisyrsln[s].values,dtype='datetime64[M]')+np.timedelta64(12,'M')
            print("start=",start)
            print("end=",end)
            sst2=sstLN.sel(time=slice(start,end),lat=slice(-5,5),lon=slice(120,280))
            #print("sst=",sst2)
            monthsLN.append(sst2)
            #print("monthsLN=",monthsLN)
        #events=xr.concat(months,dim='mo')
        eventsLN=np.stack(np.array(monthsLN),axis=0)
        meanmonthsLN=np.nanmean(eventsLN,axis=0)
        #print("meanmonths shape=",meanmonths.shape)
        memLN.append(meanmonthsLN)

    monthsmemsLN=np.stack(memLN,axis=0)

    LN_2d=np.nanmean(np.nanmean(monthsmemsLN,axis=2),axis=0)
    EN_2d=np.nanmean(np.nanmean(monthsmemsEN,axis=2),axis=0)

    return monthsmemsEN,monthsmemsLN,EN_2d,LN_2d


