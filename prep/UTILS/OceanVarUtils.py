#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 12:05:16 2023

@author: mario
"""
import sys
sys.path.insert(1, '/Users/mario/CMCC/PYTHON_SCRIPT/UTILS')
from MyUtils import find_nearest
from netCDF4 import Dataset

def read_obstat_screen(filename):
    nc=Dataset(filename)
    lon=nc.variables['LON'][:].squeeze()
    lat=nc.variables['LAT'][:].squeeze()
    time=nc.variables['TIME'][:].squeeze()
    depth=nc.variables['DEPTH'][:].squeeze()
    value=nc.variables['VALUE'][:].squeeze()
    error=nc.variables['ERROR'][:].squeeze()
    res=nc.variables['OMG'][:].squeeze()
    type_obs=nc.variables['TYPE'][:].squeeze()
    param=nc.variables['PARAM'][:].squeeze()
    inst=nc.variables['INST'][:].squeeze()
    bias=nc.variables['BIAS'][:].squeeze()
    ib=nc.variables['IB'][:].squeeze()
    jb=nc.variables['JB'][:].squeeze()
    kb=nc.variables['KB'][:].squeeze()
    bgerr=nc.variables['BGERR'][:].squeeze()
    eve=nc.variables['EVE'][:].squeeze()
    tdist=nc.variables['TDIST'][:].squeeze()
    bac=nc.variables['BAC'][:].squeeze()
    nc.close()
    return lon,lat,time,depth,value,error,res,type_obs,\
        param,inst,bias,ib,jb,kb,bgerr,eve,tdist,bac
        
def read_obstat(filename):
    #nobs = number of observations per type [INS,SLA,SST,SSS]
    #lon = longitude
    #lat = latitude
    #time = time (format?)
    
    
    nc=Dataset(filename)
    nobs=nc.variables['NOBS'][:].squeeze()
    lon=nc.variables['LON'][:].squeeze()
    lat=nc.variables['LAT'][:].squeeze()
    time=nc.variables['TIME'][:].squeeze()
    depth=nc.variables['DEPTH'][:].squeeze()
    value=nc.variables['VALUE'][:].squeeze()
    error=nc.variables['ERROR'][:].squeeze()
    resb=nc.variables['OMG'][:].squeeze()
    resa=nc.variables['OMA'][:].squeeze()
    type_obs=nc.variables['TYPE'][:].squeeze()
    param=nc.variables['PARAM'][:].squeeze()
    inst=nc.variables['INST'][:].squeeze()
    bias=nc.variables['BIAS'][:].squeeze()
    ib=nc.variables['IB'][:].squeeze()
    jb=nc.variables['JB'][:].squeeze()
    kb=nc.variables['KB'][:].squeeze()
    bgerr=nc.variables['BGERR'][:].squeeze()
    pq=nc.variables['PQ'][:].squeeze()
    inc=nc.variables['INC'][:].squeeze()
    pert=nc.variables['PERT'][:].squeeze()
    nc.close()
    return nobs,lon,lat,time,depth,value,error,resb,resa,\
        type_obs,param,inst,bias,ib,jb,kb,bgerr,pq,inc,pert
        

def read_increment(filename,klevel):
    nc=Dataset(filename)
    tem=nc.variables['INCTEMPER'][klevel,:,:].squeeze()
    sal=nc.variables['INCSALINE'][klevel,:,:].squeeze()
    try:
       ssh=nc.variables['INCSSHEIG'][:,:].squeeze()
    except KeyError:
       ssh=None
    lon=nc.variables['nav_lon'][:,:].squeeze()
    lat=nc.variables['nav_lat'][:,:].squeeze()
    nc.close()
    return tem,sal,ssh,lon,lat

def read_increment_profile(filename,lon_target,lat_target):
    nc=Dataset(filename)
    lon=nc.variables['nav_lon'][:,:].squeeze()
    lat=nc.variables['nav_lat'][:,:].squeeze()
    depth=nc.variables['deptht'][:].squeeze()
    lon=lon[0,:].squeeze()
    lat=lat[:,0].squeeze()
    idx,dummy=find_nearest(lon, lon_target)
    idy,dummy=find_nearest(lat, lat_target)
    tem=nc.variables['INCTEMPER'][:,idy,idx].squeeze()
    sal=nc.variables['INCSALINE'][:,idy,idx].squeeze()
    try:
       ssh=nc.variables['INCSSHEIG'][idy,idx].squeeze()
    except KeyError:
       ssh=None
    nc.close()
    return tem,sal,ssh,depth