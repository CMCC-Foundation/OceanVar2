#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 15:47:25 2023

@author: mario
"""

import sys
sys.path.insert(1, '/Users/mario/CMCC/PYTHON_SCRIPT/UTILS')
from OceanVarUtils import read_obstat_screen
from netCDF4 import Dataset
import numpy as np
import glob


def shuffle_obs(nind,ino,lon,lat,flg,tim,dtm,val,bac,err,res):
    nobs=len(nind)
    nind_new=np.zeros(nobs,dtype=np.int64)
    ino_new=np.zeros(nobs,dtype=np.int64)
    lon_new=np.zeros(nobs,dtype=np.float64)
    lat_new=np.zeros(nobs,dtype=np.float64)
    flg_new=np.zeros(nobs,dtype=np.int64)
    tim_new=np.zeros(nobs,dtype=np.float64)
    dtm_new=np.zeros(nobs,dtype=np.float64)
    val_new=np.zeros(nobs,dtype=np.float64)
    bac_new=np.zeros(nobs,dtype=np.float64)
    err_new=np.zeros(nobs,dtype=np.float64) 
    res_new=np.zeros(nobs,dtype=np.float64)
    for i in range(1,nobs+1):
        idx=np.where(i==nind)[0][0]
        ino_new[i-1]=ino[idx]
        lon_new[i-1]=lon[idx]
        lat_new[i-1]=lat[idx]
        flg_new[i-1]=flg[idx]
        tim_new[i-1]=tim[idx]
        dtm_new[i-1]=dtm[idx]
        val_new[i-1]=val[idx]
        bac_new[i-1]=bac[idx]
        err_new[i-1]=err[idx]
        res_new[i-1]=res[idx]
    return ino_new,lon_new,lat_new,flg_new,tim_new,\
           dtm_new,val_new,bac_new,err_new,res_new

def writebin(dirout,obs_arc,dim_arc,otype,ll_saveksat):
    if otype=='argo':
        print('Preparing ARGO ...')
        idx=np.where(obs_arc['OTYPE'][0]==831)
        ino=np.array(obs_arc['INO'][0][idx],dtype=np.int64)
        lon=np.array(obs_arc['LON'][0][idx],dtype=np.float64)
        lat=np.array(obs_arc['LAT'][0][idx],dtype=np.float64)
        flg=np.array(obs_arc['FLC'][0][idx],dtype=np.int64)
        par=np.array(obs_arc['PAR'][0][idx],dtype=np.int64)
        #swap temperature and salinity
        idxT=np.where(par==2)
        idxS=np.where(par==1)
        par[idxT]=1
        par[idxS]=2
        tim=np.array(obs_arc['TIM'][0][idx],dtype=np.float64)
        dpt=np.array(obs_arc['DPT'][0][idx],dtype=np.float64)
        val=np.array(obs_arc['VAL'][0][idx],dtype=np.float64)
        bac=np.array(obs_arc['BAC'][0][idx],dtype=np.float64)
        err=np.array(obs_arc['ERR'][0][idx],dtype=np.float64) 
        res=np.array(obs_arc['RES'][0][idx],dtype=np.float64)
        # NOT NEEDED BY 3DVAR FILL IN WITH DUMMY VALUES
        dummyi=np.zeros((len(idx[0])),dtype=np.int64)
        dummyr=np.zeros((len(idx[0])),dtype=np.float64)
        ib=dummyi
        jb=dummyi
        kb=dummyi
        pb=dummyr
        qb=dummyr
        rb=dummyr          
        # OR FILL THE WITH REAL VALUES
        # ib=np.min(obs_arc['IB'][0][:,idx],0)
        # jb=np.min(obs_arc['JB'][0][:,idx],0)
        # kb=obs_arc['KB'][0][idx]
        # pb=obs_arc['PQ'][0][:,idx] ??
        # qb=obs_arc['PQ'][0][:,idx] ??
        # rb=obs_arc['RB'][0][idx]        
        # FINALLY WRITE FILE
        fileout=dirout+'/arg_mis.dat'
        print('Writing ARGO')
        fId = open(fileout, "wb")
        np.array(8,np.int32).tofile(fId) # record size
        np.array(len(idx[0]),np.int64).tofile(fId)
        np.array(8,np.int32).tofile(fId) # record size
        np.array(8*17*len(idx[0]),np.int32).tofile(fId)# record size
        np.array(ino).tofile(fId)
        np.array(flg).tofile(fId)
        np.array(par).tofile(fId)
        np.array(lon).tofile(fId)
        np.array(lat).tofile(fId)
        np.array(dpt).tofile(fId)
        np.array(tim).tofile(fId)
        np.array(val).tofile(fId)
        np.array(bac).tofile(fId)
        np.array(err).tofile(fId)
        np.array(res).tofile(fId)
        np.array(ib).tofile(fId)
        np.array(jb).tofile(fId)
        np.array(kb).tofile(fId)
        np.array(pb).tofile(fId)
        np.array(qb).tofile(fId)
        np.array(rb).tofile(fId)
        np.array(8*17*len(idx[0]),np.int32).tofile(fId)
        fId.close()   
    elif otype=='sla': 
        print('Preparing SLA ...')
        #ino in 3dvar identifies the track
        nind=np.array(obs_arc['NIND'][0],dtype=np.int64)
        ksat=np.array(obs_arc['KSAT'][0],dtype=np.int64)
        ino=np.array(obs_arc['TRACK'][0],dtype=np.int64)
        lon=np.array(obs_arc['LON'][0],dtype=np.float64)
        lat=np.array(obs_arc['LAT'][0],dtype=np.float64)
        flg=np.array(obs_arc['FLC'][0],dtype=np.int64)
        tim=np.array(obs_arc['TIM'][0],dtype=np.float64)
        dtm=np.array(obs_arc['DPT'][0],dtype=np.float64)
        val=np.array(obs_arc['VAL'][0],dtype=np.float64)
        bac=np.array(obs_arc['BAC'][0],dtype=np.float64)
        err=np.array(obs_arc['ERR'][0],dtype=np.float64) 
        res=np.array(obs_arc['RES'][0],dtype=np.float64)
        # .... reordering obs
        ino,lon,lat,flg,tim,dtm,val,bac,err,res = shuffle_obs\
            (nind,ino,lon,lat,flg,tim,dtm,val,bac,err,res)
        # NOT NEEDED BY 3DVAR FILL IN WITH DUMMY VALUES
        dummyi=np.zeros((len(ino)),dtype=np.int64)
        dummyr=np.zeros((len(ino)),dtype=np.float64)
        ib=dummyi
        jb=dummyi
        pb=dummyr
        qb=dummyr
        # FINALLY WRITE FILE        
        fileout=dirout+'/sla_mis.dat'
        print('Writing SLA')
        fId = open(fileout, "wb")
        np.array(8,np.int32).tofile(fId) # record size
        np.array(len(ino),np.int64).tofile(fId)
        np.array(8,np.int32).tofile(fId) # record size
        if (ll_saveksat):
           np.array(8*15*len(ino),np.int32).tofile(fId)#
        else:
           np.array(8*14*len(ino),np.int32).tofile(fId)# record size
        np.array(ino).tofile(fId)
        if (ll_saveksat):
            np.array(ksat).tofile(fId)
        np.array(flg).tofile(fId)
        np.array(lon).tofile(fId)
        np.array(lat).tofile(fId)
        np.array(tim).tofile(fId)
        np.array(val).tofile(fId)
        np.array(bac).tofile(fId)
        np.array(err).tofile(fId)
        np.array(res).tofile(fId)
        np.array(ib).tofile(fId)
        np.array(jb).tofile(fId)
        np.array(pb).tofile(fId)
        np.array(qb).tofile(fId)
        np.array(dtm).tofile(fId)
        if (ll_saveksat):
           np.array(8*15*len(ino),np.int32).tofile(fId)#
        else:
           np.array(8*14*len(ino),np.int32).tofile(fId)# record size
        fId.close()        
    return

def read_grid(filename):
    nc=Dataset(filename,mode='r')
    lon=nc.variables['lon'][:]
    lat=nc.variables['lat'][:]
    nc.close()
    return lon,lat

def read_nc(filename):
    vardict={}
    dimdict={}
    nc=Dataset(filename,mode='r')
    var_list=list(nc.variables.keys())
    dim_list=list(nc.dimensions.keys())
    for var in var_list:
        x=nc.variables[var]
        dim_name=x.dimensions
        val=x[:]
        vardict[var]=[val,dim_name]
    for dim in dim_list:
        x=nc.dimensions[dim].size
        dimdict[dim]=x
    nc.close()
    return vardict,dimdict

def write_nc(fileout,obs_arc,dims_arc):
    nc=Dataset(fileout,mode='w',format='NETCDF4_CLASSIC')
    #Dimensions
    for dim in dims_arc.keys():
        if dim == 'PREDS':
           nc.createDimension(dim, None)
        else:
           nc.createDimension(dim, dims_arc[dim])
    #Variables
    for var in obs_arc.keys():
        dtype=obs_arc[var][0].dtype
        dimensions=obs_arc[var][0].shape
        dim_name=obs_arc[var][1]            
        #Create
        varvar=nc.createVariable(var, dtype, dim_name)
        #Fill in
        varvar[:]=obs_arc[var][0]
    nc.close()
    return

def subsample(obs_arc,dims_arc,index1,itype):
    if itype=='profile':
        index=list(np.where(obs_arc['PROF'][0]==index1)[0])
    elif  itype=='point':
        index=index1
    else:
        print('Type can be "Profile" or point')
        exit
    for k in obs_arc.keys():
        if obs_arc[k][1][0]=='OBS':
            obs_arc[k][0]=obs_arc[k][0][index]
        else:
            obs_arc[k][0]=obs_arc[k][0][:,index]
    dims_arc['OBS']=len(index)
    obs_arc['NIND'][0]=np.linspace(1,len(index),len(index))
    return obs_arc,dims_arc


#------------------------------------------------------------------------------
#-------BEGIN -----------------------------------------------------------------
#------------------------------------------------------------------------------
dirin='/Users/mario/CMCC/OV_rc_v0.2/tmp_shuffle_obs'

filegrid='/Users/mario/CMCC/test_mario_MedSea_red/tmp_shuffle_obs/GRID.nc'
filegrid_hr='/Users/mario/CMCC/OV_rc_v0.2/tmp_shuffle_obs/GRID.nc'
obserrorfile='/Users/mario/CMCC/DATA/ZEUS/OBSSTAT_SCREEN.NC.0000'


ll_subsample_ssh=False
ll_subsample_ins=False
ll_subsample_profile=False
ll_saveXoceanvar=False
ll_saveksat=False

if ll_saveXoceanvar:
    dirout='/Users/mario/CMCC/test_mario_MedSea_red/tmp_shuffle_obs'
else:
    dirout='/Users/mario/CMCC/3DVar/test_case_mario_MedSea_red' 
    
if (ll_subsample_ssh):
    #shallow
    #index_ssh=[250]
    #deep
    index_ssh=[275]
if (ll_subsample_ins):
    # shallow S
    # index_ins=[100]
    # shallow T
    # index_ins=[99]    
    # deep S 64.5 m
    # index_ins=[1400]
    # deep T 64.5 m
    index_ins=[1399]    
    ll_subsample_profile=False
if (ll_subsample_profile): 
    index_prof=[10]

filetype_list=['INSMIS_','SLAMIS_']
# filetype_list=['SLAMIS_']

lon,lat = read_grid(filegrid)
lon=lon[0,:]
lat=lat[:,0]
lon_hr,lat_hr = read_grid(filegrid_hr)

for filetype in filetype_list:
    print('Preparing: '+filetype)
    dirtype=dirin+'/'+filetype
    file_list=sorted(glob.glob(dirtype+"*.NC"))
    first_time=True
    for file_name in file_list:
        print(file_name)
        obs,dims=read_nc(file_name)
        if first_time:
            dims_arc=dims
            key_list=obs.keys()
            obs_arc=obs
            first_time=False
        else:
            for key in key_list:
                if obs[key][1][0]=='OBS':
                    obs_arc[key][0]=np.concatenate((obs_arc[key][0].data, obs[key][0].data),axis=0)
                else:
                    obs_arc[key][0]=np.concatenate((obs_arc[key][0].data, obs[key][0].data),axis=1)
            dims_arc['OBS']=dims_arc['OBS']+dims['OBS']
            
    #REARRANGE OBS
    print('JB = '+str(obs['JB'][0].data[:,0]))
    print('IB = '+str(obs['IB'][0].data[:,0]))
    print("Shouldn't be something like: ")
    print("      SE  SW  NE  NW")
    print('JB = '+str(obs_arc['JB'][0][:,0]))
    print('IB = '+str(obs_arc['IB'][0][:,0]))  
    print("Lon, Lat Values do not change but Gird point does!!!")
    for nobs in range(len(obs_arc['LON'][0])):
        lonobs=obs_arc['LON'][0][nobs]
        latobs=obs_arc['LAT'][0][nobs]
        
        lonmod1=lon[(lon-lonobs)<0][-1]
        lonmod2=lon[(lon-lonobs)>0][0]
        latmod1=lat[(lat-latobs)<0][-1]
        latmod2=lat[(lat-latobs)>0][0]
        # In Fortran indices start from 1.
        # In Python start from 0.
        # correct to add +1 
        # However to have same number of obs in situ it is necessary to 
        # eliminate +1.
        # To have the same number of observations for sla add +1
        i1=np.where(lon==lonmod1)[0][0]+1 #E            
        i2=np.where(lon==lonmod2)[0][0]+1 #W
        j1=np.where(lat==latmod1)[0][0]+1 #S
        j2=np.where(lat==latmod2)[0][0]+1 #N
        # i1=np.where(lon==lonmod1)[0][0] #E            
        # i2=np.where(lon==lonmod2)[0][0] #W
        # j1=np.where(lat==latmod1)[0][0] #S
        # j2=np.where(lat==latmod2)[0][0] #N
        
        obs_arc['IB'][0][:,nobs]=[i1,i2,i1,i2]
        obs_arc['JB'][0][:,nobs]=[j1,j1,j2,j2]
        
    if (ll_subsample_ssh and filetype=='SLAMIS_'):    
        obs_arc,dims_arc=subsample(obs_arc,dims_arc,index_ssh,'point')
    if (ll_subsample_ins and filetype=='INSMIS_'):    
        obs_arc,dims_arc=subsample(obs_arc,dims_arc,index_ins,'point') 
    if (ll_subsample_profile and filetype=='INSMIS_'):
        obs_arc,dims_arc=subsample(obs_arc,dims_arc,index_prof,'profile')
    
    
    if ll_saveXoceanvar:                        #WRITE NETCDF X OCEANVAR
        fileout=dirout+'/'+filetype+'0000.NC'
        write_nc(fileout,obs_arc,dims_arc)    
    else:                                       #WRITE 3DVAR
        #  Load err from OBSTAT_SCREN file
        lonERR,latERR,timeERR,depthERR,valueERR,errorERR,resERR,type_obsERR,\
            paramERR,instERR,biasERR,ibERR,jbERR,kbERR,bgerrERR,eveERR,tdistERR,bacERR=\
            read_obstat_screen(obserrorfile)
        #  Initialize obs error array 
        obs_arc['ERR']=[np.zeros(obs_arc['LON'][0].shape),\
                        obs_arc['LON'][1]]
        if filetype=='INSMIS_':
            #... fill inn obs error array for in situ..
            for i in range(len(obs_arc['LON'][0])):
                idx=np.where( (obs_arc['LON'][0][i]==lonERR) & \
                              (obs_arc['LAT'][0][i]==latERR) & \
                              (obs_arc['TIM'][0][i]==timeERR)& \
                              (obs_arc['DPT'][0][i]==depthERR)&\
                              (obs_arc['PAR'][0][i]==paramERR) ) 
                if not idx[0].size == 0:
                   obs_arc['ERR'][0][i]=errorERR[idx]
            obs_arc['ERR'][1]=obs_arc['LON'][1]  
            # finally write file...
            writebin(dirout,obs_arc,dims_arc,'argo',ll_saveksat)
            
        elif filetype=='SLAMIS_':
            #... fill inn obs error array for sla..
            for i in range(len(obs_arc['LON'][0])):
                idx=np.where( (obs_arc['LON'][0][i]==lonERR) & \
                              (obs_arc['LAT'][0][i]==latERR) & \
                              (obs_arc['TIM'][0][i]==timeERR) )
                if not idx[0].size == 0:
                   obs_arc['ERR'][0][i]=errorERR[idx]
            #adjust for not having negative error
            idx=tuple([obs_arc['ERR'][0]>0])
            ave=np.mean(obs_arc['ERR'][0][idx])
            idx=tuple([obs_arc['ERR'][0]<=0])
            obs_arc['ERR'][0][idx]=ave
            writebin(dirout,obs_arc,dims_arc,'sla',ll_saveksat)
        else:
            print('Filetype not supported')
        