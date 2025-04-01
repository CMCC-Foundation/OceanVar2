#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 17:58:50 2024

@author: mario
"""
import sys
#sys.path = ["/users_home/cmcc/ma20223/py"] + sys.path
from netCDF4 import Dataset
import numpy as np
import glob
import ttide as tt
from ttide import t_getconsts as tgc
from ttide import time as tm

#Define functions
def find_nearest(array, value):
    #it needs numpy module
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def read_tides(filename):
    nc=Dataset(filename)
    out={}
    out['M2_Amp']=nc.variables['M2_Amp'][:]
    out['K1_Amp']=nc.variables['K1_Amp'][:]
    out['O1_Amp']=nc.variables['O1_Amp'][:]
    out['S2_Amp']=nc.variables['S2_Amp'][:]
    out['P1_Amp']=nc.variables['P1_Amp'][:]
    out['N2_Amp']=nc.variables['N2_Amp'][:]
    out['Q1_Amp']=nc.variables['Q1_Amp'][:]
    out['K2_Amp']=nc.variables['K2_Amp'][:]

    out['M2_Pha']=nc.variables['M2_Pha'][:]
    out['K1_Pha']=nc.variables['K1_Pha'][:]
    out['O1_Pha']=nc.variables['O1_Pha'][:]
    out['S2_Pha']=nc.variables['S2_Pha'][:]
    out['P1_Pha']=nc.variables['P1_Pha'][:]
    out['N2_Pha']=nc.variables['N2_Pha'][:]
    out['Q1_Pha']=nc.variables['Q1_Pha'][:]
    out['K2_Pha']=nc.variables['K2_Pha'][:]
    msk=nc.variables['sossheig'][0,:,:].squeeze()
    msk[~msk.mask]=1
    out['msk']=msk
    out['lon']=nc.variables['lon'][:]
    out['lat']=nc.variables['lat'][:]

    nc.close()
    return out

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
    if otype=='INSMIS_':
        instype_list=['argo','xbt','glider']
        for instype in instype_list:
            print('Preparing '+instype+'...')
            if instype=='argo':
                idx=np.where(obs_arc['OTYPE'][0]==831)
                fileout=dirout+'/arg_mis.dat'
            elif  instype=='xbt':
                idx=np.where(obs_arc['OTYPE'][0]==401)
                fileout=dirout+'/xbt_mis.dat'
            elif instype=='glider':
                idx=np.where(obs_arc['OTYPE'][0]==820)
                fileout=dirout+'/gld_mis.dat'
            ino=np.array(obs_arc['PROF'][0][idx],dtype=np.int64)
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
            # FINALLY WRITE FILE
            print('Writing ')
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
    elif otype=='SLAMIS_':
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

#---------------------------------------------------------
# Program starts here
dirin = sys.argv[1]
dirou =  sys.argv[2]
sel_track=723

#dirin = '/Users/marioadani/SCRATCH/2020/'
#dirou = '/Users/marioadani/SCRATCH/2020/'

filetides=dirou+'/amppha2D_0_sossheig_20210701_20211231_mod_EAS7.nc'
# Read tides constituents
tides   = read_tides(filetides)
# Model tides
const, sat, shallow = tgc.t_getconsts(np.empty(0))
med_tide=tt.TTideCon()
med_tide['nameu']=np.array(['M2  ','K1  ','O1  ','S2  ',\
                            'P1  ','N2  ','Q1  ','K2  '],dtype='S4')
nconst=len(med_tide['nameu'])
med_tide['fu']=np.zeros(nconst)
for i in range(0,nconst):
    idx=np.where(med_tide['nameu'][i]==const['name'])[0][0]
    med_tide['fu'][i]=const['freq'][idx]
med_tide['ltype']='nodal'
med_tide['synth']=0
ref_time=tm.datetime(1950,1,1)

filetype_list=['SLAMIS_','INSMIS_']
#this key is for compatibility with the first version of Srdjan 3Dvar.
#it should be always True except if you are using the original Srdjan code
ll_saveksat=True
ll_detyding=True

for filetype in filetype_list:
    print('Preparing: '+filetype)
    dirtype=dirin+'/'+filetype
    file_list=sorted(glob.glob(dirtype+"*.NC"))
    first_time=True
    # aggregate all the observation in a domain
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
    #  Initialize obs error array
    obs_arc['ERR']=[np.zeros(obs_arc['LON'][0].shape),\
                    obs_arc['LON'][1]]
    # Exclude Nan... don't know why are present
    for keys in obs_arc.keys():
        try:
            idx = np.argwhere(np.isnan(obs_arc[keys][0]))
            print('keys : '+keys+' Nan found: '+str(len(idx)))
            obs_arc['FLC'][0][idx] = 0
            idx = np.argwhere(abs(obs_arc[keys][0])>1e20)
            print('keys : '+keys+' Value gt: '+str(len(idx)))
            obs_arc['FLC'][0][idx] = 0
        except TypeError:
            pass
    if filetype=='SLAMIS_' and ll_detyding:
        print('Detiding')
        count = 0
        for lon,lat,time in zip(obs_arc['LON'][0],obs_arc['LAT'][0],obs_arc['TIM'][0]):
            x,val = find_nearest(tides['lon'], lon)
            y,val = find_nearest(tides['lat'], lat)
            if tides['msk'][y,x]==1:
                med_tide['tidecon']=np.array([ [tides['M2_Amp'][y,x],0,tides['M2_Pha'][y,x],0],\
                                               [tides['K1_Amp'][y,x],0,tides['K1_Pha'][y,x],0],\
                                               [tides['O1_Amp'][y,x],0,tides['O1_Pha'][y,x],0],\
                                               [tides['S2_Amp'][y,x],0,tides['S2_Pha'][y,x],0],\
                                               [tides['P1_Amp'][y,x],0,tides['P1_Pha'][y,x],0],\
                                               [tides['N2_Amp'][y,x],0,tides['N2_Pha'][y,x],0],\
                                               [tides['Q1_Amp'][y,x],0,tides['Q1_Pha'][y,x],0],\
                                               [tides['K2_Amp'][y,x],0,tides['K2_Pha'][y,x],0] ],\
                                           dtype='>f8')
                med_tide['lat']=tides['lat'][y]
                time=tm.num2date(time+tm.date2num(ref_time))
                tide=med_tide(np.array([time]))
                # zsshcorr = obs_arc['RES'][0][count]-obs_arc['VAL'][0][count]+obs_arc['BAC'][0][count]
                obs_arc['BAC'][0][count]=obs_arc['BAC'][0][count]-tide
                # obs_arc['RES'][0][count]=obs_arc['RES'][0][count]+tide
                obs_arc['RES'][0][count]=obs_arc['VAL'][0][count]-obs_arc['BAC'][0][count]
                if tide==0:
                    obs_arc['FLC'][0][count]=0
                    print('Warning model tides = 0')
                count += 1


    # finally write file...
#    for i in range(len(obs_arc['RES'][0])):
#        print(obs_arc['RES'][0][i])
# SELECTING TRACK TO BE ASSIMILATED
        idx=obs_arc['TRACK'][0]!=sel_track
        obs_arc['FLC'][0][idx]=0
# WRITE FILE
    writebin(dirou,obs_arc,dims_arc,filetype,ll_saveksat)
