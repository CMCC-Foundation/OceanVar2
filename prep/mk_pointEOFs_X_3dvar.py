#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 12:28:47 2023

@author: mario
"""

from netCDF4 import Dataset
import numpy as np


def read_nc(filein):
    nc=Dataset(filein,mode='r')
    eigval=nc.variables['eigenvalues'][:]
    eigvec=nc.variables['eigenvectors'][:]
    regs=nc.variables['nreg'][:]
    nc.close()
    return regs,eigval,eigvec

def writenc(fileout,eva,evc):
    imt=evc.shape[3]
    jmt=evc.shape[2]
    kmt=evc.shape[1]
    neof=evc.shape[0]
    nreg=imt*jmt
    # nc=Dataset(fileout,mode='w',format='NETCDF4_CLASSIC')
    nc=Dataset(fileout,mode='w',format='NETCDF4')
    # Define Dimensions
    neof_dim = nc.createDimension('neof', neof)
    kmt_dim = nc.createDimension('nlev', kmt)
    jmt_dim = nc.createDimension('jmt', jmt) 
    imt_dim = nc.createDimension('imt', imt) 
    reg_dim = nc.createDimension('nreg', nreg)
    # Title    
    nc.title='EOF per grid point - 3dvar'
    # Create Variables
    # evavar = nc.createVariable('eva', np.single, ('neof','jmt','imt'))
    # evcvar = nc.createVariable('evc', np.single, ('neof','nlev','jmt','imt'))
    evavar = nc.createVariable('eva', np.single, ('neof','jmt','imt'),zlib=True)
    evcvar = nc.createVariable('evc', np.single, ('neof','nlev','jmt','imt'),zlib=True)    
    # Fill in variables
    evavar[:,:,:] = eva
    evcvar[:,:,:,:] = evc    
    nc.close()
    return    

#----------- begin ---------------------

filein='/Users/mario/CMCC/test_mario_MedSea_red/EOF_TSSSH.nc'
# filein='/Users/mario/CMCC/test_mario_MedSea_red/EOF_1set.nc'
fileou='/Users/mario/CMCC/3DVar/test_case_mario_MedSea_red/eofs_fullset.nc'
# fileou='/Users/mario/CMCC/3DVar/test_case_mario_MedSea_red/eofs_1set.nc'
regs,eigval,eigvec = read_nc(filein)

ll_reduce_neof=True
    

imt=regs.shape[1]
jmt=regs.shape[0]
neof=eigval.shape[0]
nlev=eigvec.shape[1]
if (ll_reduce_neof):
    neof=5
    eigval=eigval[0:neof,:]
    eigvec=eigvec[0:neof,:,:]
    
eva=np.zeros((neof,jmt, imt),dtype='float32')
evc=np.zeros((neof,nlev,jmt,imt),dtype='float32')

for reg in range(min(regs.flatten()),max(regs.flatten())):
    y,x = np.where(regs == reg)
    #From array to scalar
    y=y[0]
    x=x[0]
    eva[:,y,x]=eigval[:,reg]
    # switch sla eof from bottom to the top
    evc[:,1:,y,x]=eigvec[:,0:-1,reg]
    evc[:,0,y,x]=eigvec[:,-1,reg]
    

writenc(fileou,eva,evc)
