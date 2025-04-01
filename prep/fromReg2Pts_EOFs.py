#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 12:24:29 2024

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

def writenc(fileout,eva,evc,regs):
    imt=regs.shape[1]
    jmt=regs.shape[0]
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
    evavar = nc.createVariable('eva', np.single, ('neof','jmt','imt'),zlib=True)
    evcvar = nc.createVariable('evc', np.single, ('neof','nlev','jmt','imt'),zlib=True)    
    regvar = nc.createVariable('regs', np.single,('jmt','imt'),zlib=True)    
    # Fill in variables
    evavar[:,:,:]   = eva
    evcvar[:,:,:,:] = evc    
    regvar[:,:] = regs   
    nc.close()
    return    

#----------- begin ---------------------

filein='/work/cmcc/ma20223/MED/test_case_med/EOF.nc'
fileou='/work/cmcc/ma20223/MED/test_case_med/peofs.nc'
regs,eigval,eigvec = read_nc(filein)


imt=regs.shape[1]
jmt=regs.shape[0]
neof=eigval.shape[0]
nlev=eigvec.shape[1]

regs_new = np.zeros((jmt,imt))
eva = np.zeros((neof,jmt,imt))
evc = np.zeros((neof,nlev,jmt,imt))
tmp = np.zeros((eigvec.shape))
# switch sla eof from bottom to the top
tmp[:,1:,:]=eigvec[:,0:-1,:]
tmp[:,0,:]=eigvec[:,-1,:]

k=1

for j in range(jmt):
   for i in range(imt):
     regs_new[j,i] = k
     k = k + 1
     if regs[j,i] == 0:
         evc[:,:,j,i] = 0.0
         eva[:,j,i]   = 0.0
     else:
         eva[:,j,i]   = eigval[:,regs[j,i]]
         evc[:,:,j,i] = tmp[:,:,regs[j,i]]

writenc(fileou,eva,evc,regs_new)
