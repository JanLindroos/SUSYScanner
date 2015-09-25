# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 17:53:05 2015
@author: Jan Lindroos

Simple script to generate input dat file for gridded scan
"""
import scipy as sp

params={}
filename='input_grid.dat'#output file file

params['m0']=sp.linspace(200,3000,10)
params['m12']=sp.linspace(200,3000,10)
params['A0']=sp.linspace(-5000,5000,10)
params['tanb']=sp.linspace(2,60,10)
params['mu']=sp.array([1])

def grid_from_arrays(array_dict):
    grid_dict=dict.fromkeys(array_dict.keys())
    arrays=array_dict.values()
    arrays = [sp.asarray(a) for a in arrays]
    shape = (len(x) for x in arrays)
    ix = sp.indices(shape, dtype=int)
    ix = ix.reshape(len(arrays), -1).T
    
    for n, arr in enumerate(arrays):
        ix[:, n] = arrays[n][ix[:, n]]
    
    
    keys=grid_dict.keys()
    for i in range(len(keys)):
        grid_dict[keys[i]]=ix[:,i]
    
    return grid_dict
    
