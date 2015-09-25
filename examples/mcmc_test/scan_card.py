# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 16:38:17 2015
@author: Jan Lindroos
Default Scan Card
"""
import scipy as sp

run_name='test_run'
model='gaussian_mixture'#model to run over, a model is defined by a set of parameters and a target distribution
fileformat='dat'
alg='mcmc'
mode='MPI'#run mode, multiprocessing for single machine, MPI for cluster
cpus=20#number of cpus to use in multiprocessing (and for local multiprocessing when MPI is used)
sample_size=1000000#requested mcmc sample size

#define Parameters
#Scan parameters
scan_range={'x_1':[-10,10]}
scan_range['x_2']=[-10,10]#[-5,5,2,3], allow for multiple regions [-2,1] and [2,3]
#define constants
constants={'x_3':5}
constants['x_4']=34.1
#define parameters with functional dependence
functions={'x_5':lambda params: params['x_1']**2}

#Target distrinution
sequential=True#WTarget sequentially lnP=sum(lnP_i)
parts=['part_1','part_2']#Provide steps of target components to include
lnP=None#Provide a path to custom target distribution (default none: use model default)
lnP_min=-5#-sp.inf#Minimum allowed

#scan options
batch_size=10000#batch size, determines analysis interval for the chains

