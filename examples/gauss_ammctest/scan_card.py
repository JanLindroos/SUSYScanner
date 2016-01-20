# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 16:38:17 2015
@author: Jan Lindroos
Default Scan Card
"""
import scipy as sp
import sys
import cPickle as pickle

#Main options******************************************************************
alg='ammc'#algorithm
mode='multiprocessing'#run mode, multiprocessing for single machine, MPI for cluster
chains=10#number of cpus to use in multiprocessing (and for local multiprocessing when MPI is used)
sample_size=1000000#Number of likelihood evaluations (after initialization)
fileformat='dat'#file format for storing data
merge_files=True#Wheter to merge the output files from the different chains
run_name='ammc_test'
#******************************************************************************

#Model options*****************************************************************
scan_range=dict([('x_%i'%i,[-10,10]) for i in range(3)])#Sets the range for the scan
model='models/gaussian_mixture/model.py'#model to run over

#Generate gaussian mixture
# sys.path.append('models/gaussian_mixture')
# from dist_lib import gen_dist
# #create test distribution
# dist='gaussian'
# k=3#number of gaussians,
# sigma=0.1#Normalized variance, float or array with dim of parameter space
# dist_vars=gen_dist('gaussian',scan_range,k,sigma,sym=False,uncorr=False)

#Get distribution from file
dist_path='dist_vars.pickle'
dist_vars=pickle.load(open(dist_path,'rb'))

model_change={'dist_vars':dist_vars}
#******************************************************************************

#Alg Options*******************************************************************
sequential=False#Target sequentially lnP=sum(lnP_i)
lnP_min=-10#-sp.inf#Minimum allowed
batch_size=1000#batch size, determines analysis interval for the chains

#mcmc specific options
mcmc={}
mcmc['input_data']='path_to_file'#Path to input data for constructing proposal
#******************************************************************************

#Print Options*****************************************************************
print_level=1
#******************************************************************************