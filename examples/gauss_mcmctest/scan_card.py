# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 16:38:17 2015
@author: Jan Lindroos
Default Scan Card
"""
import scipy as sp
import sys
import cPickle as pickle

alg='mcmc'
mode='multiprocessing'#run mode, multiprocessing for single machine, MPI for cluster
chains=10#number of cpus to use in multiprocessing (and for local multiprocessing when MPI is used)
sample_size=1000000#Number of likelihood evaluations (after initialization)
fileformat='dat'
merge_files=True

scan_range=dict([('x_%i'%i,[-10,10]) for i in range(3)])
model='models/gaussian_mixture/model.py'#model to run over, a model is defined by a set of parameters and a target distribution

#Generate gaussian mixture
# sys.path.append('models/gaussian_mixture')
# from dist_lib import gen_dist
# #create test distribution
# dist='gaussian'
# k=3#number of gaussians,
# sigma=0.1#Normalized variance, float or array with dim of parameter space
# dist_vars=gen_dist('gaussian',scan_range,k,sigma,sym=False,uncorr=False)

dist_path='/home/lindroos/SUSYScanner/SUSY_runs/gauss_randtest/gaussian_rand_D3_k3_N1000000/dist_vars.pickle'
dist_vars=pickle.load(open(dist_path,'rb'))
model_change={'dist_vars':dist_vars}

#Name of run (folder for output)
run_name='%s_%s_D%i_k%i_N%i'%(dist_vars['dist'],alg,len(scan_range.keys()),len(dist_vars['w']),sample_size)

#Target distribution
sequential=False#Target sequentially lnP=sum(lnP_i)
lnP_min=-10#-sp.inf#Minimum allowed

#scan options
batch_size=1000#batch size, determines analysis interval for the chains
#print options
print_level=1