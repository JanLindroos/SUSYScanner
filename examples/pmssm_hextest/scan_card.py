# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 16:38:17 2015
@author: Jan Lindroos
Default Scan Card
"""
import scipy as sp
import os,sys


run_name='test_run'
model='models/SUSY/pMSSM_model.py'#model to run over, a model is defined by a set of parameters and a target distribution
fileformat='dat'
alg='rand'
mode='MPI'#'multiprocessing'#run mode, multiprocessing for single machine, MPI for cluster
chains=128#number of cpus to use in multiprocessing (and for local multiprocessing when MPI is used)
sample_size=1000000#requested mcmc sample size (this is the number of accepted points)

#define Parameters
#Scan parameters
scan_range={}
scan_range['tanb']=[2,60]
for key in ['M_1',"M_2","Ab","Atau","mu"]:
    scan_range[key]=[-4000,4000]
scan_range['At']=[-8000,8000]
for key in ["M_3","mA","meL","mtauL","meR","mtauR","mqL1","mqL3","muR","mtR","mdR","mbR"]:
    scan_range[key]=[100,4000]

#define constants
constants={}

#define parameters with functional dependence
functions={}
functions['mmuL']=lambda pars: pars['meL']
functions['mmuR']=lambda pars: pars['meR']
functions['mqL2']=lambda pars: pars['mqL1']
functions['mcR']=lambda pars: pars['muR']
functions['msR']=lambda pars: pars['mdR']


#Target distrinution
sequential=False#use sequentially lnP=sum(lnP_i) only calculate while constraints are fullfilled
lnP_min=-1e-10#Minimum allowed

#scan options
batch_size=10000#batch size (accepted and rejected), determines analysis interval for the chains

#Make changes to model
model_change={}
#change parts of SUSY calculations included
model_change['parts']=['softsusy','susyhit','higgsbounds','micromegas']
#Import constraints from file constraints.py
cpath='constraints_cut.py'#Path to constraints
[cpath,cmodule]=os.path.split(cpath)
sys.path.append(os.path.abspath(cpath))
exec('from '+cmodule.rstrip('.py')+' import constraints')
model_change['constraints']=constraints

print_level=2
print_params=['HB_excl','LEP','Gamma_Z','m_h0','Oh2','bsg','Bsmumu','drho','Oh2_wino']                                                                            
                                                                                                                                     
