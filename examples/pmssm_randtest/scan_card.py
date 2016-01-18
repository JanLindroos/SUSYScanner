# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 16:38:17 2015
@author: Jan Lindroos
Default Scan Card
"""
import scipy as sp
import os,sys


run_name='ATLAS_repro_100models'
model='models/SUSY/pMSSM_model.py'#model to run over, a model is defined by a set of parameters and a target distribution
fileformat='dat'
alg='rand'
mode='multiprocessing'#'multiprocessing'#run mode, multiprocessing for single machine, MPI for cluster
chains=10#number of cpus to use in multiprocessing (and for local multiprocessing when MPI is used)

sample_size=100#requested sample size

#define Parameters
#Scan parameters
scan_range={}
scan_range['tanb']=[1,60]
scan_range['mu']=[-4000,-80,80,4000]
scan_range['M_2']=[-4000,-70,70,4000]
scan_range['At']=[-8000,8000]
for key in ["M_3","mqL1","muR","mdR"]:
    scan_range[key]=[200,4000]   
for key in ["mA","mqL3","mtR","mbR"]:
    scan_range[key]=[100,4000]
for key in ["mtauR","mtauL","meL","meR"]:
    scan_range[key]=[90,4000]
for key in ['M_1',"Ab","Atau"]:
    scan_range[key]=[-4000,4000]

#define parameters with functional dependence
functions={}
functions['mqL2']=lambda pars: pars['mqL1']
functions['mcR']=lambda pars: pars['muR']
functions['msR']=lambda pars: pars['mdR']
functions['mmuL']=lambda pars: pars['meL']
functions['mmuR']=lambda pars: pars['meR']

#define constants
constants={}

#Target distrinution
sequential=False#use sequentially lnP=sum(lnP_i) only calculate while constraints are fullfilled
lnP_min=-1e10#Minimum allowed (-sp.inf means accept everything)

#scan options
batch_size=100#batch size, determines analysis interval for the chains

#Make changes to model*********************************************************
model_change={}
#change parts of SUSY calculations included
model_change['parts']=['softsusy','susyhit','higgsbounds','micromegas']
#Import constraints from file constraints.py
cpath='ATLAS_pMSSM_constraints.py'#Path to constraints
[cpath,cmodule]=os.path.split(cpath)
sys.path.append(os.path.abspath(cpath))
exec('from '+cmodule.rstrip('.py')+' import constraints')
model_change['constraints']=constraints
#LHC EvGen
model_change['evgen_ecm']=14000# [GeV] LHC center of mass energy (default 8000)
model_change['evgen_nevt']=10000# Number of LHC events to generate (default 10000)
model_change['evgen_threads']=10#number of threads for evt generation (default 1, recommended chains*threads<=cpus available)
#model_change['Analysis']=#Still need to figure out what these are
#******************************************************************************

#Print options
print_level=2#0:None,1:State,2:State+accepted, 3:State+accepted+rejected
print_params=['HB_excl','LEP','Gamma_Z','m_h0','Oh2','bsg','Bsmumu','drho']                                                                            
                                                                                                                                     
