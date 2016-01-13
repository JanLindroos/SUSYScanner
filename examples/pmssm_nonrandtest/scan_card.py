# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 16:38:17 2015
@author: Jan Lindroos
"""
import scipy as sp
import os,sys

model='models/SUSY/pMSSM_model.py'#model to run over, a model is defined by a set of parameters and a target distribution
run_name='pMSSM_bino_benchmarks'
fileformat='dat'#Determins fileformat for out files (only dat available at the moment)
merge_files=True#Option to merge dat files (default false)

alg='nonrand'
scan_range='/home/lindroos/ATLAS_comparison/ATLASpMSSM_binoII.dat'#If data points is to be drawn from file
mode='multiprocessing'#'multiprocessing'#run mode, multiprocessing for single machine, MPI for cluster
chains=5#number of chains ran in parallel, for number of threads per chain, see threads model.evgen['threads']
sample_size=20#Determine maximal number of models in the sample

#Target distrinution
sequential=True#Calculate likelihood sequentially lnP=sum(lnP_i)
lnP_min=-sp.inf#-sp.inf#Minimum allowed

#scan options
batch_size=50#batch size, determines analysis interval for the chains

#Make changes to model
model_change={}
#change parts of SUSY calculations included
model_change['parts']=['softsusy','susyhit','higgsbounds','micromegas','pythia','delphes','prospino']
model_change['keep_files']={'slha':True,'mgcard':False,'hepmc':True,'root':True}
model_change['evgen']={'mode':'SUSYMG','nevt':10000,'ecm':8000.,'threads':8,'nlo':1}

#Import constraints from file constraints.py
cpath='constraints.py'#Path to constraints
[cpath,cmodule]=os.path.split(cpath)
sys.path.append(os.path.abspath(cpath))
exec('from '+cmodule.rstrip('.py')+' import constraints')
model_change['constraints']=constraints


#Print options
print_level=2#0:None,1:State,2:State+accepted, 3:State+accepted+rejected
print_params=["Oh2","Bsmumu","a_mu","bsg","m_h0","LO_tot","pr_LO_tot","pr_NLO_tot"]                                                                            
                                                                                                                                     
