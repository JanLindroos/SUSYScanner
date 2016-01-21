# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 16:38:17 2015
@author: Jan Lindroos
"""
import scipy as sp
import os,sys

model='models/SUSY/cMSSM_model.py'#model to run over, a model is defined by a set of parameters and a target distribution
run_name='cmssm_gridtest'
fileformat='dat'#Determins fileformat for out files (only dat available at the moment)
merge_files=True#Option to merge dat files (default false)

alg='nonrand'
scan_range={'m0':sp.linspace(200,2000,10),'m12':sp.linspace(200,2000,10),'A0':sp.array([-2500,0]),'sgmu':sp.array([1]),'tanb':sp.array([10,45])}#If data points is to be drawn from file
mode='multiprocessing'#'multiprocessing'#run mode, multiprocessing for single machine, MPI for cluster
chains=5#number of cpus to use in multiprocessing (and for local multiprocessing when MPI is used)
sample_size=sp.inf#Determine maximal number of models in the sample (inf=all grid points)

#Target distrinution
sequential=True#Calculate likelihood sequentially lnP=sum(lnP_i)
lnP_min=-sp.inf#-sp.inf#Minimum allowed

#scan options
batch_size=50#batch size, determines analysis interval for the chains


#Make changes to model
model_change={}
#change parts of SUSY calculations included
model_change['parts']=['softsusy','susyhit','micromegas']#'pythia','delphes','prospino'

#Import constraints from file constraints.py (default is none (only theory errors fail models))
#cpath='constraints.py'#Path to constraints
#[cpath,cmodule]=os.path.split(cpath)
#sys.path.append(os.path.abspath(cpath))
#exec('from '+cmodule.rstrip('.py')+' import constraints')
#model_change['constraints']=constraints


#Print options
print_level=2#0:None,1:State,2:State+accepted, 3:State+accepted+rejected
print_params=["m0","m12","tanb","A0","m_h0","Oh2","LO_tot","LO_ss","LO_sb","LO_gg","LO_nn"]                                                                            
                                                                                                                                     
