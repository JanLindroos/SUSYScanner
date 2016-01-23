# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 14:32:21 2015
@author: Jan Lindroos
Base classes for algorithms
contains the methods and attributes common to all
algorithm
"""

import scipy as sp
import lib.stat_lib as stat
#mu=sum(w*x)/sum(w)=(mu_1*sum(w_1)+mu_2*sum(mu_2))/sum(w)
#State of the scan
class state(object):
    
    def __init__(self,rank,opt,model):
        #General info needed by state
        self.param_names=model.param_names#ordered names of parameteres
        self.sample_size=opt.sample_size#Final sample size requested
        
        #Properties of state       
        self.continue_sampling=True#Global sampling state
        self.lnP_max=-sp.inf#Empty arrays giving the best fit for each chain
        self.size=0#Number of likelihood evaluations
        self.accept=0#Number of accepted models
        self.mean=dict([(key,sp.nan) for key in self.param_names])#Global mean
        self.mc_error=dict([(key,sp.nan) for key in self.param_names])#Standard mc error
        
        #Properties of chains
        self.chains={'continue_sampling':sp.ones(opt.chains)}#Sampling state of chains
        self.chains['size']=sp.zeros(opt.chains)#Number of likelihood evaluations by chain
        self.chains['accept']=sp.zeros(opt.chains)#Number of accepted points by chain
        self.chains['mean']=dict([(key,sp.nan*sp.ones(opt.chains)) for key in self.param_names])#total mean for each chain
        self.chains['mc_error']=dict([(key,sp.nan*sp.ones(opt.chains)) for key in self.param_names])#batch means estimated error http://arxiv.org/abs/math/0601446
        self.chains['weight']=dict([(key,sp.nan*sp.ones(opt.chains)) for key in self.param_names])#total weights of chains
        self.chains['lnP_bf']=-sp.inf*sp.ones(opt.chains)#chain best fit likelihoods
        self.chains['params_bf']=[dict([(key,sp.nan) for key in self.param_names]) for i in range(opt.chains)]#chain best fit models
        self.chain_keys=self.chains.keys()#List of chain data to be sent to state
    
    #Update the global state based on the info in local states (dictionary)
    def update(self,chain_updates):
        #Loop over chains in update
        for chain_update in chain_updates:
            for key in self.chain_keys:
                self.chains[key][chain.rank]=chain_update[key]
            
            #Update global lnP_max and best fit point
            if chain.lnP_max>self.lnP_max:
                self.lnP_bf=chain.lnP_bf
                self.params_bf=chain.params_bf
    
        
        #Update Global statistics
        self.size=sp.sum(self.chains['size'])
        self.accept=sp.sum(self.chains['accept'])
        self.weight=sp.sum(self.chains['weight'])
        for key in self.param_names:
            #combine chain means
            self.mean[key]=sp.sum(self.chains['weight']*self.chains['mean'][key])/float(self.weight)
            #Add mc errors in (quadrature assuming no correlation)
            self.mc_error[key]=sp.sqrt(sp.sum(self.chains['weight']*self.chains['mc_error'])/self.weight)
    
        #check sampling state: terminate if total size exeeds requested sample size or all chains done
        if self.size>=self.sample_size or not self.chains['continue_sampling'].any():
            #print 'stop_sampling'
            self.continue_sampling=False
            
        return    

#local state of individual chains        
class chain(object):

    def __init__(self,rank,opt,model,state):
        #General info needed by chain
        self.rank=rank
        self.param_names=model.param_names
        self.batch_size=opt.batch_size#Number of likelihood evaluations per batch
        self.sample_size=opt.sample_size#Maximum number of likelihood evaluations 
        
        #Local properties of chain
        self.continue_sampling=True#Local state of sampling
        self.send=False#wether chain should send status to master
        self.size=0#Number of Likelihood evaluations
        self.accept=0#Number of accepted models
        self.mean=dict([(key,sp.nan) for key in self.param_names])#chain mean
        self.mc_error=dict([(key,sp.nan) for key in self.param_names])#batch means estimate
        self.weight=0#Sum of weights for full chain
        self.lnP_bf=-sp.inf#Local best fit likelihood
        self.params_bf=None#Local best fit parameter
        self.updates=dict([(key,sp.nan) for key in state.chain_keys])#Attributes to be sent to state

        #Properties of current batch
        self.batch={'data':dict([(key,sp.array([])) for key in self.param_names])}#data in current batch (only models with non-zero weight stored (accepted))
        self.batch['data']['weight']=sp.array([])#weights of points in current batch (needed for statistics)
        self.batch['mean']=dict([(key,sp.array([])) for key in self.param_names])#batch means
        self.batch['weight']=dict([(key,sp.array([])) for key in self.param_names])#batch weights
     
    #Update chain info (every loop)       
    def update(self,X):
        #Initialize batch after update
        if self.size%self.batch_size==0 and self.size>0:
            self.send=False
            self.updates=dict([(key,sp.nan) for key in state.chain_keys])
            self.batch['params']=dict([(key,sp.array([])) for key in self.param_names])#reset dictionary
            self.batch['weight']=sp.array([])
         
        #Update number of likelihood evaluations
        self.size+=1
        #if accepted update accept
        if X.accept:
            self.accept+=1
            #check best fit point
            if X.lnP>self.lnP_bf:
                self.lnP_bf=X.lnP
                self.params_bf=X.params
        
            #add accepted points to batch data (rejected models always have 0 weight)
            sp.append(self.batch['data']['weight'],X.weight)
            for key in self.batch.keys():
                sp.append(self.batch['data']['params'][key],X.params[key])
                    
        #Update at each batch length
        if self.size%self.batch_size==0 and self.size>0:
            #only send if batch check for accepted points in previous batch
            if len(self.batch['data']['weight'])>0:
                #add batch weight to array of batch weights
                sp.append(self.batch['weight'],sp.sum(self.batch['data']['weight']))
                #add batch mean to array of batch means
                for key in self.param_names:
                    sp.append(self.batch['mean'][key],stat.weighted_mean(self.batch['data'][key],self.batch['data']['weight']))
                    #calculate curren chain mean
                    self.mean[key]=stat.weighted_mean(self.batch['mean'],self.batch['weight'])
                    #Calculate current chain mc_error
                    self.mc_error[key]=stat.mc_error(self.mean,self.batch['mean'][key],self.batch['weight'])
                
                #Fill updates for sending
                self.updates=dict([(key,getattr(self,key)) for key in state.chain_keys])
                #Set send status to true    
                self.send=True

        if self.size>=self.sample_size:#Final sample size requested
            print 'Worker %i has reached max size, sampling stopped'%(self.rank)
            self.continue_sampling=False
            self.send=True

#random scan kernel for creating the chain
class kernel(object):
    
    def __init__(self,rank,opt):
        self.rank=rank
        self.lnP_min=opt.lnP_min
        #initialize weights for multi-ranged scan parameters
        self.range_weight={}
        