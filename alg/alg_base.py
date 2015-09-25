# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 14:32:21 2015
@author: Jan Lindroos
Base classes for algorithms
contains the methods and attributes common to all
algorithm
"""

import scipy as sp

#State of the scan
class state(object):
    
    def __init__(self,rank,opt,model):
        self.param_names=model.param_names#ordered names of parameteres
        self.batch_size=opt.batch_size
        self.continue_sampling=True
        self.lnP_max=-sp.inf*sp.ones(opt.chains)#Empty arrays giving the best fit for each chain
        self.params_bf=[None for i in range(opt.chains)]
        self.sample_size=0
        self.accepted=0
        self.N_prog=opt.batch_size
        self.mean=sp.nan*sp.ones(len(self.param_names))
        self.mean_error=sp.nan*sp.ones(len(self.param_names))
        self.chain_status=sp.ones(opt.chains)
        self.chain_means=[sp.nan*sp.ones(len(self.param_names)) for i in range(opt.chains)]#total mean for each chain
        self.chain_mean_errors=[sp.nan*sp.ones(len(self.param_names)) for i in range(opt.chains)]#mc error in mean
        self.chain_weights=sp.zeros(opt.chains)
        self.chain_sizes=sp.zeros(opt.chains)#Size of chains
        self.chain_accept=sp.zeros(opt.chains)#Size of chains
        self.stop_size=opt.sample_size#Final sample size requested
            
    def update_stats(self):
        self.accepted=sp.sum(self.chain_accept)
        self.sample_size=sp.sum(self.chain_sizes)
        #Update means and errors if any chain contains acceptable points
        #print 'pre state:',self.mean
        if float(sp.sum(self.chain_weights))>0:
            #Total mean
            #for i in range(len(self.chain_weights)):
                
                #print self.chain_weights[i]*self.chain_means[i]
            self.mean=sp.sum([self.chain_weights[i]*self.chain_means[i] for i in range(len(self.chain_weights)) if not sp.isnan(self.chain_means[i]).all()],0)/float(sp.sum(self.chain_weights))               
            #Total error (assumes no correlation at present)
            if sp.greater(self.chain_sizes,self.batch_size).any():
                w=sp.array(self.chain_weights)/float(sp.sum(self.chain_weights))        
                self.mean_error=sp.sqrt(sp.sum([w[i]**2*self.chain_mean_errors[i]**2 for i in range(len(w)) if self.chain_sizes[i]>self.batch_size],0))        

        #print 'post state:',self.mean    
    
    #Update states based on updates from chains
    def update(self,updates):
        #Do manipulation independent of MPI/mp
        #Update the global state based on the info in local states
        #................
        #Quit sampling if total size exceeds specified sample size
        for chain_update in updates:
            self.chain_status[chain_update.rank]=chain_update.continue_sampling
            self.chain_sizes[chain_update.rank]=chain_update.size
            self.chain_accept[chain_update.rank]=chain_update.accepted
            self.chain_means[chain_update.rank]=chain_update.mean
            self.chain_mean_errors[chain_update.rank]=chain_update.mean_error
            self.chain_weights[chain_update.rank]=chain_update.weight
            
            
            #Update global lnP_max and best fit point
            if chain_update.lnP_max>self.lnP_max[chain_update.rank]:
                self.lnP_max[chain_update.rank]=chain_update.lnP_max
                self.params_bf[chain_update.rank]=chain_update.params_bf
    
        #Update global state if any chains have acceptable points
        self.update_stats()    
    
        #print 'Inside update: %i of %i'%(self.sample_size(),self.full_size)
        if self.accepted>=self.stop_size or not self.chain_status.any():
            #print 'stop_sampling'
            self.continue_sampling=False
            
        return    

#local state of individual chains        
class chain(object):

    def __init__(self,rank,opt,model,state):
        #General info needed by chain
        self.rank=rank
        self.param_names=model.param_names
        self.batch_size=opt.batch_size
        
        #Global properties of chain
        self.continue_sampling=True
        self.send=False#wether chain shall send status to master
        self.size=0
        self.weight=0
        self.accepted=0
        self.lnP_max=-sp.inf
        self.params_bf=None
        self.batch_weight=0#Currenr batch weight
        self.batch_sum=sp.zeros(len(self.param_names))#Current batch sum
        self.batch_means=[]#list of batch means
        self.mean=sp.nan*sp.ones(len(self.param_names))
        self.mean_error=sp.nan*sp.ones(len(self.param_names))
        self.stop_size=opt.sample_size#Final sample size requested
        
    #Method for estimating MCMC error
    #NB! This method needs to be checked!
    def update_stats(self,batch_weight,batch_sum):                
        #First time, avoid add to nan error
        #print 'pre:',self.mean 
        if sp.isnan(self.mean).all():
            self.mean=self.batch_sum/float(self.batch_weight)
        else:
            self.mean=(self.weight*self.mean+batch_sum)/float(self.weight+batch_weight)
        #print 'post:',self.mean
        self.batch_means+=[batch_sum/float(batch_weight)]
        self.weight+=batch_weight
        
        if len(self.batch_means)>1:
            self.mean_error=sp.sqrt((self.batch_size/(len(self.batch_means)-1)*sp.sum((self.batch_means-self.mean)**2,0)))/float(sp.sqrt(len(self.batch_means)*self.batch_size))
            
    def update(self):
        #Update at each batch length
        if self.size%self.batch_size==0 and self.size>0:
            if self.batch_weight>0:
                self.update_stats(self.batch_weight,self.batch_sum)
                [self.batch_weight,self.batch_sum]=[0,sp.zeros(len(self.param_names))]#Initialize sums for calculaing batch mean
            self.send=True

        if self.accepted>=self.stop_size:#Final sample size requested
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
        