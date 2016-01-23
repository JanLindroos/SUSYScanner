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
        #General info needed by state
        self.param_names=model.param_names#ordered names of parameteres
        self.batch_size=opt.batch_size
        self.sample_size=opt.sample_size#Final sample size requested
        
        #Properties of state       
        self.continue_sampling=True
        self.lnP_max=-sp.inf#Empty arrays giving the best fit for each chain
        self.size=0#Number of likelihood evaluations
        self.accepted=0#Number of accepted models
        self.mean=dict([(key,sp.nan) for key in self.params_names])#Global mean
        self.mc_error=dict([(key,sp.nan) for key in self.params_names])#Standard deviation between chain means
        self.weight=0#Sum of all weights
        
        #Properties of chains
        self.chains['continue_sampling']=sp.ones(opt.chains)#Status of the different chains
        self.chains['lnP_max']=-sp.inf*sp.ones(opt.chains)#chain best fit likelihoods
        self.chains['params_bf']=[dict([(key,sp.nan) for key in self.params_names]) for i in range(opt.chains)]#chain best fit models
        self.chains['mean']=dict([(key,sp.nan*sp.ones(opt.chains)) for key in self.params_names])#total mean for each chain
        self.chains['mc_error']=dict([(key,sp.nan*sp.ones(opt.chains)) for key in self.params_names])#batch means estimated error http://arxiv.org/abs/math/0601446
        self.chains['weight']=sp.zeros(opt.chains)#total weight of chains
        self.chains['size']=sp.zeros(opt.chains)#Number of likelihood evaluations
        self.chains['accept']=sp.zeros(opt.chains)#Number of accepted points in chain
        
            
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
        for chain in updates:
            self.chains['continue_sampling'][chain.rank]=chain.continue_sampling
            self.chains['size'][chain.rank]=chain.size
            self.chains['accept'][chain.rank]=chain.accepted
            self.chains['mean'][chain.rank]=chain.mean
            self.chains['mc_error'][chain.rank]=chain.mc_error#batch meansestimate http://arxiv.org/abs/math/0601446
            self.chains['weight'][chain.rank]=chain.weight
            self.chains['lnP_max'][chain.rank]=chain.lnP_max
            self.chains['params_bf'][chain.rank]=chain.params_bf
            
            #Update global lnP_max and best fit point
            if chain.lnP_max>self.lnP_max:
                self.lnP_max=chain.lnP_max
                self.params_bf=chain.params_bf
    
        
        #Update Global statistics
        self.size=sp.sum(self.chains['size'])
        self.accepted=sp.sum(self.chains['size'])
        for key in self.param_names:
            self.mean[key]=sp.sum(self.chains['weight']*self.chains['mean'][key])/sp.sum(self.chains['mean'][key])
            self.std[key]=sp.sum()
    
        #print 'Inside update: %i of %i'%(self.sample_size(),self.full_size)
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
        self.continue_sampling=True
        self.send=False#wether chain shall send status to master
        self.size=0#Number of Likelihood evaluations
        self.weight=0#Sum of weights for full chain
        self.accepted=0#Number of accepted models
        self.lnP_max=-sp.inf#Local best fit likelihood
        self.params_bf=None#Local best fit parameter
        self.mean=dict([(key,sp.nan) for key in self.params_names])#chain mean
        self.mc_error=dict([(key,sp.nan) for key in self.params_names])#batch means estimate http://arxiv.org/abs/math/0601446
        #Properties of current batch
        self.batch={'data':dict([(key,sp.array([])) for key in self.param_names])}#batch data (only accepted points are stored)
        self.batch['data']['weight']=sp.array([])#weights of data points
        self.batch['means']=dict([(key,sp.array([])) for key in self.param_names])#list of batch means used to find mc_error
        self.batch['batch_weight']=dict([(key,sp.array([])) for key in self.param_names])#weights corresponding to params
        
        
        def mean        
        
        #calculate mc_error
        def mc_error(mu,weights):
            mu_tot=mean(mu,weights)
            V_1=sp.sum(weights)
            V_2=sp.sum(weights**2)
            sigma_2=V_1/float(V_1**2-V_2)*sp.sum(weights*(mu_b-mu_tot)**2)
            #print sp.sum(weights), len(weights)
            mc_error=sp.sqrt(sigma_2/float(len(weights)))
            return mc_error,mu_tot
    #Method for estimating MCMC error
    #NB! This method needs to be checked!
#    def update_stats(self,batch_weight,batch_sum):                
#        #First time, avoid add to nan error
#        #print 'pre:',self.mean 
#        if sp.isnan(self.mean).all():
#            self.mean=self.batch_sum/float(self.batch_weight)
#        else:
#            self.mean=(self.weight*self.mean+batch_sum)/float(self.weight+batch_weight)
#        #print 'post:',self.mean
#        self.batch_means+=[batch_sum/float(batch_weight)]
#        self.weight+=batch_weight
#        
#        if len(self.batch_means)>1:
#            self.mean_error=sp.sqrt((self.batch_size/(len(self.batch_means)-1)*sp.sum((self.batch_means-self.mean)**2,0)))/float(sp.sqrt(len(self.batch_means)*self.batch_size))
            
    def update(self,X):
        #Initialize batch after update
        if self.size%self.batch_size==0 and self.size>0:
            self.send=False
            self.batch['params']=dict([(key,sp.array([])) for key in self.param_names])#reset dictionary
            self.batch['weight']=sp.array([])
         
        #Update number of likelihood evaluations
        self.size+=1
        
        #Only add accepted points to batch data (rejected models always have 0 weight)
        if X.accept:
            sp.append(self.batch['weight'],X.weight)
            for key in self.batch.keys():
                sp.append(self.batch['params'][key],X.params[key])
                sp.append(self.batch['weight'],X.params[key])
                    
        #Update at each batch length
        if self.size%self.batch_size==0 and self.size>0:
            #check for accepted points in previous batch
            if len(self.batch.values()[0])>0:
                for key in self.params_names:
                    self.batch['weight'][key]=self.batch['weight']*self.batch['data'][key]/float(sp.sum(self.batch['weight']))
                    
                self.send=True#Set send status to true
                
                
                
            
            if self.batch_weight>0:
                self.update_stats(self.batch_weight,self.batch_sum)
                [self.batch_weight,self.batch_sum]=[0,sp.zeros(len(self.param_names))]#Initialize sums for calculaing batch mean
            #Send results to master
            

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
        