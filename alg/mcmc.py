# -*- coding: utf-8 -*-
"""
Created on Sat Jan 24 21:03:02 2015
@author: Jan Lindroos

Basic metropolis-hastings algorithm
"""

import scipy as sp
import alg_base as alg
try:
    from scipy.linalg import det, inv
except:
    from numpy.linalg import det, inv    
#State of the scan
class state(alg.state):
    
    def __init__(self,rank,opt,model):
        alg.state.__init__(self,rank,opt,model)
        if 'Q_cov' in dir(opt):
            self.Q_cov=opt.Q_cov
        else:
            self.Q_cov=0.1*sp.eye(len(opt.scan_range))
        
    #Update states based on updates from chains
    def update(self,updates):
        alg.state.update(self,updates)
    
    #Gaussian Proposal distribution    
    def lnQ(self,X_i,X_f):
        mu=sp.array([X_i.params[key] for key in X_i.model_params])
        x=sp.array([X_f.params[key] for key in X_f.model_params])
        return -0.5*sp.dot(sp.dot(sp.transpose(mu-x),inv(self.Q_cov)),(mu-x))
     
    #New sample point bsed on intitial point
    def sample_Q(self,mu):
        return sp.random.multivariate_normal(mu,self.Q_cov,1)[0]

#local state of individual chains        
class chain(alg.chain):

    def __init__(self,rank,opt,model,state):
        alg.chain.__init__(self,rank,opt,model,state)

#random scan kernel for creating the chain
class kernel(alg.kernel):
    
    def __init__(self,rank,opt):
        alg.kernel.__init__(self,rank,opt)
        self.scan_range=opt.scan_range
        self.constants=opt.constants
        self.functions=opt.functions
        
    #Method for sampling first point scan range uniformely (default)
    def initialize(self,state,chain):
        params={}
        for key in self.scan_range.keys():
            #Check for single range
            if len(self.scan_range[key])==2:
                params[key]=sp.rand()*(self.scan_range[key][1]-self.scan_range[key][0])+self.scan_range[key][0]
            else:
                #calculate weights of sub_regions
                sub_size=sp.array([])
                #Determine weights of region
                for i in range(0, len(self.scan_range[key]), 2):
                    sub_size=sp.append(sub_size,self.scan_range[key][i+1]-self.scan_range[key][i])
                    self.range_weight[key]=sub_size/float(sp.sum(sub_size))
                
                #sample region based on size
                i_sel=2*sp.searchsorted(sp.cumsum(self.range_weight[key]),sp.rand())
                #sample point
                params[key]=sp.rand()*(self.scan_range[key][i_sel+1]-self.scan_range[key][i_sel])+self.scan_range[key][i_sel]
                
        #params=dict([(key,sp.rand()*(self.scan_range[key][1]-self.scan_range[key][0])+self.scan_range[key][0]) for key in self.scan_range.keys() if type(self.scan_range[key])==list])
        
        #Add constant parameters
        for key in self.constants.keys():
            params[key]=self.constants[key]
            
        for key in self.functions.keys():
            params[key]=self.functions[key](params)
            
        modelid='%i%01i'%(self.rank,0)+'%i'%chain.accept
                
        return params,modelid
                
    #Method for proposing new point X_f given X_i
    def propose(self,X,state,chain):
        #Initialize new parameters
        mu=sp.array([X.params[key] for key in X.model_params])
        x=state.sample_Q(mu)
        params=dict((X.model_params[i],x[i]) for i in range(len(X.model_params)))
        #Add constants        
        params.update(self.constants)
        #Add functions
        for key in self.functions.keys():
            params[key]=self.functions[key](params)
            
        modelid='%i%01i'%(self.rank,0)+'%i'%chain.accept
        
        return params,modelid
        
    #Check wether proposal is accepted    
    def accept(self,X_i,X_f,state,u=sp.rand):
        #Calculate kernel
        a=X_f.lnP-X_i.lnP-state.lnQ(X_i,X_f)+state.lnQ(X_f,X_i)
        if sp.log(u)<=a:
            acc=True
        else:
            acc=False
        
        return acc
    
    def weight(self,X):
        X.weight=1
        return
    
    #Does nothing inside loop, all points are properly reweighed after
    #during postprocessing
    def reweight(self,X):
        X.weight+=1
        return
