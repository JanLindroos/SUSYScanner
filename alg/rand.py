# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 15:54:18 2014

@author: Jan Lindroos
"""
import scipy as sp
import alg_base as alg
#State of the scan
class state(alg.state):
    
    def __init__(self,rank,opt,model):
        alg.state.__init__(self,rank,opt,model)    
        
    #Update states based on updates from chains
    def update(self,updates):
        #inherit base updates from base class
        alg.state.update(self,updates)

#local state of individual chains        
class chain(alg.chain):

    def __init__(self,rank,opt,model,state):
        alg.chain.__init__(self,rank,opt,model,state)

#random scan kernel for creating the chain
class kernel(alg.kernel):
    
    def __init__(self,rank,opt,model):
        alg.kernel.__init__(self,rank,opt,model)
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
    #must be run after initialize X_i only used for
    #initialization in rand scan and state is just a dummy
    def propose(self,X_i,state,chain):
        params={}
        for key in self.scan_range.keys():
            #sample multiranged parameter
            if key in self.range_weight.keys():
                i_sel=2*sp.searchsorted(sp.cumsum(self.range_weight[key]),sp.rand())
                params[key]=sp.rand()*(self.scan_range[key][i_sel+1]-self.scan_range[key][i_sel])+self.scan_range[key][i_sel]
            else:
                params[key]=sp.rand()*(self.scan_range[key][1]-self.scan_range[key][0])+self.scan_range[key][0]
        
        #Add constant parameters
        for key in self.constants.keys(): 
            params[key]=self.constants[key]
        
        #functions
        for key in self.functions.keys():
            params[key]=self.functions[key](params)
            
        modelid='%i%01i'%(self.rank,0)+'%i'%chain.accept
            
        return params, modelid
        
    #Check wether proposal is accepted    
    def accept(self,X_f,X_i,state,u=sp.rand):
        #Calculate kernel
        acc=True
        #if self.rank==0:
        #    print X_f.lnP,self.lnP_min
        if X_f.lnP<self.lnP_min:
            acc=False
        
        return acc
    
    def weight(self,X):
        X.weight=sp.exp(X.lnP)
        return X
    
    #Does nothing inside loop, all points are properly reweighed after
    #during postprocessing
    def reweight(self,X):
        return X
        
    