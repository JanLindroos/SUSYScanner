# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 15:32:00 2015

@author: Jan Lindroos

Non-rnadom scan
Input is taken from a list of files, with either slha, dat or h5py extension
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 15:54:18 2014

@author: Jan Lindroos
"""
import scipy as sp
import alg_base as alg
import lib.io_lib as io

#Method for turning arrays into a grid of points using cartesian product
def grid_from_arrays(array_dict):
    #fix so numbers appear correctly
    for key in array_dict.keys():
        if isinstance(array_dict[key],(int, long, float, complex)):
            array_dict[key]=sp.array([array_dict[key]])
    
    arrays=array_dict.values()
    arrays = [sp.asarray(a) for a in arrays]
    
    grid_dict=dict.fromkeys(array_dict.keys())

    shape = [len(x) for x in arrays]
    ix = sp.indices(shape, dtype=int)
    ix = ix.reshape(len(arrays), -1).T
    
    for n, arr in enumerate(arrays):
        ix[:, n] = arrays[n][ix[:, n]]
    
    
    keys=grid_dict.keys()
    for i in range(len(keys)):
        grid_dict[keys[i]]=ix[:,i]
    
    return grid_dict

#State of the scan
class state(alg.state):
    
    def __init__(self,rank,opt,model):
        alg.state.__init__(self,rank,opt,model)           
        #full_paramset: None for
        if rank==0:
            #Take input from file            
            if type(opt.scan_range)==str:
                self.full_params_set=io.read_params(opt.scan_range,opt.scan_range.split('.')[-1],model.model_params)
                N=min(opt.sample_size,len(self.full_params_set.values()[0]))
                for key in self.full_params_set.keys():
                    self.full_params_set[key]=self.full_params_set[key][0:N]            
            #If grid is specified
            if type(opt.scan_range)==dict:
                self.full_params_set=grid_from_arrays(opt.scan_range)
            
            #Change stop_size if opt.sample_size>points in file
            self.sample_size=min(len(self.full_params_set.values()[0]),self.sample_size)
                    
        else:
            self.full_params_set=None
    #Update states based on updates from chains
    def update(self,updates):
        #inherit base updates from base class
        alg.state.update(self,updates)

#local state of individual chains        
class chain(alg.chain):

    def __init__(self,rank,opt,model,state):
        alg.chain.__init__(self,rank,opt,model,state)
        self.params_set=dict.fromkeys(state.full_params_set.keys())
        #reduce to subset assigned to chain
        N=len(state.full_params_set.values()[0])
        dn=N/opt.chains
        mod_N=N%dn
        for name in state.full_params_set.keys():
            if rank==opt.chains-1:
                self.params_set[name]=list(state.full_params_set[name][rank*dn:(rank+1)*dn+mod_N])
            else:
                self.params_set[name]=list(state.full_params_set[name][rank*dn:(rank+1)*dn])

        #Remove full param set from state
        state.full_params_set=None
        
        #reduce to subset assigned to chain
        #for name in full_params_set.keys():
        #    if rank==opt.chains-1:
        #        self.params_set[name]=full_params_set[name][rank*dn:(rank+1)*dn+mod_N]
        #    else:
        #        self.params_set[name]=full_params_set[name][rank*dn:(rank+1)*dn]
                
    def update(self):
        alg.chain.update(self)
        #check if all assigned points have been calculated        
        if self.size>=len(self.params_set.values()[0]):                
            self.continue_sampling=False
            self.send=True

#random scan kernel for creating the chain
class kernel(alg.kernel):
    
    def __init__(self,rank,opt):
        alg.kernel.__init__(self,rank,opt)
        
    def initialize(self,state,chain):
        params={}
        #check size of params_set      
        
        for name in chain.params_set.keys():
            params[name]=chain.params_set[name][chain.size]
         
        #print params.keys()
        if 'modelid' in params.keys():
            modelid=str(params.pop('modelid'))
            #print "%i, modelid from file:%s"%(chain.rank,modelid)
        else:
            modelid='%i%01i'%(self.rank,0)+'%i'%chain.accepted
            #print "constructed modelid:%s"%(modelid)
        #print params
        return params, modelid        
                
    #Method for proposing new point
    def propose(self,X_i,state,chain):
        params={}
        for name in chain.params_set.keys():
            params[name]=chain.params_set[name][chain.size]
            
        if 'modelid' in params.keys():
            modelid=str(params.pop('modelid'))
        else:
            modelid='%i%01i'%(self.rank,0)+'%i'%chain.accepted
        
        return params, modelid
        
    #Check wether proposal is accepted    
    def accept(self,X_i,X_f,state,u=sp.rand):
        #Calculate kernel
        acc=True
        #if self.rank==0:
        #    print X_f.lnP,self.lnP_min
        if X_f.lnP<self.lnP_min:
            acc=False
        
        return acc
    
    def weight(self,X):
        X.weight=1
        return
    
    #Does nothing inside loop, all points are properly reweighed after
    #during postprocessing
    def reweight(self,X):
        return
