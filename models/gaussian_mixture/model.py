# -*- coding: utf-8 -*-
"""
Created on Mon Jan 05 17:41:59 2015
@author: Jan Øye Lindroos
Gaussian mixture model
"""
import scipy as sp
import os
from dist_lib import lnP as lnP_fun
import cPickle as pickle

#Model class for Analytic distributions
class model(object):
    #Initialization method to integrate user changes (seen by all chains)
    @classmethod
    def initialize(cls,opt):
        cls.model_params=opt.scan_range.keys()
        cls.param_names=cls.model_params
        #Default is D-dim Gaussian Model with k=1, mu=0 and sigma=0.1
        cls.dist_vars={'dist':'gaussian','mu':[sp.zeros(len(cls.model_params))],'cov':[0.1**2*sp.eye(len(cls.model_params))],'w':[1]}
        cls.dist_vars['vars']=cls.model_params
        
        #list of allowed model changes
        allowed_changes=['dist_vars']
        for change in opt.model_change.keys():
            if change in allowed_changes:
                #Update dicts
                if change in ['dist_vars']:
                    exec('cls.%s.update(opt.model_change[change])'%change)
                        
    #Initialization method to integrate user changes (performed only by master)
    @classmethod
    def master_initialize(cls,opt):
        #store dist_vars used as pickle
        pickle.dump(cls.dist_vars,open(os.path.join(cls.datadirs['main'],'dist_vars.pickle'), 'wb'))
        print 'pickled to %s'%(os.path.join(cls.datadirs['main'],'dist_vars.pickle'))
    
    #Read covariance from file        
    def __init__(self,modelid,params):
        self.modelid=modelid
        self.params=params#All parameters
        self.weight=0
        self.lnP=0
            
    def calculate(self):
        if self.error>0:
            self.accept=False
            self.lnP=-sp.inf
            return
        self.lnP=lnP_fun(self.dist_vars,self.params)          
        return
