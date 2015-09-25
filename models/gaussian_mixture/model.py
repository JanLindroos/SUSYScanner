# -*- coding: utf-8 -*-
"""
Created on Mon Jan 05 17:41:59 2015
@author: Jan Øye Lindroos
Gaussian mixture model
"""
import scipy as sp
from scipy.linalg import inv, det
from scipy.misc import logsumexp

#Model baseclass for Gaussian mixture model
#Default is given
class model(object):
    model_params=['x_1','x_2','x_3','x_4','x_5']#Full set of names for model parameters (This is mandatory)
    param_names=model_params#Full set of parameter names (This is what will be written to file)    
    cov=0.1*sp.eye(2)
    mu=sp.zeros(2)    
    
    #Read covariance from file        
    def __init__(self,modelid,params):
        self.modelid=modelid
        self.params=params#All parameters
        self.weight=0
        self.lnP=0
            
    def calculate(self):
        #print self.params.keys(),self.sample_params
        x=sp.array([self.params[key] for key in ['x_1','x_2']])
        self.lnP=-0.5*sp.dot(sp.dot(sp.transpose(self.mu-x),inv(self.cov)),(self.mu-x))            
        return
