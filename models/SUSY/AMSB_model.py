# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 18:42:12 2015

@author: Jan Lindroos

pMSSM Model
"""

import base_model as base
import shutil,os
import scipy as sp
from susy_lib import init_slha
    
class model(base.model):
    model_params=["m0","m32","tanb","sgmu","cQ","cuR","cdR","cL","ceR","cHu","cHd"]
    
    param_names=model_params+base.model.obs_params

    @classmethod
    def initialize(cls,opt):
        super(model, cls).initialize(opt)
        cls.param_names=cls.model_params+cls.obs_params
    
    def __init__(self,modelid,params):
        base.model.__init__(self,modelid,params)
        model.slha_tmpl='amsb_slha.tmpl'
        #Initialize model parameters
        self.params=dict([(name,sp.nan) for name in self.param_names])
        self.params.update(params)#All parameters
        init_slha(self)
            