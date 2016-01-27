# -*- coding: utf-8 -*-
"""
Created on Sat Jan 24 21:03:02 2015
@author: Jan Lindroos

Addaptive mcmc algorithm
"""

import lib.cluster_lib as cl
import lib.stat_lib as stat
import lib.io_lib as io 
import scipy as sp
import alg_base as alg
from scipy.linalg import det, inv

   
#State of the scan
class state(alg.state):
    
    def __init__(self,rank,opt,model):
        alg.state.__init__(self,rank,opt,model)
        #Import mcmc options
        self.mcmc={'input_data':None}
        self.mcmc['k']=opt.chains#Number of clusters
        self.mcmc['alpha_Q']=sp.ones(self.mcmc['k'])#Covariance scaling factors
        self.mcmc['beta_Q']=0.1#fraction of global sampling
        if 'mcmc' in dir(opt):
            self.mcmc.update(opt.mcmc)
        self.mcmc['n_acc']=sp.zeros(self.mcmc['k'])
        self.mcmc['n_try']=sp.zeros(self.mcmc['k'])
        #Import data for clustering
        if self.mcmc['input_data']==None:
            clusters={'norm':None,'i_p':None,'weight':None,'mean':None,'cov':None,'size':None}
            pass
        else:
            cluster_data=io.import_data(self.mcmc['input_data'],self.mcmc['params']+[self.mcmc['lnw']],norm=False)
            #Get weight
            lnw=cluster_data.pop(self.mcmc['lnw'])
            clusters=cl.cluster(cluster_data,self.mcmc['k'],lnw=lnw)
        
        self.mcmc['clusters']=clusters
        
    #Update states based on updates from chains
    def update(self,chain_updates):
        alg.state.update(self,chain_updates)
        for chain_update in chain_updates:
            #update global trial/error count 
            state.mcmc['n_try']+=chain.batch['n_try']
            state.mcmc['n_acc']+=chain.batch['n_acc']
            #get cluster data from chains
            #Need to figure out how to do this in practise...Maybe do all updates through chain.update?
            cluster_data=chain.batch['data']['params']
        
        for i in range(state.mcmc['k']):
            if state.mcmc['n_try'][i]>=state.batch_size:
                r=state.mcmc['n_acc'][i]/float(state.mcmc['n_try'][i])
                if r<0.1:
                    sigma=0.04659906017846561
                    scale=sp.exp(-(0.1-r)**2/float(2*sigma**2))
                if r>0.3:
                    a=18.367346938775512
                    scale=a*(r-0.3)**2+1

                state.mcmc['alpha_Q'][i]=scale*state.mcmc['alpha_Q'][i]

            lnw=cluster_data.pop('lnP')
            #Update proposal based on chain data
            clusters=cl.cluster(cluster_data,self.mcmc['k'],lnw=lnw,centroids=self.clusters['mean'],norm=False)
            #Update proposal distribution
            if self.mcmc['clusters']['weight']==None:
                self.mcmc['clusters']=clusters
            else:
                #Go through clusters and update covariances and means
                for i in range(self.mcmc['k']):
                    means=[self.mcmc['clusters']['mean'][i],clusters['mean'][i]]
                    covs=[self.mcmc['clusters']['cov'][i],clusters['cov'][i]]
                    sizes=[self.mcmc['clusters']['size']*self.mcmc['clusters']['weight'][i],clusters['size']*clusters['weight'][i]]
                    self.mcmc['clusters']['mean'][i]=stat.combine_means(means,sizes)
                    self.mcmc['clusters']['cov'][i]=stat.combine_covs(covs,means,sizes)
                    self.mcmc['clusters']['weight'][i]=(sizes[0]+sizes[1])/float(self.mcmc['clusters']['size']+clusters['size'])
                
                self.mcmc['clusters']['size']+=clusters['size']

#local state of individual chains        
class chain(alg.chain):

    def __init__(self,rank,opt,model,state):
        alg.chain.__init__(self,rank,opt,model,state)
        self.cind=None#Initialize current cluster index
        self.batch['n_acc']=sp.zeros(state.mcmc['k'])
        self.batch['n_try']=sp.zeros(state.mcmc['k'])
        
    #Update states based on updates from chains
    def update(self,X):
        alg.chain.update(self,X)
        
        #Update n_acc,n_try

#random scan kernel for creating the chain
class kernel(alg.kernel):
    
    def __init__(self,rank,opt,model):
        alg.kernel.__init__(self,rank,opt,model)
        self.scan_range=opt.scan_range
        self.constants=opt.constants
        self.functions=opt.functions

    #Gaussian Proposal distribution    
    def lnQ(self,X_i,X_f,state):
        x_i,x_f=[sp.zeros(range(len(state.mcmc['i_p'].keys()))),sp.zeros(range(len(state.mcmc['i_p'].keys())))]
        for key in state.mcmc['i_p'].keys():
            x_i[clusters['i_p'][key]]=X_i.params[key]
            x_f[clusters['i_p'][key]]=X_f.params[key]
            
        cind=cl.nearest_centroid(x_i,mcmc['mean'])
        
        #Add local contribution
        Q_k=[-0.5*sp.dot(sp.dot(sp.transpose(x_i-x_f),inv(state.mcmc['alpha_Q'][cind]*state.mcmc['cov'][cind])),(x_i-x_f))]
        Q_w=[(1-state.mcmc['beta_Q'])*sp.array([1/sp.sqrt(det(state.mcmc['alpha_Q'][cind]*state.mcmc['cov'][cind]))])]
        
        #Add non-local uniform contribution (Volume normalized to 1)
        Q_k+=[0]
        Q_w+=1
        
        #Add non-local gaussian contribution
        for i in range(state.mcmc['k']):
            #exclude current cluster from large jump
            if i!=cind:
                Q_k+=[-0.5*sp.dot(sp.dot(sp.transpose(state.mcmc['mean'][i]-x),inv(state.mcmc['cov'][i])),(state.mcmc['mean'][i]-x))]
                Q_w+=[state.mcmc['beta_Q']*state.mcmc['weight']/float(sp.sqrt(det(inv(state.mcmc['cov'][i]))))]

        lnQ=logsumexp(Q_k,b=Q_w)
        return lnQ
     
    #New sample point based on intitial point
    def sample_Q(self,state,x_i=None,cind=None,):
        #Generate random number to determine wether to do local or global step
        u=sp.rand()
        if u<state.mcmc['beta_Q'] or x_i==None:
            #Generate completely random step
            u=sp.rand()
            if u<state.mcmc['beta_Q']:
                x_f=(mcmc['range'][1]-mcmc['range'][0])*sp.rand()+state.mcmc['range'][0]
            #otherwise sample from mixture
            else:
                cind=cl.nearest_centroid(x_i,state.mcmc['centroids'])
                x_f=sp.random.multivariate_normal(state.mcmc['mean'],state.mcmc['cov'][cind],1)[0]
        #do local step
        else:
            cind=cl.nearest_centroid(x_i,state.mcmc['mean'])
            x_f=sp.random.multivariate_normal(x_i,state.alpha_Q[cind]*state.Q_cov[cind],1)[0]
            
        return x_f
        
    #Method for sampling first point (default)
    def initialize(self,state,chain):
        #If no initial distributions is given, sample randomly:
        if state.mcmc['input_data']==None:
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
        else:
            #Sample globally according to initial proposal
            x=self.sample_Q(state)
            params=dict([(key,x[state.mcmc['i_p'][key]]) for key in state.mcmc['i_p'].keys()])
        
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
        x_i=sp.array([X.params[key] for key in X.model_params])
        x_f=state.sample_Q(x_i,chain.cind,state)
        params=dict((X.model_params[i],x_f[i]) for i in range(len(X.model_params)))
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
        a=X_f.lnP-X_i.lnP-state.lnQ(X_i,X_f)+self.lnQ(X_f,X_i,state.mcmc)
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
