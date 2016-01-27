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
from scipy.misc import logsumexp
import cPickle as pickle

   
#State of the scan
class state(alg.state):
    
    def __init__(self,rank,opt,model):
        alg.state.__init__(self,rank,opt,model)
        #Import mcmc options
        #if 'mcmc' in dir(opt):
        #    self.mcmc.update(opt.mcmc)
        self.mcmc={'input':None}
        self.mcmc['k']=opt.chains#Number of clusters
        self.mcmc['alpha_Q']=sp.ones(self.mcmc['k'])#Covariance scaling factors
        self.mcmc['beta_Q']=0.1#fraction of global sampling
        self.mcmc['n_acc']=sp.zeros(self.mcmc['k'])
        self.mcmc['n_try']=sp.zeros(self.mcmc['k'])
        #Ordinary mh covariance for non cluster sampling
        dims=sp.array([opt.scan_range[key][1]-opt.scan_range[key][0] for key in opt.scan_range.keys()])
        self.mcmc['Q_cov']=dims*(sp.eye(len(opt.scan_range.keys()))*dims).T
        #Import data for clustering
        if self.mcmc['input']==None:
            self.mcmc['clusters']=None
            pass
        else:
            #if input is a dat file
            if self.mcmc['input'].split('.')[1]=='dat':
                cluster_data=io.read_params(self.mcmc['input'],'dat',params=opt.scan_range.keys()+['weight'])
                #Get weight
                lnw=sp.log(cluster_data.pop('weight'))
                self.mcmc['clusters']=cl.cluster(cluster_data,self.mcmc['k'],lnw=lnw)
            #if input is pickle
            if self.mcmc['input'].split('.')[1]=='pickle':
                self.mcmc['clusters']=pickle.load(open(self.mcmc['input'],'rb'))    
        
    #Update states based on updates from chains
    def update(self,chain_updates):
        alg.state.update(self,chain_updates)
        for chain_update in chain_updates:
            #update global trial/error count 
            self.mcmc['n_try']+=chain_update['n_try']
            self.mcmc['n_acc']+=chain_update['n_acc']
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
            
            #Update proposal distribution
            if self.mcmc['clusters']==None:
                #Update proposal based on chain data
                clusters=cl.cluster(cluster_data,self.mcmc['k'],lnw=lnw,norm=False)
                self.mcmc['clusters']=clusters
            else:
                clusters=cl.cluster(cluster_data,self.mcmc['k'],lnw=lnw,centroids=self.mcmc['clusters']['mean'],norm=False)
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
        
    def calculate(self,X,state,X_last=None):
        X=alg.kernel.calculate(self,X,state,X_last)
        #Only calculate closest centroid if clusters exist 
        if state.mcmc['clusters']==None:
            pass
        else:            
            N_p=len(state.mcmc['clusters']['i_p'].keys())
            x=sp.zeros(N_p)
            for key in state.mcmc['i_p'].keys():
                x[state.mcmc['clusters']['i_p'][key]]=X.params[key]
            X.cind=cl.nearest_centroids(x,state.mcmc['clusters']['mean'])

        return X

    #Ammc Proposal distribution value   
    def lnQ(self,X_i,X_f,state):
        #number of parameters
        N_p=len(self.scan_range.keys())                
        #if no clusters simply normal mcmc with spherical gaussian
        if state.mcmc['clusters']==None:
            #Order not relevant due to spherical symmetry
            x_i=sp.array([X_i.params[key] for key in self.scan_range.keys()])
            x_f=sp.array([X_f.params[key] for key in self.scan_range.keys()])
            lnQ=-0.5*sp.dot(sp.dot(sp.transpose(x_i-x_f),inv(state.mcmc['alpha_Q'][0]*state.mcmc['Q_cov'])),(x_i-x_f))
        #If clusters exist include global term
        else:
            x_i,x_f=[sp.zeros(N_p),sp.zeros(N_p)]
            for key in state.mcmc['clusters']['i_p'].keys():
                x_i[state['clusters']['i_p'][key]]=X_i.params[key]
                x_f[state['clusters']['i_p'][key]]=X_f.params[key]
            
            #Add local contribution
            Q_k=[-0.5*sp.dot(sp.dot(sp.transpose(x_i-x_f),inv(state.mcmc['alpha_Q'][X_i.cind]*state.mcmc['cov'][X_i.cind])),(x_i-x_f))]
            Q_w=[(1-state.mcmc['beta_Q'])*sp.array([1/sp.sqrt(det(state.mcmc['alpha_Q'][X_i.cind]*state.mcmc['cov'][X_i.cind]))])]
            
            #Add non-local gaussian contribution
            for i in range(state.mcmc['k']):
                #exclude current cluster from large jump
                if i!=X_i.cind:
                    Q_k+=[-0.5*sp.dot(sp.dot(sp.transpose(state.mcmc['mean'][i]-x_f),inv(state.mcmc['cov'][i])),(state.mcmc['mean'][i]-x_f))]
                    Q_w+=[state.mcmc['beta_Q']*state.mcmc['weight']/float(sp.sqrt(det(state.mcmc['cov'][i])))]
    
            lnQ=logsumexp(Q_k,b=Q_w)
        return lnQ
     
    #Sample point based on Proposal distribution
    def sample_Q(self,state,X_i=None):
        #Number of parameters
        N_p=len(self.scan_range.keys())
        #if no clusters simply normal mcmc with spherical gaussian
        if state.mcmc['clusters']==None:
            #order of parameters irrelevant due to spherical symmetry
            x_i=sp.array([X_i.params[key] for key in self.scan_range.keys()])
            x_f=sp.random.multivariate_normal(x_i,state.mcmc['alpha_Q'][0]*state.mcmc['Q_cov'],1)[0]
            params_f=dict([(self.scan_range.keys()[i],x_f[i]) for i in range(N_p)])
        #Sample from clusters
        else:
            u=sp.rand()#Random number for selecting local/global jump
            #Global step
            if u<state.mcmc['beta_Q'] or X_i==None:
                #Select cluster based on cluster weights
                cind=sp.searchsorted(sp.cumsum(state.mcmc['clusters']['weight']),sp.rand())
                x_f=sp.random.multivariate_normal(state.mcmc['clusters']['mean'][cind],state.mcmc['clusters']['cov'][cind],1)[0]
            #Local step
            else:
                x_i=sp.zeros(N_p)
                #ensure order consistent with clusters
                for key in state.mcmc['clusters']['i_p'].keys():
                    x_i[state['clusters']['i_p'][key]]=X_i.params[key]
                x_f=sp.random.multivariate_normal(x_i,state.mcmc['alpha_Q'][X_i.cind]*state.mcmc['clusters']['cov'][X_i.cind],1)[0]
                
            params_f=dict([(key,x_f[state.mcmc['clusters']['i_p'][key]]) for key in self.scan_range.keys()])
            
        return params_f
        
    #Method for sampling first point (default)
    def initialize(self,state,chain):
        #If no initial distributions is given, sample randomly:
        if state.mcmc['clusters']==None:
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
        #Sample globally according to initial proposal
        else:
            print 'sampling from clusters'
            params=self.sample_Q(state)
            
        #Add constant parameters
        params.update(self.constants)
        #Add functional parameters    
        for key in self.functions.keys():
            params[key]=self.functions[key](params)
            
        modelid='%i%01i'%(self.rank,0)+'%i'%chain.accept
                
        return params,modelid
                
    #Method for proposing new point X_f given X_i
    def propose(self,X,state,chain):
        #sample new parameters based on current point and state
        params=self.sample_Q(state,X)
        #Add constants        
        params.update(self.constants)
        #Add functions
        for key in self.functions.keys():
            params[key]=self.functions[key](params)
            
        modelid='%i%01i'%(self.rank,0)+'%i'%chain.accept
        
        return params,modelid
        
    #Check wether proposal is accepted    
    def accept(self,X_i,X_f,state,u=sp.rand()):
        #Calculate kernel
        a=X_f.lnP-X_i.lnP-self.lnQ(X_i,X_f,state)+self.lnQ(X_f,X_i,state)
        if sp.log(u)<=a:
            acc=True
        else:
            acc=False
        
        return acc
    
    def weight(self,X):
        X.weight=1
        return X
    
    #Does nothing inside loop, all points are properly reweighed after
    #during postprocessing
    def reweight(self,X):
        X.weight+=1
        return X
