# -*- coding: utf-8 -*-
"""
Created on Tue May 13 20:45:25 2014

@author: Jan

Library for clustering, uses cython for speedup
"""

import scipy as sp
from stat_lib import weighted_cov
from math import *

##Method for clustering data
#Input: 
#data: dictionary data['params']=numpy.array
#k: number of clusters
#maxit: maximum number of iterations
#Output: list of normalized cluster indices, centroids, cluster weights and covariances 
def cluster(data_dict,k,norm=False,maxit=100,lnw=None,dcmin=1e-5,do_pp=1,centroids=None,return_data=False):
    #print 'Setting up cluster data...'
    clusters={}
    #set up data in array form
    N_d,N_p=[len(data_dict.values()[0]),len(data_dict.keys())]
    
    #data weights, based on likelihood
    if lnw==None:
        weight=sp.ones(N_d)/float(N_d)
    else:
        lnw_max=lnw.max()
        weight=sp.exp(lnw-lnw_max-sp.log(sp.sum(sp.exp(lnw-lnw_max))))    

    data=sp.zeros((N_p,N_d))
    #data_norm=data[0:l_d,:].copy()
    #clusters['vars']={}
    if norm:
        clusters['norm']={}
    clusters['i_p']={}
    for i in range(len(data_dict.keys())):
        #Normalize data
        key=data_dict.keys()[i]
        if norm:
            data_i=data_dict[key]
            clusters['norm'][key]={'max':data_i.max(),'min':data_i.min()}
            data[i,:]=(data_i-clusters['norm'][key]['min'])/float(clusters['norm'][key]['max']-clusters['norm'][key]['min'])
        else:
            data[i,:]=data_dict[key]
        clusters['i_p'][key]=i
          

    passed=0
    while passed==0:
        #initialize centroids
        if centroids==None:
            print 'Initializing centroids...'
            centroids=init_centroids(data,weight,k,do_pp)
            print centroids

        for i in range(maxit):
            #Find nearest centroid
            #print '%i. iter: Partitioning data...'%i
            cind=nearest_centroids(data,centroids)
            #print '%i. iter: New centroids...'%i
            c_means,dc,w_s=new_centroids(data,weight,cind,centroids)
            #print c_means
            dc_max=dc.max()
            #print 'max change:%f'%dc_max

            if dc_max<=dcmin:
                break;
            else:
                centroids=c_means

        #Normalize weights
        clusters['weight']=sp.array(w_s)/float(sum(w_s))
        #print 'weight:\n',clusters['weight']
        clusters['mean']=[centroids[:,i] for i in range(len(centroids[0,:]))]
        #print 'mean:\n',clusters['mean']
        clusters['cov']=[weighted_cov(data[:,cind==c_i],weight[cind==c_i],clusters['mean'][c_i]) for c_i in range(k)]
        clusters['size']=N_d

        #for i in range(len(clusters['weight'])):
        #    print '\ncluster %i: weight=%f'%(i,clusters['weight'][i])
        #    print 'cov:\n',clusters['cov'][i]
        #    print 'mean:\n',clusters['mean'][i]
            
        #print 'data points:',N_d

        passed=1
    
    if return_data:    
        cluster_data={'cind':cind,'data':data,'weight':weight}
        return clusters,cluster_data
        
    return clusters

## Method for initializing centroids    
def init_centroids(x,weight,k,do_pp):
    y=sp.cumsum(weight)/float(sp.sum(weight))
    for i in range(k):
        r=sp.rand()
        #Find random initial centroid from weigthed sample
        if i==0:
            i_c=sp.searchsorted(y,r)
            x_s=x[:,i_c].copy().reshape(-1,1)
                
        else:
            if do_pp:
                if i>1:
                    cind=nearest_centroids(x,x_s)
                    #Calculate distance to centroid from all points
                    dc_2=sp.sum((x-x_s[:,cind])**2,0)
                else:
                    cind=sp.zeros(len(x[0,:]))
                    dc_2=sp.sum((x-x_s)**2,0)

                y=sp.cumsum(weight*dc_2)/float(sp.sum(weight*dc_2))

            
            i_c=sp.searchsorted(y,r)
            x_snew=x[:,i_c].copy().reshape(-1,1)

            x_s=sp.hstack((x_s,x_snew))

    return x_s
    
def nearest_centroids(x,c):    
    #only one point
    if len(sp.shape(x))==1:
        c=sp.array(c)
        dc=sp.sum((c-x)**2,1)
        cind=sp.argmin(dc)
        return cind

    
    dc=sp.zeros((len(c[0,:]),len(x[0,:])))
    for i in xrange(len(c[0,:])):
        dc[i,:]=sp.sum(sp.subtract(x.T,c[:,i])**2,1)

    cind=sp.argmin(dc,0)
    return cind
    
def new_centroids(x,weight,cind,c):
    w_s=[]
    w_s=sp.array([])
    x_s=c.copy()
    cl_ind=sp.unique(cind)
    dc=sp.zeros(len(cl_ind))
    for i in range(len(cl_ind)):
        #w_s+=[sp.sum(weight[cind==i])]
        weight_sum=sp.sum(weight[cind==i])
        w_s=sp.append(w_s,weight_sum)
        P=weight[cind==i]/float(weight_sum)
        for j in range(len(x[:,0])):
            x_s[j,i]=sp.sum(P*x[j,cind==i])
            dc[i]+=(x_s[j,i]-c[j,i])**2
        dc[i]=sp.sqrt(dc[i])

    return x_s,dc,w_s
