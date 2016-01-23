# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 19:14:14 2016

@author: Jan Lindroos
#Statistics library
"""
import scipy as sp

#weighted mean
def weighted_mean(data,weight):
    mean=sp.sum(weight*data/sp.sum(weight))
    return mean
    
def weighted_cov(data,weight,mean):
    w_sum=sp.sum(weight)
    weight=weight/w_sum
    #weight=1/float(len(data[0,:]))*sp.ones(len(data[0,:]))
    cov=0
    for j in range(len(data[0,:])):
        dmu=sp.matrix(data[:,j]-mean)
        #print dmu
        dmu_T=sp.transpose(dmu)
        #print dmu_T
        #print weight_i[j]
        cov+=weight[j]*sp.dot(dmu_T,dmu)
        #print cov

    w_sum=sp.sum(weight)
    try:
        C=1/float((1-sp.sum(weight**2)))
        cov=sp.array(C*cov)
    except:
        print 'Weighting is singular, reweighting with max point removed'
        cov=weighted_cov(data[:,weight!=weight.max()],weight[weight!=weight.max()],mean)        

    return cov

#calculate mc_error
def mc_error_rand(batch_means,batch_weights,alg='rand'):
    total_mean=weighted_mean(batch_means,batch_weights)
    if alg=='rand':
        V_1=sp.sum(batch_weights)
        V_2=sp.sum(batch_weights**2)
        sigma_2=V_1/float(V_1**2-V_2)*sp.sum(batch_weights*(batch_means-total_mean)**2)
        #print sp.sum(weights), len(weights)
        mc_error=sp.sqrt(sigma_2/float(len(batch_weights)))
    #mcmc requires same size batches
    if alg=='mcmc':
        batch_size=batch_weights#integer: All batches have an equal number of models
        n_batches=len(batch_means)#number of batches
        mc_error=batch_size/float(n_batches-1)*sp.sum((batch_means-total_mean)**2)#mcmc error according to http://arxiv.org/abs/math/0601446      
        
    return mc_error,mu_tot