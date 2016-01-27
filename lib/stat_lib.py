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

#combine means of two samples weighted by sample size    
def combine_means(means,sizes):
    new_mean=(means[0]*sizes[0]+means[1]*sizes[1])/float(sizes[0]+sizes[1])
    return new_mean

#combines two covariances     
def combine_covs(covs,means,sizes):
    N_d=len(means[0])
    size_tot=sp.sum(sizes)
    new_cov=sp.zeros((N_d,N_d))
    for i in range(N_d):
        for j in range(N_d):
            new_cov[i,j]=(sizes[0]*covs[0][i,j]+sizes[1]*covs[1][i,j])/float(size_tot)+(means[0][i]-means[1][i])*(means[0][j]-means[1][j])*(sizes[0]*sizes[1])/float(size_tot**2)
    return new_cov

#standard mc_error for weighted batches (tested for uniform sampling on gaussian distribution)
def mc_error(total_mean,batch_means,batch_weights):
    #batch_weights>1 to make sense
    if len(batch_weights)<2:
        return sp.nan
    
    V_1=sp.sum(batch_weights)
    V_2=sp.sum(batch_weights**2)
    sigma_2=V_1/float(V_1**2-V_2)*sp.sum(batch_weights*(batch_means-total_mean)**2)
    mc_error=sp.sqrt(sigma_2/float(len(batch_weights)))         
    return mc_error