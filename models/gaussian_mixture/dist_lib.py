# -*- coding: utf-8 -*-
"""
Created on Thu May 29 12:58:44 2014

@author: Jan
library for test distributions
"""
#import matplotlib
#matplotlib.use('Agg')
import scipy as sp
import sys
#Workaround for faulty python install   
try:
    from scipy.misc import logsumexp
    from scipy.linalg import inv, det
except Exception,e:
    print 'Something is wrong with your scipy:'
    print str(e)
    sys.exit() 
#from matplotlib import pyplot as plt, figure

#Generate lnP function
def lnP(dist_vars,params,debug=False):
    #Gaussian distribution
    if dist_vars['dist']=='gaussian':
        mu=dist_vars['mu']
        cov=dist_vars['cov']
        w=dist_vars['w']
        #Construct subset included in gaussian
        #print mu[0]
        #print dist_vars['vars']
        #print params
        x=sp.array([params[dist_vars['vars'][i]] for i in range(len(mu[0]))])
        lnP_k=sp.array([-0.5*(sp.dot((x-mu[i]),sp.dot(inv(cov[i]),(x-mu[i])))) for i in range(len(w))])
        #lnp=sp.log(sp.sum(w*sp.exp(lnP_k)))
        lnp=logsumexp(lnP_k,b=sp.array(w))
        
        if debug==True:
            lnp_t=0
            print 'x:',x
            for i in range(len(w)):
                lnp_t+=w[i]*sp.exp(lnP_k[i])
                print 'w_%i,lnP_%i,P_%i:'%(i,i,i),w[i],lnP_k[i],lnp_t
            lnp_t=sp.log(lnp_t)
            print 'lnP: ',lnp,lnp_t
    return lnp

#Generate random Gaussian mixture model with k gaussians within range given by span
def gen_dist(dist,scan_range,k,sigma,w=None,sym=True,uncorr=True,norm=False):    
    if dist=='gaussian':
        x_max,x_min=[[],[]]
        dist_vars={'vars':[]}
        for key in scan_range.keys():
                dist_vars['vars']+=[key]
                x_min.append(scan_range[key][0])
                x_max.append(scan_range[key][1])
    
        x_max=sp.array(x_max)
        x_min=sp.array(x_min)
        N_d=len(x_max)
        
        if w==None:
            w=sp.rand(k)
        dist_vars['dist']='gaussian'
        dist_vars['w']=w/sp.sum(w)
        dist_vars['mu']=[scl((1-3*sigma)*sp.rand(N_d)+3*sigma,x_min,x_max) for i in range(k)]
        
        SD=sp.outer((x_max-x_min),(x_max-x_min))
        if uncorr:
            if sym:
                SIGMA=sigma**2*sp.ones(N_d)*(x_max-x_min)**2
            else:
                SIGMA=sigma**2*sp.rand(N_d)*(x_max-x_min)**2
            #A=[sp.matrix(SIGMA*sp.eye(N_d)) for i  in range(k)]
            #dist_vars['cov']=[sp.multiply(A[i]*A[i].T,SD) for i  in range(k)]
            cov_N=[SIGMA*sp.eye(N_d) for i  in range(k)]
            dist_vars['cov']=[cov_N[i]*SD for i  in range(k)]
        else:
            #A=[sigma*sp.matrix(sp.rand(N_d,N_d)) for i  in range(k)]
            dist_vars['cov']=[rnd_cov(N_d,sigma)*SD for i  in range(k)]
            
    #adjust weights according to covariance
    if norm:
        dist_vars['w']=[dist_vars['w'][i]/float(sp.sqrt(det(dist_vars['cov'][i]))) for i in range(len(dist_vars['w']))]
    return dist_vars

#Method to scale unit coordinates
def scl(x,x_min,x_max):
    x=x*(x_max-x_min)+x_min
    return x

#Construct Random Covariance matrix    
def rnd_cov(N_d,sigma=0.1):
    if type(sigma)==int or type(sigma)==float:
        sigma=sigma*sp.rand(N_d)
    else:
        pass

    cov=sp.outer(sigma,sigma)
    #correlation matrix
    r=2*sp.rand(N_d)-1
    rho=sp.outer(r,r)+(1-r**2)*sp.eye(N_d)
    cov=rho*cov
    return cov

#Integrate out coordinates [marg_vars]    
def marginalize(dist_vars,marg_vars):
    #Initialize marginal dict, same for all dists
    margdist_vars={}
    margdist_vars['dist']=dist_vars['dist']
    #Gaussian
    if dist_vars['dist']=='gaussian':
        N_k=len(dist_vars['w'])#Number of gaussians
        N_D=len(dist_vars['mu'][0])#Dim of orgiginal parameter space
        
        #Initialize remaining components of marg dict, before any marginalization        
        margdist_vars['mu']=dist_vars['mu'][:]
        margdist_vars['cov']=dist_vars['cov'][:]
        margdist_vars['w']=dist_vars['w'][:]
        margdist_vars['vars']=dist_vars['vars'][:]
        
        for marg_var in marg_vars:
            #Get indices of marginalized var in current gaussian
            i_m=margdist_vars['vars'].index(marg_var)
            #Create list of current indices
            i_old=list(range(N_D))
            #remove index of marg_var
            i_old.remove(i_m)
            
            
            #remove marg_var from list of vars
            margdist_vars['vars'].remove(marg_var)
        
            margdist_vars['mu']=[sp.delete(margdist_vars['mu'][i],i_m,0) for i in range(len(margdist_vars['w']))]
            
            #For testing
#            for i in range(N_k):
#                margdist_vars['w'][i]=dist_vars['w'][i]
#                margdist_vars['cov'][i]=sp.delete(sp.delete(margdist_vars['cov'][i],i_m,0),i_m,1)
            
            #Loop over components in mixture
            #marg cov:T_M=L_m-T_m
            #marg weight:w_m=sp.sqrt(2*pi/L_mm)
            for i in range(N_k):
                #invert original covariance matrix
                Lambda=inv(sp.matrix(margdist_vars['cov'][i]))
                #Store marg compononent of 
                L_mm=Lambda[i_m,i_m]
                #Remove marginal component from Lambda
                L_m=sp.delete(sp.delete(Lambda,i_m,0),i_m,1)
                #Construct skew matrix
                l_m=sp.matrix(Lambda[i_m,i_old]+Lambda[i_old,i_m])
                T_m=l_m.T*l_m/(4*L_mm)
                #Construct marginalized covariance matrix
                margdist_vars['cov'][i]=sp.asarray(inv(L_m-T_m))
                #Scale weight
                margdist_vars['w'][i]=sp.sqrt(2*sp.pi/L_mm)*dist_vars['w'][i]
            
            #Update dimensions of marginalized parameter space
            N_D=N_D-1
         
        return margdist_vars
            
                
                