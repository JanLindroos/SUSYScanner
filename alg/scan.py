# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 16:46:12 2014
@author: Jan Øye Lindroos

Generic scan algorithm based on the notion that all
scan algorithms included can be formulated as markov chains C={X_1,...,X_N}
detmined by an initial distribution P(X_i) for X_1 and a transition amplitude 
of going from X_i to X_i+1, Q(X_i+1|X_i,C_i), where C represents the whole state
of the chain and is only relevant for addaptive algorithms 

uniform random: Q

"""
import scipy as sp
import os
import lib.io_lib as io
import lib.para_lib as prl
#from lib.prepost_lib import init_model
#import dill
import time

# change some text
def scan(rank,alg,model,opt,connection=None):
    #Initialize random seed
    sp.random.seed(int(time.time())+rank)   
    #change directory
    orgdir=os.getcwd()
    os.chdir(os.path.join('tmp','tmp_%i'%rank))
    comm=prl.communicate(rank,opt,connection) 
    #Initialize Markov chain kernel and markov chain   
    kernel=alg.kernel(rank,opt)
    #Initialize global state
    state=alg.state(rank,opt,model)

    if rank==0:
        comm.send_state(state)         
    else:
        state=comm.worker_check(state,blocked=True)
   
    #Initialize chain (local state)
    chain=alg.chain(rank,opt,model,state)
    #Set up writer and printer
    writer=io.writer(rank,opt)
    printer=io.printer(rank,opt)
    
    if rank==0:
        updates=[]#initialize list of chain updates  
    
    #check if model has sequantial target distribution
    if 'parts' in dir(model):
        partial_lnP=True
    else:
        partial_lnP=False
    
    #Sample initial point
    init=True
    #print 'initializing'
    while init:
        params,modelid=kernel.initialize(state,chain)#Sample from initial distribution
        X_i=model(modelid,params)#Construct model from parameters
        chain.size+=1
        
        #Check if likelihood should be calculated sequentially
        if partial_lnP:                                  
            for part in X_i.parts:
                #Calculate a target component and update
                X_i.calculate(part)
 
                #check wether proposal already fails
                if opt.sequential:
                    if X_i.accept:
                        X_i.accept=X_i.lnP>=opt.lnP_min
                    if not X_i.accept:
                        break
            
            if not opt.sequential and X_i.accept:
                X_i.accept=X_i.lnP>=opt.lnP_min
                        
        else:
            #print "Calculating X_i"
            X_i.calculate()
            #print "X_i.accept:%i X_i.lnP: %f"%(X_i.accept,X_i.lnP)
            if X_i.accept:
                X_i.accept=X_i.lnP>=opt.lnP_min
                #print 'X_i accepted'
            
        #If accept update
        if X_i.accept:
            #print 'Stop initialization'
            init=False
            #chain.size+=1
            chain.accepted+=1
            #Assign initial weight, differs for different algorithms
            kernel.weight(X_i)
            if X_i.lnP>chain.lnP_max:
                chain.lnP_max=X_i.lnP
                chain.params_bf=X_i.params
        else:
            printer.print_model(X_i)
                        
            
        #finalize model if method exists
        if 'finalize' in dir(X_i):
            X_i.finalize()
                    
    #Loop while global and local state permits it
    while state.continue_sampling and chain.continue_sampling:   
        #sample proposal
        params,modelid=kernel.propose(X_i,state,chain)
        X_f=model(modelid,params)        
        
        chain.size+=1
        
        #generate random number for checking, dummyfor non-mcmc 
        u=sp.rand()        
        #Check if likelihood should be calculated sequentially
        if partial_lnP:                                  
            for part in X_f.parts:
                #Calculate a target component and update the total target
                X_f.calculate(part)                
                
                #check wether proposal already fails
                if opt.sequential:
                    if X_f.accept:
                        X_f.accept=kernel.accept(X_i,X_f,state,u)
                    if not X_f.accept:
                        break
                    
            if not opt.sequential and X_f.accept:
                X_f.accept=kernel.accept(X_i,X_f,state,u)
                        
        else:
            X_f.calculate()
            if X_f.accept:
                X_f.accept=kernel.accept(X_i,X_f,state,u)            
        
        #finalize model if method exists
        if 'finalize' in dir(X_f):
            X_f.finalize()
        
        
        
        #If accept update
        if X_f.accept:
#            chain.accepted+=1
#            chain.batch_sum+=X_i.weight*sp.array(X_i.params.values())
#            chain.batch_weight+=X_i.weight
#            if X_i.lnP>chain.lnP_max:
#                chain.lnP_max=X_i.lnP
#                chain.params_bf=X_i.params
            
            #Write data to file and screen        
            writer.add(X_i)
            printer.print_model(X_i)
            #Put X_i in chain according to alg
            chain.update(X_i)
            #Move to new point
            X_i=X_f
            #Assign initial weight, differs for different algorithms
            kernel.weight(X_i)
        else:
            #Write data to screen
            printer.print_model(X_f)
            #Reweigh X_i according to alg
            kernel.reweight(X_i)
            #Put X_f in chain according to alg
            chain.update(X_f)
        
        #Send state if send=True
        if chain.send:
            #No need to send to itself, just update state directly
            if rank==0:
                updates+=[chain]
            #Send updated chain state to master
            else:
                comm.send_chain(chain)
        
        if rank==0:
            #Master: check for and get updated chain_state of workers
            while True:
                updates+=comm.master_check()
                #Master: update global state
                if len(updates)>0:
                    state.update(updates)
                    updates=[]#new empty updates list
                
                    #Send updated state to workerss
                    comm.send_state(state)
                    
                    #only print state if all chains have reached batch size, ignore master (which is slower)
                    if all(state.chain_sizes[1:][state.chain_status[1:]==1]>=state.N_prog):
                        state.N_prog+=state.batch_size    
                        printer.print_state(state)
                    
                #If master is finished before worker then wait
                if chain.continue_sampling==False and state.continue_sampling==True:
                    pass
                else:
                    break
                
        else:
            #Worker: check for status update from master
            try:
                state=comm.worker_check(state) 
            except:
                state.continue_sampling=False
        
        #Write last point and print results
        if chain.continue_sampling==False or state.continue_sampling==False:
            writer.add(X_i)
            printer.print_model(X_i)
            if rank==0:
                io.print_finish(opt)
        
    #close file
    #print 'closing writer'
    writer.close()
    os.chdir(orgdir)
    
