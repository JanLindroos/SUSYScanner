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
import time

# change some text
def scan(rank,alg,model,opt,connection=None):
    #Initialize random seed
    sp.random.seed(int(time.time())+rank)   
    #change to temporary working directory
    orgdir=os.getcwd()
    os.chdir(os.path.join('tmp','tmp_%i'%rank))
    
    #Setup communicator
    comm=prl.communicate(rank,opt,connection) 
    #Initialize Markov chain kernel and markov chain   
    kernel=alg.kernel(rank,opt)
    #Initialize global state
    state=alg.state(rank,opt,model)

    if rank==0:
        #Sync states to master state
        comm.send_state(state)
        #initialize list of chain updates
        updates=[]        
    else:
        state=comm.worker_check(state,blocked=True)
   
    #Initialize chain (local state)
    chain=alg.chain(rank,opt,model,state)
    #Set up writer and printer
    writer=io.writer(rank,opt)
    printer=io.printer(rank,opt)          
    
    #check if model has sequantial target distribution
    if 'parts' in dir(model):
        partial_lnP=True
    else:
        partial_lnP=False
    
    #Initialize scan according to alg
    init=True
    while init:
        #sample initial parameters
        params,modelid=kernel.initialize(state,chain)#Sample from initial distribution
        #Construct model from parameters
        X_i=model(modelid,params)
        
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
        
        #If not simply calculate likelihood            
        else:
            X_i.calculate()
            if X_i.accept:
                X_i.accept=X_i.lnP>=opt.lnP_min
            
        #If accept update
        if X_i.accept:
            #Set init flag
            init=False
            #Assign initial weight, differs for different algorithms
            kernel.weight(X_i)
            chain.update(X_i)
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
            #Write data to file and screen        
            writer.add(X_i)
            printer.print_model(X_i)
            #Update chain
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
            #Put X_f in chain according
            chain.update(X_f)
        
        #Send state if send=True
        if chain.send:
            #No need to send to itself, just update state directly
            print 'sending data:',rank,chain.send
            if rank==0:
                updates+=[chain.updates]
            #Send updated chain state to master
            else:
                comm.send_chain(chain.updates)
        
        if rank==0:
            #Master: check for and get updated chain_state of workers
            while True:
                updates+=comm.master_check()
                #Master: update global state
                if len(updates)==opt.chains:
                    state.update(updates)
                    updates=[]
                
                    #Send updated state to workerss
                    comm.send_state(state)    
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
    
