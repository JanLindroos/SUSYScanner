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
    #turn on debugging
    opt.debug=False
    #Initialize random seed
    sp.random.seed(int(time.time())+rank)   
    #change to temporary working directory
    orgdir=os.getcwd()
    os.chdir(os.path.join('tmp','tmp_%i'%rank))
    
    #Setup communicator
    comm=prl.communicate(rank,opt,connection) 
    #Initialize Markov chain kernel and markov chain   
    kernel=alg.kernel(rank,opt,model)
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
    
    #Initialize scan according to alg
    if opt.debug:
        print rank,'before init'

    while True:
        #sample initial parameters
        params,modelid=kernel.initialize(state,chain)#Sample from initial distribution
        #Construct model from parameters
        X_i=model(modelid,params)

        #New****************************************
        X_i=kernel.calculate(X_i,state)
        #New****************************************
            
        #If accept update
        if X_i.accept:
            #Assign initial weight, differs for different algorithms
            X_i=kernel.weight(X_i)
            #should I update???            
            #chain.update(X_i)
            break            
            if opt.debug:
                print rank,'initialized...'
        else:
            if opt.debug:
                print rank,'initialization failed...'
            printer.print_model(X_i)
                          
        #finalize model if method exists
        if 'finalize' in dir(X_i):
            X_i.finalize()
     
    if opt.debug:
        print rank,'before main loop'               
    #Loop while global and local state permits it
    while state.continue_sampling and chain.continue_sampling:
        #sample proposal
        #print 'X_i before propose:',dir(X_i)
        params,modelid=kernel.propose(X_i,state,chain)
        X_f=model(modelid,params)        
        
        X_f=kernel.calculate(X_f,state,X_i)
        #print 'X_f after calc:',dir(X_f)
        #print 'X_i after calc:',dir(X_i)
        
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
            X_i=kernel.weight(X_i)
            #if opt.debug:
                #print rank,'model accepted', chain.accept, chain.size            
        else:
            #Write data to screen
            printer.print_model(X_f)
            #Reweigh X_i according to alg
            kernel.reweight(X_i)
            #Put X_f in chain according
            chain.update(X_f)
            if opt.debug:
                print rank,'model failed', chain.accept, chain.size
            
        #Send state if send=True
        if chain.send:
            if opt.debug:
                print rank,'sent chain update'
            #if sp.isnan(chain.updates['weight']):
            #    print rank,'has sp.nan in update, but send=', chain.send
            #    print chain.updates
            #No need to send to itself, just update state directly
            if rank==0:
                updates+=[chain.updates]
            #Send updated chain state to master
            else:
                comm.send_chain(chain.updates)
        
            if opt.debug:
                print rank,'updates sent:', chain.updates
                
        if rank==0:
            #Master: check for and get updated chain_state of workers
            while True:
                updates+=comm.master_check()
                #Master: update global state
                if opt.debug:
                    print rank,'checking for updates:',len(updates),opt.chains
                    print state.chains['continue_sampling']
                if len(updates)>=opt.chains:
                    if opt.debug:
                        print rank,'master recieved all updates:',updates
                    state.update(updates)
                    updates=[]
                
                    #Send updated state to workerss
                    comm.send_state(state)    
                    printer.print_state(state)
                    
                #If master is finished before worker then wait
                if chain.continue_sampling==False and state.continue_sampling==True:
                    if opt.debug:
                        print rank,'master waiting for workers to finish:',len(updates)
                    pass
                else:
                    break
                
        else:
            #Worker: check for status update from master
            try:
                #if opt.debug:
                    #print rank,'worker checking for state'
                state=comm.worker_check(state) 
            except:
                if opt.debug:
                    print rank,'worker finished'
                state.continue_sampling=False
        
        #Write last point and print results
        if chain.continue_sampling==False or state.continue_sampling==False:
            writer.add(X_i)
            printer.print_model(X_i)
            if rank==0:
                if opt.debug:
                    print rank,'master finishing up'
                io.print_finish(opt)
        
    #close file
    #print 'closing writer'
    if opt.debug:
        print rank,'worker exiting scan'
    writer.close()
    os.chdir(orgdir)
    
