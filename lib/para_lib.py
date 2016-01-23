# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 20:15:35 2015
@author: Jan Øye lindroos

Communication library which handles the communication between
the chains
"""
import lib.io_lib as io
from lib.prepost_lib import master_init
import time
import os,sys
#MPI/multiprocessing shit
class communicate(object):
    
    def __init__(self,rank,opt,connection=None):
        if opt.mode=='multiprocessing':
            #import multiprocessing as mp
            self.rank=rank
            self.mode='multiprocessing'
            self.comm=connection
            self.n_chains=opt.chains
            
        if opt.mode=='MPI':
            from mpi4py import MPI
            if rank==0:
                self.req=[None for i in range(opt.chains-1)]
            else:
                self.req=None
            self.rank=rank
            self.mode='MPI'
            self.comm = MPI.COMM_WORLD
            self.n_chains= self.comm.Get_size()
            
     
    #Worker and Master: Send state of chain to master (master itself is a worker)
    def send_chain(self,chain):
        #MPI uses send and recv to communicate
        if self.mode=='MPI':
            print '%s: Worker %i sending chain'%(time.strftime("%Y-%m-%d %H:%M:%S"),chain.rank)
            #Make sure last update was recieved by master
            if self.req!=None:
                print '%s: Worker %i waiting for master'%(time.strftime("%Y-%m-%d %H:%M:%S"),chain.rank)
                self.req.Wait()
                print '%s: Worker %i finished waiting for master'%(time.strftime("%Y-%m-%d %H:%M:%S"),chain.rank)
            self.req=self.comm.isend(chain,dest=0,tag=0)
        #Multiprocess uses a shared queue
        if self.mode=='multiprocessing':
            #queue
            self.comm.send(chain)
        
        return 
       
    #Master methods
    def send_state(self,state):
        #print '%i: Sending state'%self.rank,state.continue_sampling           
        for i in range(self.n_chains-1):
            if self.mode=='MPI':
                #Make sure last update has been recieved by workers
                if self.req[i]!=None:
                    print "%s: Master waiting for worker %i"%(time.strftime("%Y-%m-%d %H:%M:%S"),i+1)
                    self.req[i].Wait()
                    print "%s: Master finished waiting%i"%(time.strftime("%Y-%m-%d %H:%M:%S"),i+1)
                self.req[i]=self.comm.isend(state,dest=i+1,tag=1)
                print "%s: Master sent state to %i"%(time.strftime("%Y-%m-%d %H:%M:%S"),i+1)
            if self.mode=='multiprocessing':
                if state.chains['continue_sampling'][i+1]:
                    self.comm[i].send(state)
                    #print "send state to %i"%(i+1)
                #else:
                #    print "%i completed sampling"%(i+1)
                    
        return
       
    #Master: check for and get state of workers     
    def master_check(self):
        updates=[]         
        for i in range(self.n_chains-1):
            #Test for incoming message
            if self.mode=='MPI': 
                if self.comm.Iprobe(source=i+1,tag=0):
                    #Add state to list of local states
                    updates+=[self.comm.recv(source=i+1,tag=0)]
                    print "%s: Master recieved update from %i"%(time.strftime("%Y-%m-%d %H:%M:%S"),i+1)
            if self.mode=='multiprocessing':
                if self.comm[i].poll():
                    updates+=[self.comm[i].recv()]                                
        return updates        
    
    #Worker methods
    #Worker: Check for change in state of scan    
    def worker_check(self,state,blocked=False):
        #print '%i: inside worker check:'%self.rank,state.continue_sampling
        if not blocked:
            if self.mode=='MPI':
                if self.comm.Iprobe(source=0,tag=1):
                    #state=self.comm.recv(source=0,tag=1)
                    state=self.comm.recv(source=0,tag=1)
                    print '%s: Worker %i received state'%(time.strftime("%Y-%m-%d %H:%M:%S"),self.rank)
                    
            if self.mode=='multiprocessing':
                #print '%i: Polling for state'%self.rank,state.continue_sampling
                if self.comm.poll():
                    state=self.comm.recv()
                    #print '%i: State recieved'%self.rank,state.continue_sampling
        else:
            while blocked:
                if self.mode=='MPI':
                    if self.comm.Iprobe(source=0,tag=1):
                        #state=self.comm.recv(source=0,tag=1)
                        state=self.comm.recv(source=0,tag=1)
                        blocked=False
                        print '%s: Worker %i received state'%(time.strftime("%Y-%m-%d %H:%M:%S"),self.rank)
                    
                if self.mode=='multiprocessing':
                    #print '%i: Polling for state'%self.rank,state.continue_sampling
                    if self.comm.poll():
                        state=self.comm.recv()
                        blocked=False
                        #print '%i: State recieved'%self.rank,state.continue_sampling
            
        return state
    
#Launch the scan in parallell    
def launch(scan,alg,model,opt):      
    #Set up MPI
    if opt.mode=='MPI':
        from mpi4py import MPI
        comm = MPI.COMM_WORLD 
        rank = comm.Get_rank()
    
        #let mother process write some initial info
        if rank==0:
            #print headers
            io.print_init(opt)
            #Set up directories            
            master_init(opt)            
        
        #Wait for mother before starting
        comm.barrier()
        
        scan(rank,alg,model,opt)
        
        if rank==0:
            io.print_finish(opt)
            if opt.merge_files:
                infiles=[os.path.join(opt.path,opt.run_name,'%s_%i.dat'%(opt.run_name,i)) for i in range(opt.chains)]
                outfile=os.path.join(opt.path,opt.run_name,'%s.dat'%(opt.run_name))
                io.merge_files(infiles,outfile,delete_files=True)
    
    #Otherwise use multiprocessing    
    if opt.mode=='multiprocessing':
        import multiprocessing as mp
        
        #print headers and initialize run_log
        io.print_init(opt)
        #master initialize: Sets up folders, and performs model.master_init if it exists
        master_init(opt,model)
        
        #queue for passing updates to master
        #Create Pipes for communication 
        master_ends,worker_ends=[[],[]]
        for rank in range(opt.chains-1):
            master,worker=mp.Pipe(duplex=True)
            master_ends+=[master]
            worker_ends+=[worker]
        
        processes=[]
        for rank in range(opt.chains):
            if rank==0:
                connection=master_ends
            else:
                connection=worker_ends[rank-1]    

            proc=mp.Process(target=scan, args=(rank,alg,model,opt,connection))
            proc.start()
    
            processes.append(proc)
    
        for proc in processes:
            proc.join()
            
        
        #Check if output from different files should be merged
        if opt.merge_files:
            infiles=[os.path.join(opt.path,opt.run_name,'%s_%i.dat'%(opt.run_name,i)) for i in range(opt.chains)]
            outfile=os.path.join(opt.path,opt.run_name,'%s.dat'%(opt.run_name))
            io.merge_files(infiles,outfile,delete_files=True)
            
    
    return

