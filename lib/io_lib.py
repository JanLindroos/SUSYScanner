
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 04 16:42:15 2015
@author: Jan Øye Lindroos

IO library methods for reading and writing to file
and printing to screen
"""
import os,sys
import scipy as sp
sp.seterr(divide='ignore')#Inore errors related to log giving inf
import pickle
import time

#Writer class for writing models to file
class writer(object):    
    
    def __init__(self,rank,opt):       
        self.rank=rank
        self.path=opt.path
        self.format=opt.fileformat
        if self.format=='hdf5':
            import h5py
            self.filename=opt.run_name+'_'+str(rank)+'.hdf5'
            
        else:
            self.filename=opt.run_name+'_'+str(rank)+'.dat'
            self.f=open(os.path.join(opt.path,opt.run_name,self.filename),'w',0)
            self.header=True#Determines weather the header shoul be written to file
            
        #self.f=h5py.File(os.path.join(self.path,self.filename), 'w')
    
    #write model to file    
    def add(self,model):
        if self.format=='dat':
            
            #Write header data
            if self.header:
                h_str="".join(['%15s'%key for key in ['modelid','error','lnP','weight']+model.param_names])
                self.f.write(h_str+'\n')                
                self.header=False
            d_str=''.join(['%15s'%model.modelid]+['%15.6e'%value for value in [model.error,model.lnP,model.weight]+[model.params[key] for key in model.param_names]])
            self.f.write(d_str+'\n') 
            
        return
    
    #close file
    def close(self):
        self.f.close()
        return

class printer(object):
    
    def __init__(self,rank,opt):
        self.rank=rank
        self.print_params=opt.print_params
        self.print_level=opt.print_level
             
        
    def print_model(self,X):
        if self.print_level>1 and X.accept:
            #print time.strftime("%Y-%m-%d %H:%M:%S")
            d_str=''.join(['%7i'%X.accept,'%7i'%self.rank,'%10i'%X.error,'%15s'%X.modelid]+['%15.6e'%value for value in [X.weight,X.lnP]+[X.params[name] for name in self.print_params]])
            sys.stdout.write(d_str + '\n')    
        if self.print_level>2 and not X.accept:
            d_str=''.join(['%7i'%X.accept,'%7i'%self.rank,'%10i'%X.error,'%15s'%X.modelid]+['%15.6e'%value for value in [X.weight,X.lnP]+[X.params[name] for name in self.print_params]])
            sys.stdout.write(d_str + '\n')
                        
            
        return
        
    def print_state(self,state):
        #Set up parameters to print
        if self.print_level>0:
            i_max=sp.argmax(state.lnP_max)
            print '\n%s: %i of %i models accepted (%4.1f %%), global best fit point:'%(time.strftime("%Y-%m-%d %H:%M:%S"),state.accepted,state.stop_size,100*state.accepted/float(state.sample_size))
            
            #Print global best fit and means
            
            h_str=''.join(['%15s'%key for key in ['','chain','lnP_max']+self.print_params])
            d_str=''.join(['%15i'%i_max]+['%15.6e'%value for value in [state.lnP_max[i_max]]+[state.params_bf[i_max][name] for name in self.print_params]])
            print '%s\n%15s%s'%(h_str,'best fit:',d_str)
        
            imap=[state.params_bf[i_max].keys().index(name) for name in self.print_params]
        
            #print 'print:',state.mean
            mu_str=''.join(['%15s'%'***','%15s'%'***']+['%15.6e'%state.mean[i] for i in imap])
            print '%15s%s'%('mean:',mu_str)
            std_str=''.join(['%15s'%'***','%15s'%'***']+['%15.6e'%state.mean_error[i] for i in imap])
            print '%15s%s'%('mean_error:',std_str)
            
#            if len(state.chain_means)>1:       
#                print 'chains:'
#                for i in range(len(state.chain_means)):
#                    mu_str=''.join(['%15s'%'***','%15s'%'***']+['%15.6e'%state.chain_means[i][j] for j in imap])
#                    print '%15s%s'%('%i mean:'%(i),mu_str)
#                    std_str=''.join(['%15s'%'***','%15s'%'***']+['%15.6e'%state.chain_mean_errors[i][j] for j in imap])
#                    print '%15s%s'%('%i mean_error:'%(i),std_str)
        if self.print_level>1:
            #print header
            h_str=''.join(['%15s'%key for key in ['chain','modelid','weight','lnP']+self.print_params])
            sys.stdout.write('\n'+h_str+'\n')
            
        
        return    

def read_params(filename,filetype,params=None,skip_header=0,skip_footer=0):    
    if filetype=='dat':
        #print "Reading data from %s..."%(filename)
        #skip_header : int, optional
        #The number of lines to skip at the beginning of the file.

        #skip_footer : int, optional
        #The number of lines to skip at the end of the file.
        p=sp.genfromtxt(filename,names=True)
                
        #transform into normal dict
        p_dict={}        
        for key in p.dtype.names:
            if p[key].size==1:
                x=p[key].tolist()
                x=sp.array([x])
            else:
                x=p[key].copy()
                
            if key.lower()=='modelid':
                p_dict['modelid']=map(int,x)
                continue
                    
            if params==None:
                p_dict[key]=x
            else:
                if key in params:
                    p_dict[key]=x    
                
        p=p_dict
     
    if filetype=='pickle':
        p=pickle.load(open(filename, 'rb' ))
        #check for modelid:
        for key in p.keys():
            if key.lower=='modelid':
                p['modelid']=p[key]
        
        if params==None:
            pass
        else:
            p=dict([(name,p[name]) for name in p.keys() if name in params])
     
    if filetype=='h5py':
        import h5py

    if filetype=='slha':
        f.open(filename,'r')
        f.close()        
        
    return p    
    

#Print initial stuff about model    
def print_init(opt):
    print ' %s '%(os.path.split(os.getcwd())[1]).center(40,'*')
    print '\nModel: %s\nAlgorithm: %s\nMode: %s'%(opt.model,opt.alg,opt.mode)
    if isinstance(opt.scan_range,basestring):
        print "models imported from from %s\n"%(opt.scan_range)
    if isinstance(opt.scan_range,dict):
        print "Scan range:"
        for key in opt.scan_range.keys():
            print "%s: "%key,opt.scan_range[key]
        print "\n"
    
    if opt.print_level>1:
        #print header
        h_str=''.join(['%7s%7s%10s'%('accept','chain','error')]+['%15s'%key for key in ['modelid','weight','lnP']+opt.print_params])
        sys.stdout.write(h_str + '\n')
    return
    
def merge_files(infiles,outfile,delete_files=False):
    f_out=open(outfile,'w')
    write_header=True
    for infile in infiles:
        f_in=open(infile,'r')
        header_parsed=False
        lines=f_in.readlines()
        f_in.close()
        if delete_files:
            os.remove(infile)
        for line in lines:
            if not header_parsed:
                header_parsed=True
                if write_header:
                    f_out.write(line)
                    write_header=False
                else:
                    continue
            
            else:
                f_out.write(line)
    
    f_out.close()    
    return
    
def print_finish(state):
    print '\nscan finished in %s seconds'%('(to appear)')
    print 'stats:\n (to appear)'    
    return

#Write info about run to file    
def write_opt(opt):
    return

#Write some basic run statistics    
def write_stats(state):
    return
    