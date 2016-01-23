# -*- coding: utf-8 -*-
"""
Created on Mon Jan 05 15:40:58 2015
@author: Jan Øye Lindroos
pre and post-processing:
initialization and clean up methods routines
"""
import sys,os,shutil
import scipy as sp

##turn module into
#class cls(object):
#        def __init__(self,module):
#            for cls_att in dir(module):
#                if '__' not in cls_att:
#                    #Turn lambda function into proper function
#                    if cls_att=='functions':
#                        
#                    self.__setattr__(cls_att,getattr(module,cls_att))
#
#def mod_to_obj(mod):
#
#    obj_obj=cls(mod)            
#    
#    return obj_obj
   

#takes the scan_card and sets up the run
def init_run(args):
    #Add path to scan card to system path
    [input_path,module]=os.path.split(args[1])
    sys.path.append(os.path.abspath(input_path))
    #import scan_card
    import_opt='import '+module.rstrip('.py')+' as opt'    
    exec(import_opt)
    
    #Add path to options file to opt
    opt.options_path=os.path.join(os.path.abspath(input_path),module)  
    
    #set up path if not specified in options
    if 'path' not in dir(opt):
        opt.path=os.path.abspath(input_path)
    
    [model_path,module]=os.path.split(opt.model)
    sys.path.append(os.path.abspath(model_path))
    
    #import algorithm       
    exec('import alg.'+opt.alg+' as alg')
        
    #Add default print options if not specified
    #Initialize dummy functions and constants if not present in model
    if 'model_change' not in dir(opt):
        opt.model_change={} 
    if 'constants' not in dir(opt):
        opt.constants={}   
    if 'functions' not in dir(opt):
        opt.functions={}
    if 'merge_files' not in dir(opt):
        opt.merge_files=False
    if 'print_level' not in dir(opt):
        opt.print_level=1
    if 'write_level' not in dir(opt):
        opt.write_level=1
    if 'sequential' not in dir(opt):
        opt.sequential=True
    if 'sample_size' not in dir(opt):
        if opt.alg=="nonrand":
            opt.sample_size=sp.inf
        else:
            opt.sample_size=1000
    if 'chains' not in dir(opt):
        opt.chains=1
    if 'mode' not in dir(opt):
        opt.mode='multiprocessing'
    if 'run_name' not in dir(opt):
        run_nr=0
        subdirs=[name for name in os.listdir(opt.path) if os.path.isdir(name)]
        while True:
            run_name='run_%i'%run_nr
            if run_name not in subdirs:
                opt.run_name=run_name
                break
            run_nr+=1
        
    
    #check that chains make sense compared to sample size
    if opt.chains>opt.sample_size:
        opt.chains=opt.sample_size
        
    #import model
    exec('from '+module.rstrip('.py')+' import model')
    model=init_model(opt,model)
    print 'after init_model:'
    print dir(model)
    
    if 'print_params' not in dir(opt):
        opt.print_params=model.param_names[0:5]
    
    #check print_params matches param_names after model updates
    new_print_params=[]
    for name in opt.print_params:
        if name not in model.param_names:
            continue
        else:            
            new_print_params+=[name]
    opt.print_params=new_print_params
    
    print 'returning from init_run'
    print dir(model)
    return alg,model,opt

#Modify model based on options    
def init_model(opt,model):
    #Set up data directory to be made by init_dirs
    model.datadirs={'main':os.path.join(opt.path,opt.run_name)}           
    
    #Initialize model if initialization method included
    if 'initialize' in dir(model):
        model.initialize(opt)
        print 'after initialize:'
        print dir(model)        
        
    #Add accept attribute
    model.accept=True
    #Add error attribute
    model.error=0    
    
    return model

#Sets up directories to be used    
def master_init(opt,model):    
    #Delete run directory if it already exists
    if os.path.exists(model.datadirs['main']):
        shutil.rmtree(model.datadirs['main'])
        
    #create datadirs
    os.makedirs(model.datadirs['main'])
    dirkeys=[key for key in model.datadirs.keys() if key!='main']
    for key in dirkeys:       
        os.makedirs(model.datadirs[key])
        
    #Add option file to run dir
    shutil.copy(opt.options_path,os.path.join(opt.path,opt.run_name))
    
    #Create tmp dir for each chain
    if os.path.exists('tmp'):
        shutil.rmtree('tmp')
    os.mkdir('tmp')
    #make working directories for chains
    for i in range(opt.chains):
        #print 'making directory %s'%os.path.join('tmp','tmp_%i'%i)
        os.mkdir(os.path.join('tmp','tmp_%i'%i))
        #copy model inputfiles to working directories
        if 'inputfiles' in dir(model):
            for inputfile in model.inputfiles.values():
                if os.path.isfile(inputfile):
                    shutil.copy(inputfile,os.path.join('tmp','tmp_%i'%i,os.path.split(inputfile)[1]))
                if os.path.isdir(inputfile):
                    shutil.copytree(inputfile,os.path.join('tmp','tmp_%i'%i,os.path.split(inputfile)[1]))
                    
    #Initialize model if master_initialize is included     
    if 'master_initialize' in dir(model):
        model.master_initialize(opt)
        
    return
            
            
def clean_run(opt):
    #remove tmp directory
    #if os.path.exists('tmp'):
    #    shutil.rmtree('tmp')    
    return
    

    

