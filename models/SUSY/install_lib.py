# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 12:42:36 2015
@author: Jan Lindroos

SUSY install lib
"""
import os, shutil, sys
import subprocess
import tarfile
import types
from multiprocessing import cpu_count


#Method for extracting tools    
def extract(tool,source_file):
    #Create out folder
    [out,err]=['',0]
    source_path=os.path.join(os.getcwd(),'src',source_file)
    tool_path=os.path.join(os.getcwd(),'tools',tool)
    #extract source    
    try:
        tar_file=tarfile.open(source_path,'r:gz')
        tar_file.extractall(tool_path)
        tar_file.close()
    except:
        out='Unable to extract from %s'%(source_path)
        err=1     
        
    return out,err
    
#Method for running bash commands    
def run_cmd(cmd):
    try:
        proc=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
        [out, err] = proc.communicate()
    except:
        [out,err]=['','unable to run cmd']
    return out,err

#Method for installing SUSY tools
def install_tool(tool,options):
    output,status=['',0]
    if tool=='softsusy':
        output,status=install_softsusy(options)  
    if tool=='susyhit':
        output,status=install_susyhit(options) 
    if tool=='feynhiggs':
        output,status=install_feynhiggs(options) 
    if tool=='higgsbounds':
        output,status=install_higgsbounds(options) 
    if tool=='micromegas':
        output,status=install_micromegas(options) 
    if tool=='pythia':
        output,status=install_pythia(options)
    if tool=='prospino':
        output,status=install_prospino(options)
    if tool=='delphes':        
        output,status=install_delphes(options)
    #Auxiliiary tools
    if tool=='hepmc':
        output,status=install_hepmc(options)
    if tool=='root':
        output,status=install_root(options)
    
    return output,status

#Method for testing the tools    
def test_results(test_params,bench_params):
    results=''
    
    #print model data (10 parameters per line)
    N_perline=10
    N_params=len(test_params.keys())
    N_lines=N_params/N_perline+1
    pnames=sorted(test_params.keys())
    print pnames
    for i in range(N_lines):
        j_range=range(N_perline*i,min(N_perline*i+9,N_params))
        results+=''.join(['\n%15s'%('')]+['%15s'%(pnames[j]) for j in j_range])
        results+=''.join(['\n%15s'%('benchmark')]+['%15s'%('%.3e'%(bench_params[pnames[j]])) for j in j_range])
        results+=''.join(['\n%15s'%('test')]+['%15s'%('%.3e'%(test_params[pnames[j]])) for j in j_range])
        #results+=''.join(['\n%15s'%('benchmark')]+['%15s'%(str(round(bench_params[test_params.keys()[j]],3))) for j in j_range])
        #results+=''.join(['\n%15s'%('test')]+['%15s'%(str(round(test_params[test_params.keys()[j]],3))) for j in j_range])
        results+='\n'
        
    return results

#*********************Tool installation functions******************************      
def install_softsusy(options):
    #install softsusy
    output,status=['',1]
    
    #configure
    output+='\nconfiguring....\n'
    config_flags=['--enable-fast-install']
    cmd=' '.join(['./configure']+config_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
        return output,status
    #make
    output+='\nmaking....\n'
    make_flags=['programs','-j']
    cmd=' '.join(['make']+make_flags)   
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
    
    return output,status

#Tool installation functions       
def install_susyhit(options):
    #install susyhit
    output,status=['',1]
    output+='\nmaking....\n'
    make_flags=['-j']
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
    
    return output,status
    
#Tool installation functions       
def install_feynhiggs(options):
    output,status=['OK',1]
    #configure
    output+='\nconfiguring....\n'
    config_flags=['--prefix=install']
    cmd=' '.join(['./configure']+config_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
        return output,status
        
    #make
    output+='\nmaking....\n'
    make_flags=['-j']
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
        return output,status
        
    #make install
    output+='\nmake install....\n'
    cmd='make install'
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0    
    
    return output,status
    
#Tool installation functions       
def install_higgsbounds(options):
    output,status=['',1]
    #Set compiler flags and location of feynhiggs
    flags={'F90C':options['F90C'],'F77C':options['F77C'],
    'F90FLAGS':'-fcheck=bounds -ffixed-line-length-none',
    'FHINCLUDE':'-I'+options['FHINCLUDE'],
    'FHLIBS':'-L'+options['FHLIBS']+' -lFH'}
    output+='set compilers and paths for feynhiggs:\n%s'%('\n'.join(['%s=%s'%(key,flags[key]) for key in flags.keys()]))
    
    output+='\nchanging configure file...\n'
    f=open('configure','r')
    lines=f.readlines()
    f.close()
    
    f=open('configure','w')
    for line in lines:
        for key in flags.keys():
            if line[0]!='#' and key in line:
                line=key+'='+flags[key]
        f.write(line.strip()+'\n')     
    f.close()
    
    #configure
    output+='\nconfiguring....'
    config_flags=[]
    cmd=' '.join(['./configure']+config_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
        return output,status
    
    #make
    output+='\nmaking....\n'    
    make_flags=[]
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0    
    
    return output,status
    
#Tool installation functions       
def install_micromegas(options):
    output,status=['',1]
    
    #make main micromegas
    output+='\nmaking....\n' 
    make_flags=[]
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
        return output,status     
    
    #make main_SS
    output+='\nmaking....\n'
    os.chdir(os.path.join(os.getcwd(),'MSSM'))
    make_flags=['main=main_SS.c']
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0   
    
    return output,status
    
#Tool installation functions (Requisites: HEPMC)      
def install_pythia(options):
    output,status=['',1]
    main_path=os.getcwd()
    #configure    
    output+='\nconfiguring....'
    config_flags=['--with-hepmc2=%s'%(options['HEPMC'])]
    cmd=' '.join(['./configure']+config_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
        os.chdir(main_path)
        return output,status
    
    #make pythia 
    output+='\nmaking....\n'
    make_flags=[]
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
        os.chdir(main_path)
        return output,status
    
    #make pMSSM processes
    output+='\nmaking MG processes.\n'
    os.chdir(os.path.join(main_path,'Processes_mssm'))
    make_flags=[]
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
        os.chdir(main_path)
        return output,status

    #make SUSY programs
    output+='\nmaking main_SUSY.\n'
    os.chdir(os.path.join(main_path,'share','Pythia8','examples'))
    #Program using built-in SUSY processes
    make_flags=['main_SUSY']
    cmd=' '.join(['make']+make_flags)
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    #Program using Madgraph SUSY processes
    output+='\nmaking main_SUSYMG.\n'
    make_flags=['main_SUSYMG']
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
     
    if err=='unable to run cmd':
        status=0
        os.chdir(main_path)
        return output,status

    #Set pythia environment variable
    os.chdir(main_path)   
    os.system('export PYTHIA8DATA=%s'%(os.path.join(main_path,'share','Pythia8','xmldoc')))    
    
    return output,status
    
#Tool installation functions       
def install_delphes(options):
    #install softsusy
    output,status=['',1]
#    base_path=os.getcwd().split('/delphes')[0]
#    os.system('export PYTHIA8=%s'%(os.path.join(base_path,'pythia','pythia_8.2.01')))
#    os.system('export PYTHIA8DATA=%s'%(os.path.join(base_path,'pythia','pythia_8.2.01','share','Pythia8','xmldoc')))
    #configure
    output+='\nconfiguring....\n'
    config_flags=['']
    cmd=' '.join(['./configure']+config_flags)
    out,err=run_cmd(cmd)#install command
    print out,err
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
        return output,status
    #make
    output+='\nmaking....\n'
    make_flags=['-j']
    cmd=' '.join(['make']+make_flags)   
    out,err=run_cmd(cmd)#install command
    print out,err
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err
    if err=='unable to run cmd':
        status=0
    
    return output,status
    
def install_prospino(options):
    output,status=['',1]
    #get path
    pro_path=os.path.join(os.getcwd(),"")    
    #Set path and compiler in Makefile
    output+='Changing Makefile'
    f=open('Makefile','r')
    makelines=f.readlines()
    f.close()
    os.remove('Makefile')
    f=open('Makefile','w')
    for line in makelines:
        print line,line.split('=')[0].strip()
        if line.split('=')[0].strip()=='DIRECT':
            line='DIRECT=%s\n'%(pro_path)
        f.write(line)
    f.close()
    
    #make main micromegas
    output+='\nmaking....\n'
    make_flags=[]
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    output+='*****************out**********:\n%s'%out
    output+='\n****************err*********:\n%s'%err     
    
    return output,status
#******************************************************************************
    
    
#Auxilliary tools**************************************************************        
def check_root():
    output,status=['',1]
    #try running root from command line
    output+='checking ROOT cmd...:\n'
    try:
        proc=subprocess.Popen(['root','-l','-b','-q'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        out,err=proc.communicate()
        if err==None and len(out.strip())==0:
            output+='Root testrun successful\n'
        else:
            output+='Unable to run ROOT cmd\n'
            status=0
            return output,status
    except:
        output+='Unable to run ROOT cmd\n'
        status=0
        return output,status
        
    #check pyroot
    output+='checking pyroot...:\n'
    try:
        import ROOT
        output+='ROOT module successfully loaded'
    except:
        output+='Unable to load ROOT module\n'
        status=0
        return output,status
        
    #checking environmen variable
    try:
        output+='Root detected in %s\n'%(os.environ['ROOTSYS'])
    except:
        output+='Root environment variables not set, please source thisroot.sh\n'
        status=0
     
    return output,status

def install_root(options):
    output,status=['',1]
    #install    
    print 'Install ROOT does nothing at the moment'
    #test    
    
    return output,status
      
def install_hepmc(options):
    output,status=['',1]
    #create install and build directories (outside main directory recommended)
    main_path=os.getcwd()
    
    #move to source dir
    os.chdir(os.path.join(main_path,'src'))    
    #run ./bootsrap to do autconf/make    
    cmd='./bootstrap'
    out,err=run_cmd(cmd)
    
    #move to build directory for configuration
    os.chdir(os.path.join(main_path,'build'))   
    #configure
    config_flags=['--prefix=%s'%(os.path.join(main_path,'install')),'--with-momentum=GEV','--with-length=MM']
    cmd=' '.join([os.path.join(main_path,'src','./configure')]+config_flags)
    out,err=run_cmd(cmd)
    print '*****************out**********:\n',out
    print '\n****************err*********:\n',err
    #make
    make_flags=['']
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    print '*****************out**********:\n',out
    print '\n****************err*********:\n',err
    #make install
    make_flags=['install']
    cmd=' '.join(['make']+make_flags)
    out,err=run_cmd(cmd)#install command
    print '*****************out**********:\n',out
    print '\n****************err*********:\n',err
        
    return output,status    
    

