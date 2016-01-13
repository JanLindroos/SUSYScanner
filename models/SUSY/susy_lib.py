# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 12:37:34 2015

@author: Jan Lindroos
library for Susy tools
"""
import subprocess,os,shutil
import scipy as sp
from time import time
from slha_converter import convert_to_mg5card
import multiprocessing

run_method='proc_shell'
#Get root path to SUSYScanner
#root_path=os.path.join(os.getcwd().split('SUSYScanner')[0],'SUSYScanner')

#import paths to tools from paths.txt
libpath=os.path.dirname(os.path.abspath(__file__))
f=open(os.path.join(libpath,'paths.txt'),'r')
tool_paths=dict([(line.split()[0],line.split()[1]) for line in f.readlines()])
f.close()
    
#set up slha import dict
f=open(os.path.join(libpath,'pdgids.txt'),'r')
pdgids=dict([(line.split()[1],line.split()[0]) for line in f.readlines() if line[0]!='#'])
f.close()

#different run commands
def run_tool(cmd,mode='proc_shell',result_queue=None,tag=''):
    #cmd is a list
    if isinstance(cmd,basestring):
        cmd=cmd.split()
        
    if mode=='proc_shell':
        cmd=' '.join(cmd)
        shell=True
        stdin=subprocess.PIPE
        stdout=subprocess.PIPE
        stderr=subprocess.PIPE

        
    if mode=='mp_proc_noshell':
        shell=False
        stdin=subprocess.PIPE
        if tag!='':
            tag='_%s'%tag
        stderr=open('err%s.txt'%tag,'w+')
        stdout=open('out%s.txt'%tag,'w+')
        

    tries=0
    while tries<=3:
        try:
            proc=subprocess.Popen(cmd,shell=shell,stdin=stdin,stdout=stdout,stderr=stderr,close_fds=True)
            if mode=='mp_proc_noshell':
                [out, err] = [stdout.name,stderr.name]
                proc.wait()
                result_queue.put([out, err])
                return
            if mode=='proc_shell':
                [out, err] = proc.communicate()
            break
        except Exception as err:
            tries+=1
            print "%s did not execute properly on %i attempt..."%(cmd,tries)
            print err
            out,err=["",str(err)]
            
    if mode=='mp_proc_noshell':
        return
    
    return out,err
        

#Construct 8 digit erro with first numbers giving the tool and last specifying error
def error_code(N_tool,N_error):
    error_code=N_tool*10**(8-int(sp.log10(N_tool)))+N_error    
    return error_code

##Method for initialising slha file
def init_slha(model):
    #print os.path.join(os.getcwd(),model.slhafile)
    f_new=open(os.path.join(os.getcwd(),model.slhafile),'w')
    f_temp=open(os.path.join(libpath,'model_tmpl',model.slha_tmpl),'r')

    lines=f_temp.readlines()

    block=''
    for line in lines:
        line_list=line.split()
        if 'Block' in line:
            block=line_list[1]
        if block in ['MINPAR','EXTPAR']:
            for param_name in model.model_params:
                if param_name==line_list[3].split('(')[0]:
                    if block=="MINPAR":
                        line=line.replace(line_list[1],'%.9e'%model.params[param_name])
                    if block=='EXTPAR':
                        line=line.replace(line_list[1],'%.15e'%model.params[param_name])
        #print line
        f_new.write(line)

    f_new.close()
    f_temp.close()
    
    return
    
#Fix higgs slha for madgraph
def merge_slha(slhafile):

    H_keys=['25','35','36','37']
    H_blocks=['higgsboundsinputhiggscouplingsbosons','higgsboundsinputhiggscouplingsfermions','higgsboundsresults']
    
    #open and read SDECAY slha
    f=open(slhafile,'r')
    slines=f.read().splitlines()
    f.close()
    
    #open and read Higgsbounds slha
    f=open(slhafile+'.fh','r')
    hlines=f.read().splitlines()
    f.close()
    
    #parse relevant higgs lines
    block=''
    decay=''
    hdecays=''
    hblocks=''
    dmass=''
    dalpha=''
    extpar,mass,hmix=[{},{},{}]
    for line in hlines:
        if line.strip().split()[0]=='#':
            continue
        if line.strip().lower().split()[0]=='block':
            block=line.strip().lower().split()[1]
            decay=''
            if block in H_blocks:
                hblocks+='#\n'+line+'\n'
            continue
        if line.strip().lower().split()[0]=='decay':
            block=''
            decay=line.strip().split()[1]
            if decay in H_keys:
                hdecays+='#\n#         PDG            Width\n'+line+'\n#          BR         NDA      ID1       ID2\n'
            continue
        
        #add higgs decays
        if decay in H_keys:
                hdecays+=line+'\n'
        #add hb blocks
        if block in H_blocks:
                hblocks+=line+'\n'
                
                
        if block=='extpar' and line.strip().split()[0]=='0':
            extpar[0]=line.strip().split()[1]
        if block=='mass' and line.strip().split()[0] in H_keys:
            mass[int(line.strip().split()[0])]=line.strip().split()[1]
        if block=='dmass':
            if 'Q' in line:
                dmass+=' '
                line=line.replace(line.strip().split()[1],' '+line.strip().split()[1])
            dmass+='\t\t'+line.strip()+'\n'
        if block=='alpha':
            alpha=line.strip().split()[0]
        if block=='dalpha':
            dalpha+='\t\t   '+line.strip()+'\n'
        if block=='hmix':
            hmix[int(line.strip().split()[0])]=line.strip().split()[1]
        
            
    
    #new merged slha file
    f=open(slhafile,'w')
    block=''
    for line in slines:
        if line.strip().lower().split()[0]=='#':
            f.write(line+'\n')
            continue
        
        if line.strip().lower().split()[0]=='block':
            if block=='mass':
                f.write('BLOCK DMASS  #Higgs mass uncertainties\n'+dmass+'#\n')
            if block=='alpha':
                f.write('BLOCK DALPHA  #Higgs mixing uncertainty\n'+dalpha+'#\n')
            block=line.strip().lower().split()[1]
        
        #extract top, tau and Z mass from SM inputs
        if block=='sminputs' and line.strip().split()[0]=='6':
            M_top=line.strip().split()[1]
        if block=='sminputs' and line.strip().split()[0]=='7':
            M_tau=line.strip().split()[1]
        if block=='sminputs' and line.strip().split()[0]=='4':
            M_Z=line.strip().split()[1]
        
        if block=='extpar' and line.strip().split()[0]=='0':
            line=line.replace(line.strip().split()[1],' '+extpar[0])
            
        if block=='mass' and line.strip().split()[0] in H_keys:
            line=line.replace(line.strip().split()[1],mass[int(line.strip().split()[0])])
        #add top mass after b mass
        if block=='mass' and line.strip().split()[0]=='5':
            new_lines='         6     '+M_top+'   # t-quark pole mass\n'
            new_lines+='        15     '+M_tau+'   # tau pole mass\n'
            new_lines+='        23     '+M_Z+'   # Z pole mass'
            line=line+'\n'+new_lines
        
        if block=='hmix' and 'block' not in line.lower():
            line=line.replace(line.strip().split()[1],hmix[int(line.strip().split()[0])])
            
        if block=='alpha' and 'block' not in line.lower():
            line=line.replace(line.strip().split()[0],alpha)
            
        f.write(line+'\n')    
     
    #append decays
    #f.write('#\n# Decays from FeynHiggs version 2.10\n')
    #f.write(hdecays)
    
    #write higgsbounds blocks
    f.write('#\n# Blocks from HiggsBounds\n')
    f.write(hblocks)
    f.close()

    return 

#Method for parsing slhafiles
def parse_slha(model,blocks):
    f=open(model.slhafile,'r')
    lines=f.readlines()
    f.close()
    
    #current block
    curr_block=''
    for line in lines:
        line=line.strip()
        #skip line if comment
        if line[0]=='#':
            continue
        #change block if new block included        
        if line.split()[0].lower()=='block':
            #set current block
            curr_block=line.split()[1].lower()
            continue
        
        if curr_block in blocks:
            if curr_block=='mass':
                if 'm_%s'%pdgids[line.split()[0]] in model.param_names:
                    model.params['m_%s'%pdgids[line.split()[0]]]=abs(float(line.split()[1]))
                    
            if curr_block=='nmix':
                ij= map(int,line.split()[0:2])
                if 'N_%s'%(''.join(map(str,ij))) in model.param_names:
                    model.params['N_%s'%(''.join(map(str,ij)))]=float(line.split()[2])
            if curr_block=='gauge':
                i=int(line.split()[0])
                if 'g_%i'%i in model.param_names:
                    model.params['g_%i'%i]=float(line.split()[1])
            if curr_block=='dmass':
                if line.split()[0] in pdgids.keys():
                    if 'dm_%s'%pdgids[line.split()[0]] in model.param_names: 
                        model.params['dm_%s'%pdgids[line.split()[0]]]=abs(float(line.split()[1]))
                    
            if curr_block=='higgsboundsresults':
                if int(line.split()[0])==1:
                    if int(line.split()[1])==2:
                        hb_ok=int(line.split()[2])
                        #Change from 0,1 to 1,0 (1 meaning exclusion)
                        model.params['HB_excl']=(hb_ok+1)-2*hb_ok  
    
    return model

#method for reading in process dictionary for pythia    
def read_proccard(proccard):
    f=open(proccard,'r')
    lines=f.readlines()
    f.close()
    
    proc_dict={}
    for line in lines:
        if line[0]=='#':
            continue
        [proc_type,proc_nrs]=[line.split(':')[0],line.split(':')[1].split()]

        for proc_nr in proc_nrs:
            proc_dict[int(proc_nr)]=proc_type
            
    return proc_dict
    
#Method for creating prospino process strings
def pros_proc(state):
    procs =[]

    if state=='nn':
        subprocs=['1 1', '1 2', '1 3', '1 4', '1 5', '1 6', '1 7', '1 8',
                  '2 2', '2 3', '2 4', '2 5', '2 6', '2 7', '2 8', '3 3',
                  '3 4', '3 5', '3 6', '3 7', '3 8', '4 4', '4 5', '4 6',
                  '4 7', '4 8', '5 7', '5 8', '6 7', '6 8']
    elif state=='ll':
        subprocs=['%i 1'%(i) for i in range(1,14)]
    elif state in ['ns','ng']:
        subprocs=['1 1', '2 1', '3 1', '4 1', '5 1', '6 1', '7 1', '8 1']
    elif state in ['tb','bb']:
        subprocs=['1 1','2 1']
    else:
        subprocs=['1 1']

    for subproc in subprocs:
        procs+=['%s %s'%(state,subproc)]

    return procs
    
def gammaZ(params):
    #Leading order SUSY contributions from arXiv:hep-ph/0109283v2
    m_Z=9.11876000e+01
    s_W=0.22333#On-shell 0.23126#MSbar
    c_W=sp.sqrt(1-s_W**2)
    g_2=0.6352
    if 2*params['m_Neu1']<=m_Z:
        beta_Z=sp.sqrt((1-4*params['m_Neu1']**2/m_Z**2))
        #c_W=params['g_1']/sp.sqrt(params['g_2']**2+params['g_1']**2)
        
        g_ZXX=[]
        g_ZXX+=[g_2/(2*c_W)*(params['N_13']**2-params['N_14']**2)]
        
        gamma_Z=beta_Z**3*m_Z/(24*sp.pi)*sp.sum(g_ZXX)**2
    else:
        gamma_Z=0.0
        
    return gamma_Z


def run_softsusy(model):                          
    #command for running
    tool_nr=1    
    cmd='%s/./softpoint.x leshouches < %s'%(tool_paths['softsusy'],model.slhafile)
    out,err=run_tool(cmd,run_method)
    if len(err)>0:
        print 'err:',err
        print 'out:',out
    
    #Check output
    #Check if slha output is produced
    if 'SLHA compliant output' in out:
        #Parse slha output
        for line in out.splitlines():
            #Check for problems with point        
            if "# SOFTSUSY problem with point:" in line:
                err_string=line.split(':')[1].strip()
                ss_err=0
                for key in model.error_dict['softsusy'].keys():
                    if key in err_string:
                        ss_err+=model.error_dict['softsusy'][key]
                
                model.error=error_code(tool_nr,ss_err)
                return model
                
    else:
        if "SOFTSUSY problem" in out:
            model.error=error_code(tool_nr,model.error_dict['softsusy']["NaNs present"])
        else:
            model.error=error_code(tool_nr,model.error_dict['softsusy']["Unknown error"])
            print 'err:',err
            print 'out:',out
            print tool_nr, model.error
        return model

       
    f=open(os.path.join(model.slhafile),'w')
    f.write(out)
    f.close()
        
    if os.path.exists(model.slhafile):
        model=parse_slha(model,['mass','nmix','gauge'])
    else:
        model.error=error_code(tool_nr,model.error_dict['softsusy']["No slha"])
        print 'no slha created\n', out
        return model
        
    #check LSP
    LSP,M_LSP=0,sp.inf
    susyids=[pdgid for pdgid in pdgids.keys() if (int(pdgid)>1000000 and 'm_%s'%pdgids[pdgid] in model.params.keys())]   
    for susyid in susyids:
        M_i=model.params['m_%s'%pdgids[susyid]]
        if M_i<M_LSP:
            LSP,M_LSP=int(susyid),model.params['m_%s'%pdgids[susyid]]
    
    if LSP!=1000022:
        model.error=error_code(tool_nr,model.error_dict['softsusy']["Wrong LSP"])
        return model
    
    #Calculate Invisible Z-width and Btaunu
    model.params['Gamma_Z']=gammaZ(model.params)
    
    #Check for nans in softsusy
    #if sp.isnan([model.params[key] for key in model.part_params['softsusy']]).any():
        #print out
           
    return model    
    
def run_susyhit(model):
    tool_nr=2
    cmd=tool_paths['susyhit']+"/./run %s"%(model.slhafile)
    out,err=run_tool(cmd,run_method)
    
    if os.path.exists('%s.sh'%(model.slhafile)) and len(err)==0:
        shutil.copy('%s.sh'%(model.slhafile),'%s'%(model.slhafile))
        os.remove('%s.sh'%(model.slhafile))
    else:
        model.error=error_code(tool_nr,model.error_dict['susyhit']["Unknown error"])
        print 'susyhit out:',out
        print 'susyhit err:',err,'\n'
        
    #Check for nans in decays
    if sp.isnan([model.params[key] for key in model.part_params['susyhit']]).any():
        model.error=error_code(tool_nr,model.error_dict['susyhit']["NaN decay"])
        print out        
    
    return model

#Method for running feynhiggs   
def run_feynhiggs(model):
    tool_nr=3
    #Flags for running with stop/top resummation    
    flags="400203110"
    cmd=os.path.join(tool_paths['feynhiggs'],'build',"./FeynHiggs ")+model.slhafile+" "+flags
    out,err=run_tool(cmd,run_method)
    
    if len(err)>0:
        model.error=error_code(tool_nr,model.error_dict['feynhiggs']["Unknown error"])
        print 'FH_out:',out
        print 'FH_err:',err
        return model
        
    if not os.path.exists(model.slhafile+'.fh-001'):
        model.error=error_code(tool_nr,model.error_dict['feynhiggs']["No slha"])
        print 'no slhafile, FH_out:',out
        return model
    else:
        shutil.move(model.slhafile+'.fh-001',model.slhafile+'.fh')
    
    merge_slha(model.slhafile)
    model=parse_slha(model,['mass','dmass'])
    os.remove(model.slhafile+'.fh')      
    
    #Check for nans in higgsbounds
    if sp.isnan([model.params[key] for key in model.part_params['feynhiggs']]).any():
        print 'Feynhiggss NaNs:'
        print 'out:',out
        print 'err:',err
        model.error=error_code(tool_nr,model.error_dict['feynhiggs']["NaN output"])
        return model 
    
    return model
#    #Flags for running with stop/top resummation    
#    flags="400203110"
#    cmd=model.toolpaths['feynhiggs']+"./FeynHiggs "+model.slhafile+" "+flags
#    proc=subprocess.Popen(cmd,shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#    [out, err] = proc.communicate()
#    
#    #Merging slha files
#
#    return

#Method for running higgsbounds   
def run_higgsbounds(model):
    tool_nr=4
    cmd=tool_paths['higgsbounds']+"/example_programs/./HBSLHAinputblocksfromFH "+model.slhafile
    out,err=run_tool(cmd,run_method)
    
    for line in out.splitlines():
        if 'error' in line:
            for key in model.error_dict["higgsbounds"].keys():
                if key in line:
                    model.error=error_code(tool_nr,model.error_dict['higgsbounds'][key]) 
                    return model
            
            if model.error==0:
                model.error=error_code(tool_nr,model.error_dict['higgsbounds']["Unknown error"])
                print 'higgsbounds:',line
                return model  

    if not os.path.exists(model.slhafile+'.fh'):
        model.error=error_code(tool_nr,model.error_dict['higgsbounds']["No slha"])
        print 'no slhafile, HBSLHA_out:',out
        return model                   
               
    #run higgsbounds
    #print 'running higgsbounds' 
    cmd=tool_paths['higgsbounds']+"/./HiggsBounds LandH SLHA 3 1 "+model.slhafile+'.fh'
    #print 'finnished higgsbounds' 
    out,err=run_tool(cmd,run_method)
    if len(err)>0:
        model.error=error_code(tool_nr,model.error_dict['higgsbounds']["Unknown error"])
        print 'HB_out:',out
        print 'HB_err:',err
        return model
        
    if not os.path.exists(model.slhafile+'.fh'):
        model.error=error_code(tool_nr,model.error_dict['higgsbounds']["No slha"])
        print 'no slhafile, HB_out:',out
        return model
    
    merge_slha(model.slhafile)
    model=parse_slha(model,['mass','dmass','higgsboundsresults'])
    os.remove(model.slhafile+'.fh')      
    
    #Check for nans in higgsbounds
    if sp.isnan([model.params[key] for key in model.part_params['higgsbounds']]).any():
        print 'Higgbounds NaNs:'
        print 'out:',out
        print 'err:',err
        model.error=error_code(tool_nr,model.error_dict['higgsbounds']["NaN output"])
        return model    
    
    return model

#Method for running Micromegas   
def run_micromegas(model):
    tool_nr=5
    problems=['Dark Matter has electric charge','Dark Matter is a color particle','~o1 is not CDM']
    cmd=tool_paths['micromegas']+"/MSSM/./main_SS "+model.slhafile
    out,err=run_tool(cmd,run_method)
    
    if 'mMout' not in out:
        for problem in problems:
            if problem in out:
                model.error=error_code(tool_nr,model.error_dict['micromegas'][problem])
                return model

        model.error=error_code(tool_nr,model.error_dict['micromegas']['Unknown error'])        
        print 'out:',out
        print 'err:',err
        return model
    
    if len(err)>0:
        model.error=error_code(tool_nr,model.error_dict['micromegas']['Unknown error'])
        print 'out:',out
        print 'err:',err
        return model
    
    #Btaunu_SM found from tune to ATLAS sample
    Br_SM=1.032e-4#PDG (Same as Mahmoudi formula 0808.3144) 1.058e-4#From fit to ATLAS
    for line in out.splitlines():
        if 'mMout' in line:
            [key,value]=line.split(':')[1].strip().split('=')
            model.params[key]=float(value)
            if key=='Btaunu':
                model.params[key]=model.params[key]*Br_SM
                
            #check for NaNs in output
            if sp.isnan(model.params[key]):
                model.error=error_code(tool_nr,model.error_dict['micromegas']['NaN output'])
                print 'NaN in MM'
                print out
                return model
                
    #Remove higgs exclusion from micromegas  
    if model.params["LEP"]>=128:
        x=(model.params["LEP"]-sp.mod(model.params["LEP"],128))/128
        if x%2:
            model.params["LEP"]=model.params["LEP"]-128
      
    return model


def run_pythia(model):
    tool_nr=6
    #t_start=time()
    procs=['gg','sg','ss','sb','bb','tb','nn','ng','ns','ll','nl','btb']
    #Initialize cross sections
    for proc in procs:
        model.params['LO_%s'%proc]=0
        model.params['LO_%s_err'%proc]=0    
    
    #convert to MG process compatible slha2 input
    if model.evgen['mode']=='SUSYMG':
        out=convert_to_mg5card(model.slhafile,model.mgcard)
        model_card=model.mgcard
        pycmd='main_SUSYMG'
        keyword='mssm'
        proc_dict=read_proccard(os.path.join(libpath,'src','cards','pyprocs_SUSYMG.txt'))
    else:
        model_card=model.slhafile
        pycmd='main_SUSY'
        keyword='->'
        proc_dict=read_proccard(os.path.join(libpath,'src','cards','pyprocs_SUSY.txt'))
    
    pycard=os.path.split(model.inputfiles['pythia'])[-1]    
    nevt=model.evgen['nevt']/model.evgen['threads']
    dnevt=model.evgen['nevt']%model.evgen['threads']
    #t_init=time()
    processes=[]
    result_queue = multiprocessing.Queue()
    #set PYTHIAXML PATH    
    #print 'PYTHIA8DATA:',env['PYTHIA8DATA']

    xmlpath=os.path.join(tool_paths['pythia'],'share','Pythia8','xmldoc')

    for thread in range(model.evgen['threads']):#single thread
        nevt_i=nevt+int(thread<dnevt)
        rnd_seed=(int(time())+thread)%900000000       
        
        opt=[model_card,str(model.evgen['ecm']),str(nevt_i),model.hepmcfiles[thread],pycard,str(rnd_seed),xmlpath]
        cmd=[os.path.join(tool_paths['pythia'],'share','Pythia8','examples','./%s'%(pycmd))]+opt
        #out,err=run_tool(cmd,run_method)
        p=multiprocessing.Process(target=run_tool, args=(cmd,'mp_proc_noshell',result_queue,str(thread)))
        p.start()
        processes.append(p)
       
    for p in processes:
        p.join()
        
        #Initialise Cross sections for thread
        LO_part=dict([(proc,sp.zeros(2)) for proc in procs])
        [outfile,errfile]=result_queue.get()
        parse=True
        while parse:
            if os.path.exists(outfile):
                out=open(outfile,'r').read()
                err=open(errfile,'r').read()
                os.remove(outfile)
                os.remove(errfile)
                parse=False
            else:
                print '%s doesnt exist'%outfile
            
        
        #Parse output
        #t_parse=time()
        cs_parse=False
        for line in out.splitlines():
            if 'PYTHIA Event and Cross Section Statistics' in line:
                if 'End' in line:
                    cs_parse=False
                else:
                    cs_parse=True
                
            if keyword in line and cs_parse:
                proc_nr=int(line.split('|')[1].split()[-1])
                sigma=sp.array(map(float,line.split('|')[3].split()))
                #print proc_nr,sigma
                LO_part[proc_dict[proc_nr]]+=sigma
                #LO_part[model.evgen['proc_dict'][proc_nr]]+=sigma
            
        #Add to totals
        for proc in procs:
            model.params['LO_%s'%proc]+=LO_part[proc][0]*1e9
            model.params['LO_%s_err'%proc]+=(LO_part[proc][1]*1e9)**2
     
    model.params['LO_tot'],model.params['LO_tot_err']=[0,0]
     
    for proc in procs:
        #Take average over threads and add errors in quadrature
        model.params['LO_%s'%proc]=model.params['LO_%s'%proc]/float(model.evgen['threads'])
        model.params['LO_%s_err'%proc]=sp.sqrt(model.params['LO_%s_err'%proc])/float(model.evgen['threads'])
        #print 'LO_%s: %e+-%e'%(proc,model.params['LO_%s'%proc],model.params['LO_%s_err'%proc])
        #Add to total cross section and error
        model.params['LO_tot']+=model.params['LO_%s'%proc]
        model.params['LO_tot_err']+=model.params['LO_%s_err'%proc]**2
    
    model.params['LO_tot_err']=sp.sqrt(model.params['LO_tot_err'])
    #print '\nLO_tot: %e+-%e'%(model.params['LO_tot'],model.params['LO_tot_err'])
    
    #t_end=time()
    #print 'pythia finished in %s seconds: \n%s seconds init\n%s seconds parsing'%(t_end-t_start,t_init-t_start,t_end-t_parse)    
    
    return model
    
def run_prospino(model):
    tool_nr=7
    #initialize prospino params
    for pname in model.part_params['prospino']:
        model.params[pname]=0.
        
    procs=['gg','sg','ss','sb','bb','tb','nn','ng','ns','ll']
    pro_instates=[]
    for proc in procs:
        #Only add modes with appreciable cross section from pythia ()
        if model.params['LO_%s'%proc]>0.01*model.params['LO_tot'] or sp.isnan(model.params['LO_tot']):
            pro_instates+=pros_proc(proc)
    

    #set up pool of workers
    if len(pro_instates)>0:
        pool=multiprocessing.Pool(model.evgen['threads'])
        common_cmd=[os.path.join(tool_paths['prospino'],'./prospino_2.run'),model.slhafile]
        #print common_cmd+pro_instates
        results=[pool.apply_async(run_tool,args=(common_cmd+[state,str(model.evgen['ecm']),str(model.evgen['nlo'])],'proc_shell')) for state in pro_instates]
        for p in results:
            out,err=p.get()            
            #print out,err
        
        #parse output
        f_all=open(model.profile,'w')
        do_header=True
        for state in pro_instates:
            [proc,in_1,in_2]=[state.split()[0],int(state.split()[1]),int(state.split()[2])]
            f_proc=open('%s_%s%02i%02i_pros.dat'%(model.modelid,proc,in_1,in_2),'r')
            lines=f_proc.readlines()
            for line in lines:
                if not line.strip():
                    continue

                if proc==line.split()[0]:
                    d_line=line
                    [LO_err,NLO_err,LO,NLO]=[float(line.split()[i]) for i in [10,12,14,15]]
                    model.params['pr_LO_tot']+=LO
                    model.params['pr_LO_tot_err']+=LO_err**2
                    model.params['pr_NLO_tot']+=NLO
                    model.params['pr_NLO_tot_err']+=NLO**2
                    model.params['pr_LO_%s'%proc]+=LO
                    model.params['pr_LO_%s_err'%proc]+=LO_err**2
                    model.params['pr_NLO_%s'%proc]+=NLO
                    model.params['pr_NLO_%s_err'%proc]+=NLO_err**2
                if 'i1' in line:
                    h_line=line                    

            if do_header:
                f_all.write(h_line)
                do_header=False
            f_all.write(d_line)
            f_proc.close()
            #remove single process datafiles
            os.remove('%s_%s%02i%02i_pros.dat'%(model.modelid,proc,in_1,in_2))
            os.remove('%s_%s%02i%02i_pros2.dat'%(model.modelid,proc,in_1,in_2))
            os.remove('%s_%s%02i%02i_pros3.dat'%(model.modelid,proc,in_1,in_2))
          
        f_all.close()
        #Take root of errors
        model.params['pr_LO_tot_err']=sp.sqrt(model.params['pr_LO_tot_err'])
        model.params['pr_NLO_tot_err']=sp.sqrt(model.params['pr_NLO_tot_err'])
        for proc in procs:
            model.params['pr_LO_%s_err'%proc]=sp.sqrt(model.params['pr_LO_%s_err'%proc])
            model.params['pr_NLO_%s_err'%proc]=sp.sqrt(model.params['pr_NLO_%s_err'%proc])
                
    else:
        print 'no processes selected for NLO'               
        
    return model

def run_delphes(model):
    tool_nr=8   
    delcard=os.path.split(model.inputfiles['delphes'])[-1]
    cmd=[os.path.join(tool_paths['delphes'],"./DelphesHepMC"),delcard,model.rootfile]+model.hepmcfiles
    out,err=run_tool(cmd)
    for line in out.splitlines():
        if "ERROR" in line:
            print line
            
    if not err.strip():
        print "No results found, out:"
        model.error=error_code(tool_nr,model.error_dict['delphes']["Unknown error"])
        print out
        
    return model