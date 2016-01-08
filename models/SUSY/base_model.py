# -*- coding: utf-8 -*-
"""
Created on Sat Dec 27 15:51:37 2014

@author: Jan
This file contains the SUSY model class
"""
import scipy as sp
sp.seterr(divide='ignore')#Inore errors related to log giving inf
import susy_lib as sl
from collections import OrderedDict
import os,shutil,sys


#Model class
#SUSY base class
class model(object):   
    
    #parameter names for Observables (These are the same in all SUSY models)
    #For discrete constraints 1 means ok, 0 means failed
    #Theory should be classified (ESWB, Tachyons, etc)
    part_params={'softsusy':['m_W','m_h0','m_A0', 'm_H0', 'm_Hp','g_1', 'g_2', 'Gamma_Z',
                'm_Neu1', 'm_Neu2', 'm_Neu3', 'm_Neu4', 'm_Ch1', 'm_Ch2', 
                'N_13', 'N_14', 'm_gluino', 'm_st_1', 'm_st_2', 'm_sb_1', 'm_sb_2',  
                'm_ss_L', 'm_ss_R', 'm_sc_L', 'm_sc_R', 'm_su_L', 'm_su_R', 'm_sd_L', 'm_sd_R', 
                'm_stau_1', 'm_stau_2', 'm_smu_L', 'm_smu_R', 'm_se_L', 'm_se_R', 
                'm_snu_t', 'm_snu_mL', 'm_snu_eL'], 
                'susyhit':[],
                'feynhiggs':['m_h0','dm_h0','m_A0','dm_A0','m_H0','dm_H0','m_Hp','dm_Hp'],
                'higgsbounds':['HB_excl','m_h0','dm_h0','m_A0','dm_A0','m_H0','dm_H0','m_Hp','dm_Hp'],
                'micromegas':['LEP','Oh2', 'Oh2_bino', 'Oh2_hsino1', 
                'Oh2_hsino2', 'Oh2_wino','a_mu','drho', 'bsg','Bsmumu','Btaunu','sigma_nSD', 'sigma_nSI', 'sigma_pSD', 'sigma_pSI'],
                'pythia':['LO_tot', 'LO_tot_err', 'LO_gg', 'LO_gg_err', 'LO_sg', 'LO_sg_err',
                'LO_ss', 'LO_ss_err', 'LO_sb', 'LO_sb_err', 'LO_bb', 'LO_bb_err','LO_tb', 'LO_tb_err',
                'LO_nn', 'LO_nn_err', 'LO_ng', 'LO_ng_err', 'LO_ns', 'LO_ns_err',
                'LO_ll', 'LO_ll_err','LO_nl', 'LO_nl_err','LO_btb', 'LO_btb_err'],
                'delphes':[],
                'prospino':['pr_LO_tot', 'pr_NLO_tot', 'pr_LO_tot_err', 'pr_NLO_tot_err', 'pr_LO_gg', 'pr_NLO_gg', 'pr_LO_gg_err',
                            'pr_NLO_gg_err', 'pr_LO_sg', 'pr_NLO_sg', 'pr_LO_sg_err', 'pr_NLO_sg_err', 'pr_LO_ss', 'pr_NLO_ss',
                            'pr_LO_ss_err', 'pr_NLO_ss_err', 'pr_LO_sb', 'pr_NLO_sb', 'pr_LO_sb_err', 'pr_NLO_sb_err', 'pr_LO_bb',
                            'pr_NLO_bb', 'pr_LO_bb_err', 'pr_NLO_bb_err', 'pr_LO_tb', 'pr_NLO_tb', 'pr_LO_tb_err', 'pr_NLO_tb_err',
                            'pr_LO_nn', 'pr_NLO_nn', 'pr_LO_nn_err', 'pr_NLO_nn_err', 'pr_LO_ng', 'pr_NLO_ng', 'pr_LO_ng_err',
                            'pr_NLO_ng_err', 'pr_LO_ns', 'pr_NLO_ns', 'pr_LO_ns_err', 'pr_NLO_ns_err', 'pr_LO_ll', 'pr_NLO_ll',
                            'pr_LO_ll_err', 'pr_NLO_ll_err', 'pr_LO_nl', 'pr_NLO_nl', 'pr_LO_nl_err', 'pr_NLO_nl_err', 'pr_LO_btb',
                            'pr_NLO_btb', 'pr_LO_btb_err', 'pr_NLO_btb_err']}
                
    error_dict={'softsusy':{"Higgs potential ufb": 2,"tachyon":4,"m3sq-problem":8,"Inaccurate Higgs mass":16,"No convergence":32,"Numerical problemThrown":64,"Wrong LSP":128, "No slha":256,"NaNs present":512,"Unknown error":1024},
                'susyhit':{'NaN decay': 2,'Unknown error':4},
                'feynhiggs':{'negative sbottom mass squares':2,'DRbartoOS transition failed':4,'error computing Mh1':8,'No slha':16,'NaN output':32,'Unknown error':64},
                'higgsbounds':{'negative sbottom mass squares':2,'DRbartoOS transition failed':4,'error computing Mh1':8,'No slha':16,'NaN output':32,'Unknown error':64},
                'micromegas':{'Dark Matter has electric charge':2,'Dark Matter is a color particle':4,'~o1 is not CDM':8,'NaN output':16,'Unknown error':32},
                'pythia':{'Unknown error':2},
                'delphes':{'Unknown error':2},
                'prospino':{'Unknown error':2}}
    
    #The different likelihood components for SUSY
    parts=['softsusy','susyhit','feynhiggs','higgsbounds','micromegas','pythia','delphes','prospino']          
    #import constraints
    
    #Create the list of observables included
    obs_params=list(OrderedDict.fromkeys(sum([part_params[key] for key in parts],[])))
                                                          
    #Constraints: construct components of lnP****************************************************************************                  
    constraints=dict([(part,[]) for part in parts])
    #constraints['higgsbounds']+=[lambda pars: sp.log(int(pars['HB_excl']==0))]
    #mh0_mu,mh0_sigma=[125.5 , 1.7]
    #constraints['higgsbounds']+=[lambda pars: -(pars['m_h0']-mh0_mu)**2/(2*mh0_sigma**2)]
    #constraints['micromegas']+=[lambda pars: sp.log(int(pars['LEP']==0))]
    #Oh2_mu,Oh2_sigma=[0.1199 , 0.0027]
    #constraints['micromegas']+=[lambda pars: -int(pars['Oh2']>Oh2_mu)*(pars['Oh2']-Oh2_mu)**2/(2*Oh2_sigma**2)+sp.log(int(pars['Oh2']>0))]
    #Bsmu_mu,Bsmu_sigma=[3.2*10**(-9),1.5*10**(-9)]
    #constraints['micromegas']+=[lambda pars: -(pars['Bsmumu']-Bsmu_mu)**2/(2*Bsmu_sigma**2)]
    #bsg_mu,bsg_sigma=[3.55*10**(-4) , 0.42*10**(-4)]
    #constraints['micromegas']+=[lambda pars:  -(pars['bsg']-bsg_mu)**2/(2*bsg_sigma**2)]
    
    #LHC event generation defaults
    #modes: SUSYMG (all processes), SUSY (internal pythia, customizable processes)
    evgen={'mode':'SUSYMG','nevt':10000,'ecm':8000.,'threads':1,'nlo':1}
    #read in process dictionary
    #evgen['proc_dict']=sl.read_proccard(os.path.join(sl.libpath,'src','cards','pyprocs_%s.txt'%evgen['mode']))
    #wether or not to keep output files
    keep_files={'slha':False,'mgcard':False,'hepmc':False,'root':False,'prospino':False}
    #input files
    inputfiles={'susyhit':os.path.join(sl.libpath,'src','cards','susyhit.in'),
                'pythia':os.path.join(sl.libpath,'src','cards','pycard_default.txt'),
                'prospino':os.path.join(sl.libpath,'src','cards','Pro2_subroutines'),
                'delphes':os.path.join(sl.libpath,'src','cards','delcard_default.dat')}

    #Initialization class method (applies to all processes)
    @classmethod
    def initialize(cls,opt):
        #list of allowed model changes
        allowed_changes=['parts','constraints','evgen','keep_files','inputfiles']
        for change in opt.model_change.keys():
            if change in allowed_changes:
                #Update dicts
                if change in ['keep_files','evgen','inputfiles']:
                    exec('cls.%s.update(opt.model_change[change])'%change)
                else:            
                    exec('cls.%s=opt.model_change[change]'%change)
                #update obs_params if parts are changed               
                if change=='parts':
                    cls.obs_params=list(OrderedDict.fromkeys(sum([cls.part_params[key] for key in cls.parts],[])))
            
        #set up file dirs
        for filetype in cls.keep_files.keys():
            if cls.keep_files[filetype]:
                dat_dir=os.path.join(cls.datadirs['main'],filetype)
                cls.datadirs[filetype]=dat_dir
    
    #Read covariance from file        
    def __init__(self,modelid,params):
        self.modelid=modelid
        self.weight=0
        self.lnP=0
        #Set up error dictionary
        self.error=0
        self.accept=True
        self.slhafile='%s.slha'%(self.modelid)
        self.mgcard='%s_mg.dat'%(self.modelid)
        self.hepmcfiles=['%s_%s.hepmc'%(self.modelid,i) for i in range(self.evgen['threads'])]
        self.profile='%s_pro.dat'%(self.modelid)
        self.rootfile='%s.root'%(self.modelid)
        
    #clean up method
    def finalize(self):
        #move/delete slha
        if os.path.exists(self.slhafile):
            if self.accept:
                if self.keep_files['slha']:
                    shutil.copy(self.slhafile,os.path.join(self.datadirs['slha'],self.slhafile))
                    
                os.remove(self.slhafile)
            else:
                os.remove(self.slhafile)
                
        #move/delete mgcard        
        if os.path.exists(self.mgcard):
            if self.accept:
                if self.keep_files['mgcard']:
                    shutil.copy(self.mgcard,os.path.join(self.datadirs['mgcard'],self.mgcard))

                os.remove(self.mgcard)
            else:
                os.remove(self.mgcard)
        
        #move/delete hepmc files
        for hepmcfile in self.hepmcfiles:        
            if os.path.exists(hepmcfile):
                if self.accept:
                    if self.keep_files['hepmc']:
                        shutil.copy(hepmcfile,os.path.join(self.datadirs['hepmc'],hepmcfile))
                        
                    os.remove(hepmcfile)

                else:
                    os.remove(hepmcfile)
                    
        #move/delete prospino
        if os.path.exists(self.profile):
            if self.accept:
                if self.keep_files['prospino']:
                    shutil.copy(self.rootfile,os.path.join(self.datadirs['prospino'],self.profile))
                    
                os.remove(self.profile)
            else:
                os.remove(self.profile)
                
        #move/delete root
        if os.path.exists(self.rootfile):
            if self.accept:
                if self.keep_files['root']:
                    shutil.copy(self.rootfile,os.path.join(self.datadirs['root'],self.rootfile))
                    
                os.remove(self.rootfile)
            else:
                os.remove(self.rootfile)
                
        #merge Datafile
                
            
    def calculate(self,part):
        
        #check if spectrum is OK for further calculation else return
        if self.error>0:
            #print self.errors
            self.accept=False
            self.lnP=-sp.inf
            #remove files created by tools
            if os.path.exists(self.slhafile):
                os.remove(self.slhafile)
            if os.path.exists(self.slhafile+'.fh'):
                os.remove(self.slhafile+'.fh')
            if os.path.exists(self.mgcard):
                os.remove(self.mgcard)
            for hepmcfile in self.hepmcfiles:
                if os.path.exists(hepmcfile):
                    os.remove(hepmcfile)
            if os.path.exists(self.profile):
                os.remove(self.profile)
            if os.path.exists(self.rootfile):
                os.remove(self.rootfile)
                
            return        
        
        #Calculate likelihood components related to tool
        if 'softsusy' in part:
            self=sl.run_softsusy(self)
        if 'susyhit' in part: 
            self=sl.run_susyhit(self)
        if 'feynhiggs' in part:
            self=sl.run_feynhiggs(self)
        if 'higgsbounds' in part:
            self=sl.run_higgsbounds(self)
        if 'micromegas' in part:
            self=sl.run_micromegas(self)
        if 'pythia' in part:
            self=sl.run_pythia(self)
        if 'delphes' in part:
            self=sl.run_delphes(self)
        if 'prospino' in part:
            self=sl.run_prospino(self)
            
        #Calculate likelihood components related to tool               
        for constraint in self.constraints[part]:
            try:
                self.lnP+=constraint(self.params)
            except:
                print "Error in %s constraints... terminating"%(part)
                sys.exit()
            #print self.lnP,constraint(self.params)
            
        return
        
            
    
