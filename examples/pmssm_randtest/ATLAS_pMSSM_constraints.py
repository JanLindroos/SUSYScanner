# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 17:15:45 2015

@author: Jan Lindroos
Replica of ATLAS pMSSM constraints 
"""
import scipy as sp
import os
import cPickle as pickle
                                                      
#Constraints: construct components of lnP****************************************************************************                  
constraints=dict([(part,[]) for part in ['softsusy','susyhit','higgsbounds','micromegas','pythia','delphes','prospino']])
#softsusy: 'Gamma_Z', SUSY masses
#GammaZ (Calculated at Tree level)
GammaZ_max=2*10**(-3)
constraints['softsusy']+=[lambda pars: sp.log(int(pars['Gamma_Z']<GammaZ_max))]
#Higgsbounds:'HB_excl','m_h0','dm_h0','m_A0','dm_A0','m_H0','dm_H0','m_Hp','dm_Hp'
#HB exclusion
constraints['higgsbounds']+=[lambda pars: sp.log(int(pars['HB_excl']==0))]
#Higgs mass: (limited by theory uncertainty. 3-loop could improve this ATLAS)
mh0_range=[123.,129.]
constraints['higgsbounds']+=[lambda pars:sp.log(int(pars['m_h0']>mh0_range[0] and pars['m_h0']<mh0_range[1]))]

#Micromegas:'LEP','Oh2','Oh2_bino','Oh2_hsino1','Oh2_hsino2','Oh2_wino','a_mu','drho','bsg','Bsmumu','Btaunu','sigma_nSD', 'sigma_nSI', 'sigma_pSD', 'sigma_pSI'
#LEP limits
constraints['micromegas']+=[lambda pars: sp.log(int(pars['LEP']==0))]
#Relic density (5 sigma upper limit combined Planck+10% theory error: arXiv:1303.5076)
Oh2_max=0.1208 
constraints['micromegas']+=[lambda pars: sp.log(int(pars['Oh2']<Oh2_max))]
#delta_rho (arXiv:1209.2716)
drho_range=[-0.0005,0.0017]
constraints['micromegas']+=[lambda pars: sp.log(int(pars['drho']>drho_range[0] and pars['drho']<drho_range[1]))]
#Bsmumu (microMegas corrected by 8.8% per arXiv:1204.1737)
Bsmumu_range=[1.1*10**(-9),6.4*10**(-9)]
corr=1/(1-0.088)
constraints['micromegas']+=[lambda pars: sp.log(int(pars['Bsmumu']*corr>Bsmumu_range[0] and pars['Bsmumu']*corr<Bsmumu_range[1]))]
#Btaunu (union of 2sigma exp and theory, 1211.1976, 1210.8443, Nakao ICHEP 2012, and CKM Fitter )
Btaunu_range=[64*10**(-6),161*10**(-6)]
constraints['micromegas']+=[lambda pars: sp.log(int(pars['Btaunu']>Btaunu_range[0] and pars['Btaunu']<Btaunu_range[1]))]
#bsg (union of 2sigma exp and theory, no ref)
bsg_range=[0.269*10**(-3),0.387*10**(-3) ]
constraints['micromegas']+=[lambda pars: sp.log(int(pars['bsg']>bsg_range[0] and pars['bsg']<bsg_range[1]))]
#a_mu (union of 3sigma exp and theory, no ref)
amu_range=[-1.77*10**(-9),4.38*10**(-9)]
constraints['micromegas']+=[lambda pars: sp.log(int(pars['a_mu']>amu_range[0] and pars['a_mu']<amu_range[1]))]
#LUX limit (4*LUX lim on relic density corrected spin-independent Xenon (N=124, A=54) cs
LUX_path=os.path.join(os.path.dirname(__file__),'LUX_lim.pickle')
LUX_lim=pickle.load(open(LUX_path,'rb'))
constraints['micromegas']+=[lambda pars: sp.log(int(pars['Oh2']/float(0.1188)*(54/float(124)*pars['sigma_pSI']+70/float(124)*pars['sigma_nSI'])<4*LUX_lim(pars['m_Neu1'])))]
