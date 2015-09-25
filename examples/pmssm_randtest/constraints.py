# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 17:15:45 2015

@author: Jan Lindroos
Replica of ATLAS pMSSM constraints 
"""
import scipy as sp

#The different likelihood components for SUSY
parts=['softsusy','susyhit','higgsbounds','micromegas','pythia','delphes','prospino']          
#import constraints
                                                      
#Constraints: construct components of lnP****************************************************************************                  
constraints={part:[] for part in parts}
#softsusy: Gamma_Z, SUSY masses
GammaZ_max=2*10**(-3)
constraints['softsusy']+=[lambda pars: sp.log(int(pars['Gamma_Z']<GammaZ_max))]
#Higgsbounds:'HB_excl','m_h0','dm_h0','m_A0','dm_A0','m_H0','dm_H0','m_Hp','dm_Hp'
#HB exclusion
constraints['higgsbounds']+=[lambda pars: sp.log(int(pars['HB_excl']==0))]
#Higgs mass: (limited by theory uncertainty. 3-loop could improve this ATLAS)
mh0_range=[123.,129.]
constraints['higgsbounds']+=[lambda pars:sp.log(int(pars['mh0']>mh0_range[0] and pars['mh0']<mh0_range[1]))]

#Micromegas:'LEP','Oh2','Oh2_bino','Oh2_hsino1','Oh2_hsino2','Oh2_wino','a_mu','drho','bsg','Bsmumu','Btaunu','sigma_nSD', 'sigma_nSI', 'sigma_pSD', 'sigma_pSI'
#LEP limits
constraints['micromegas']+=[lambda pars: sp.log(int(pars['LEP']==0))]
#Relic density (5 sigma upper limit combined Planck+10% theory error: arXiv:1303.5076)
Oh2_max=0.1248 
constraints['micromegas']+=[lambda pars: sp.log(int(pars['Oh2']<Oh2_max))]
#delta_rho (arXiv:1209.2716)
drho_range=[-0.0005,0.0017]
constraints['micromegas']+=[lambda pars: sp.log(int(pars['drho']>drho_range[0] and pars['drho']<drho_range[1]))]
#Bsmumu (microMegas corrected by 8.8% per arXiv:1204.1737)
Bsmumu_range=[1.1*10**(-9),6.4*10**(-9)]
corr=1
constraints['micromegas']+=[lambda pars: sp.log(int(pars['Bsmumu']*corr>Bsmumu_range[0] and pars['Bsmumu']*corr<Bsmumu_range[1]))]
#bsg (union of 2sigma exp and theory, no ref)
bsg_range=[0.269*10**(-3),0.387*10**(-3) ]
constraints['micromegas']+=[lambda pars: sp.log(int(pars['bsg']>bsg_range[0] and pars['bsg']<bgs_range[1]))]
