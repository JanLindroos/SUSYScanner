# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 17:15:45 2015

@author: Jan Lindroos
Custom SUSY Constraints file
"""
import scipy as sp

#The different likelihood components for SUSY
parts=['softsusy','susyhit','feynhiggs','higgsbounds','micromegas','pythia','delphes','prospino']          
#import constraints
                                                      
#Constraints: construct components of lnP****************************************************************************                  
constraints=dict([(part,[]) for part in parts])

constraints['higgsbounds']+=[lambda pars: sp.log(int(pars['HB_excl']==0))]
mh0_mu,mh0_sigma=[125.5 , 1.7]
constraints['higgsbounds']+=[lambda pars: -(pars['m_h0']-mh0_mu)**2/(2*mh0_sigma**2)]
constraints['micromegas']+=[lambda pars: sp.log(int(pars['LEP']==0))]
Oh2_mu,Oh2_sigma=[0.1199 , 0.0027]
constraints['micromegas']+=[lambda pars: -int(pars['Oh2']>Oh2_mu)*(pars['Oh2']-Oh2_mu)**2/(2*Oh2_sigma**2)+sp.log(int(pars['Oh2']>0))]
Bsmu_mu,Bsmu_sigma=[3.2*10**(-9),1.5*10**(-9)]
constraints['micromegas']+=[lambda pars: -(pars['Bsmumu']-Bsmu_mu)**2/(2*Bsmu_sigma**2)]
bsg_mu,bsg_sigma=[3.55*10**(-4) , 0.42*10**(-4)]
constraints['micromegas']+=[lambda pars:  -(pars['bsg']-bsg_mu)**2/(2*bsg_sigma**2)]
