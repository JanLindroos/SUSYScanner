# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 17:50:29 2014
@author: Jan øye Lindroos

This is the main file which parses the user options
and launches the scan based on chioce of paralellism
"""
import sys
from lib.prepost_lib import init_run, clean_run
import lib.para_lib as prl
from alg.scan import scan

if __name__=='__main__':
    #let all workers parse user provided options
    alg,model,opt=init_run(sys.argv)
    print 'after init_run'
    #Launch parallell scan based on options
    #opt='hello'
    print 'before launch'
    print dir(model)
    prl.launch(scan,alg,model,opt)

    #Clean up after run, empty tmp dir etc
    clean_run(opt)               