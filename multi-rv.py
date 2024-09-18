#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 14:23:26 2024

@author: heitor
"""

import os
import numpy as np
import sys



#---------------




folder='../out/helio/'


#os.system('ls '+folder+'*.csv > sample_helio.txt')
specs = np.loadtxt(folder+"sample_helio.txt", dtype=str) 


for s in specs:
    print('------')


    print(s)
    
    os.system('python3  rv-calc.py '+ s)


    print('------')



































#