# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 14:23:12 2015

@author: Richard
"""

import itertools
import cPickle

from cfg_sympy_solver import SEMIPRIME_FILENAME, FACTOR_DICT_FILENAME, FACTORS_FILENAME

factor_dict = {}

sp_file = open(SEMIPRIME_FILENAME, 'r')
f_file = open(FACTORS_FILENAME, 'r')

for semiprime, factors in itertools.izip(sp_file, f_file):
    semiprime = int(semiprime)
    f1, f2 = map(int, factors.split(', '))    
    
    if f1 > f2:
        temp = f1
        f1 = f2
        f2 = temp
    
    assert f1 * f2 == semiprime    
    assert f1 <= f2    
        
    factor_dict[semiprime] = (f1, f2)

dict_file = open(FACTOR_DICT_FILENAME, 'w')
cPickle.dump(factor_dict, dict_file)
dict_file.close()