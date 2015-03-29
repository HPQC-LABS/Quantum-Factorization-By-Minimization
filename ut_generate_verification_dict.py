# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 14:23:12 2015

@author: Richard
"""

import itertools
import cPickle

from cfg_sympy_solver import (SEMIPRIME_FILENAME, FACTOR_DICT_FILENAME, 
                              FACTORS_FILENAME, SEMIPRIMES_HAMMING_TEMPLATE)


def get_factors_large_semiprimes():
    factor_dict = {}

    # First get the large semiprimes
    sp_file = open(SEMIPRIME_FILENAME, 'r')
    f_file = open(FACTORS_FILENAME, 'r')
    
    for semiprime, factors in itertools.izip(sp_file, f_file):
        semiprime = int(semiprime)
        f1, f2 = map(int, factors.split(', '))    
        
        if f1 > f2:
            temp = f1
            f1 = f2
            f2 = temp
        
        check_factorisation(semiprime, f1, f2)
            
        factor_dict[semiprime] = (f1, f2)
    return factor_dict


def get_factors_hamming(dim='20'):
    factors = {}

    dim = str(dim)
    
    filename = SEMIPRIMES_HAMMING_TEMPLATE.format(dim, dim)
    exps = open(filename, 'r').read()

    # Turn the text into tuples
    exps = exps.split('\n')
    exps = filter(None, exps)

    for exp in exps:
        prod, dist, p, q = map(int, exp.strip('(').strip(')').split(', '))

        if p > q:
            temp = p
            p = q
            q = temp

        check_factorisation(prod, p, q)
        factors[prod] = (p, q)
    return factors

def check_factorisation(prod, f1, f2):
    assert f1 * f2 == prod
    assert f1 <= f2   

DIMENSIONS = range(20, 530, 10)

factor_dict = {}

factor_dict.update(get_factors_large_semiprimes())

for dim in DIMENSIONS:
    factor_dict.update(get_factors_hamming(dim))

dict_file = open(FACTOR_DICT_FILENAME, 'w')
cPickle.dump(factor_dict, dict_file)
dict_file.close()