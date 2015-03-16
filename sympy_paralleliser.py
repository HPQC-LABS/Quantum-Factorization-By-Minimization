# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 21:09:13 2015

@author: Richard
"""

import itertools
import multiprocessing

DEFAULT_PROCESSES = 3

# Minimum number of equations needed before parallelisation is used.
PARALLELISED_SUBS_EQN_LIMIT = 30

def get_pool(processes=DEFAULT_PROCESSES):
    ''' Return a pool object '''
    return multiprocessing.Pool(processes=processes)

def _subs_wrap((eqn, sub_tuple)):
    ''' Wrapper function to allow parallelisation '''
    return eqn.subs(sub_tuple)

def paralellised_subs(equations, to_sub, processes=DEFAULT_PROCESSES, pool=None):
    ''' Take a list of equations and a substitution dict and parallelise the
        substitution.
        If a pool is given, ignore processes
    '''
    if len(equations) < PARALLELISED_SUBS_EQN_LIMIT:
        return [eqn.subs(to_sub) for eqn in equations]

    if isinstance(to_sub, dict):
        to_sub = to_sub.items()
    
    if pool is None:
        pool = get_pool(processes=processes)
    sub_copies = itertools.repeat(to_sub, len(equations))
    result = pool.map(_subs_wrap, zip(equations, sub_copies))
    return result