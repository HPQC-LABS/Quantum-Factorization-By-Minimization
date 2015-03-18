# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 21:09:13 2015

@author: Richard
"""

import itertools
import multiprocessing

from sympy.core.cache import clear_cache

DEFAULT_NUM_PROCESSES = 3

def get_pool(processes=DEFAULT_NUM_PROCESSES):
    ''' Return a pool object '''
    return multiprocessing.Pool(processes=processes)

def _subs_wrap((eqns, sub_tuple)):
    ''' Wrapper function to allow parallelisation '''
    subbed = [eqn.subs(sub_tuple) for eqn in eqns]
    clear_cache()
    return subbed
    

def paralellised_subs(equations, to_sub, processes=DEFAULT_NUM_PROCESSES, pool=None):
    ''' Take a list of equations and a substitution dict and parallelise the
        substitution.
        If a pool is given, ignore processes

        >>> x, y, z = sympy.symbols('x y z')
        >>> eqns = [x + y - 1,
        ...         x*z - 1,
        ...         x - 1]
        >>> eqns = map(balance_terms, map(sympy.Eq, eqns))
        >>> eqns = [(eqn,) for eqn in eqns]
        >>> subs = {x: 1, z: 2}
        >>> subbed = paralellised_subs(eqns, subs)
        >>> for eqn in subbed: print eqn
        [y + 1 == 1]
        [False]
        [True]
    '''
    if isinstance(to_sub, dict):
        to_sub = to_sub.items()
    
    if pool is None:
        pool = get_pool(processes=processes)
        close_on_exit = True
    else:
        close_on_exit = False

    sub_copies = itertools.repeat(to_sub, len(equations))
    result = pool.map(_subs_wrap, zip(equations, sub_copies))

    # If we made the pool, shut it down again
    if close_on_exit:
        pool.close()
        pool.join()
        pool = None
    return result
    
if __name__ == "__main__":
    import doctest
    import sympy
    from sympy_helper_fns import balance_terms
    doctest.testmod()