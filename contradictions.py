# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 12:26:24 2015

@author: Richard
"""

import itertools

import sympy

from contradiction_exception import ContradictionException
from sympy_helper_fns import min_value, max_value, parity, str_eqns_to_sympy_eqns
from sympy_subs import subs

def apply_contradictions(equations):
    ''' Now look for contradictions in the equations '''
    for eqn in equations:
        contradiction_1(eqn)
        contradiction_2(eqn)

## Look for contradictions
def contradiction_1(eqn):
    ''' Check the values could be equal 
    
    >>> x, y, z = sympy.symbols('x y z')

    >>> eqn = sympy.Eq(x*y*z, 2)
    >>> contradiction_1(eqn)
    Traceback (most recent call last):
        ...
    ContradictionException: contradiction_1: x*y*z == 2


    >>> eqn = sympy.Eq(x*y*z + 1)
    >>> contradiction_1(eqn)
    Traceback (most recent call last):
        ...
    ContradictionException: contradiction_1: x*y*z + 1 == 0
    '''
    def _helper(lhs, rhs):
        if min_value(lhs) > max_value(rhs):
            raise ContradictionException('contradiction_1: {}'.format(eqn))

    _helper(eqn.lhs, eqn.rhs)
    _helper(eqn.rhs, eqn.lhs)

def contradiction_2(eqn):
    ''' Check the parity 
    
    >>> x, y, z = sympy.symbols('x y z')
    
    >>> eqn = sympy.Eq(2*x*y + 4*z + 1)
    >>> contradiction_2(eqn)
    Traceback (most recent call last):
        ...
    ContradictionException: contradiction_2: 2*x*y + 4*z + 1 == 0
    
    >>> eqn = sympy.Eq(2*x*y * 4*z - 2)
    >>> contradiction_2(eqn)

    >>> eqn = sympy.Eq(2*x*y + z - 1)
    >>> contradiction_2(eqn)
    '''
    l_parity = parity(eqn.lhs)
    if l_parity is not None:
        r_parity = parity(eqn.rhs)
        if (r_parity is not None) and (l_parity != r_parity):
            raise ContradictionException('contradiction_2: {}'.format(eqn))

def contradiction_mini_assump(eqn, atom_limit=4):
    ''' If we have a small equation, see if any valid solutions exist at all
    
        >>> eqns = ['x + y == 3', 'x + x*y + y == 2', 'x + y == 2', 'x*y == 1']
        >>> eqns = str_eqns_to_sympy_eqns(eqns)
        
        >>> for e in eqns:
        ...     try:
        ...         contradiction_mini_assump(e, atom_limit=3)
        ...         print '{} is fine'.format(e)
        ...     except ContradictionException as contradiction:
        ...         print contradiction
        No valid combinations for x + y == 3
        No valid combinations for x*y + x + y == 2
        x + y == 2 is fine
        x*y == 1 is fine
    '''
    atoms = list(eqn.atoms(sympy.Symbol))
    if len(atoms) > atom_limit:
        return
    
    vals = itertools.product(range(2), repeat=len(atoms))
    for val in vals:
        to_sub = dict(zip(atoms, val))
        evaluated = subs(eqn, to_sub)
        if evaluated:
            return
    
    raise ContradictionException('No valid combinations for {}'.format(eqn))

if __name__ == "__main__":
    import doctest
    doctest.testmod()