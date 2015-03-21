# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 12:26:24 2015

@author: Richard
"""

import sympy

from contradiction_exception import ContradictionException
from sympy_helper_fns import min_value, max_value, parity

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

if __name__ == "__main__":
    import doctest
    doctest.testmod()