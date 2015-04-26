# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 10:44:19 2015

@author: Richard
"""

import math

import sympy
from semiprime_tools import num_to_factor_num_qubit
from sympy_helper_fns import num_add_terms, is_equation

def generate_carry_equations(num_dig1=None, num_dig2=None, product=None):
    ''' Generate the carry equations for a given factorisation 
    
        >>> product = 143
        >>> eqns = generate_carry_equations(product=product)
        >>> for e in eqns: print e
        p0*q0 == 1
        p0*q1 + p1*q0 == 2*z12 + 1
        p0*q2 + p1*q1 + p2*q0 + z12 == 2*z23 + 4*z24 + 1
        p0*q3 + p1*q2 + p2*q1 + p3*q0 + z23 == 2*z34 + 4*z35 + 1
        p1*q3 + p2*q2 + p3*q1 + z24 + z34 == 2*z45 + 4*z46
        p2*q3 + p3*q2 + z35 + z45 == 2*z56 + 4*z57
        p3*q3 + z46 + z56 == 2*z67
        z57 + z67 == 2*z78 + 1
        z78 == 0
    '''
    if product is None:
        raise ValueError('generate_carry_equations must be given a product')
    if num_dig1 is None:
        assert num_dig2 is None
        num_dig1, num_dig2 = num_to_factor_num_qubit(product)
    
    eqns_rhs = [int(digit) for digit in bin(product)[2:][::-1]]
    eqns_lhs = [0 for _ in eqns_rhs]
    
    # Now pad them
    for i in xrange(5):
        eqns_lhs.append(0)
        eqns_rhs.append(0)
    
    ## Now add the contributions from the actual factors
    for pi in xrange(num_dig1):
        for qi in xrange(num_dig2):
            pq_str = 'p{} * q{}'.format(pi, qi)
            eqns_lhs[pi + qi] += sympy.sympify(pq_str)
    
    ## Now loop over and add the carry variables
    for column_ind, sum_ in enumerate(eqns_lhs):
        num_terms = num_add_terms(sum_)
        max_pow_2 = int(math.floor(math.log(num_terms, 2)))
        for i in xrange(1, max_pow_2 + 1):
            z = sympy.Symbol('z{}{}'.format(column_ind, column_ind + i))
            eqns_rhs[column_ind] += (2 ** i) * z
            eqns_lhs[column_ind + i] += z
    
    eqns = [sympy.Eq(lhs, rhs) for lhs, rhs in zip(eqns_lhs, eqns_rhs)]
    
    eqns = filter(is_equation, eqns)
    return eqns

if __name__ == "__main__":
    import doctest
    doctest.testmod()