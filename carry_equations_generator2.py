# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 10:44:19 2015

@author: Richard
"""

import math

import sympy
from semiprime_tools import num_to_factor_num_qubit
from sympy_helper_fns import max_value, is_equation
from sympy_subs import subs_many

def generate_carry_equations(num_dig1=None, num_dig2=None, product=None):
    ''' Generate the carry equations for a given factorisation

        >>> product = 25
        >>> eqns = generate_carry_equations(product=product)
        >>> for e in eqns: print e
        p1 + q1 == 2*z12
        p1*q1 + z12 + 2 == 2*z23 + 4*z24
        p1 + q1 + z23 == 2*z34 + 1
        z24 + z34 + 1 == 2*z45 + 1
        z45 == 0

        >>> product = 143
        >>> eqns = generate_carry_equations(product=product)
        >>> for e in eqns: print e
        p1 + q1 == 2*z12 + 1
        p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
        p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
        p1 + p2*q2 + q1 + z24 + z34 == 2*z45 + 4*z46
        p2 + q2 + z35 + z45 == 2*z56 + 4*z57
        z46 + z56 + 1 == 2*z67
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
        if pi in [0, num_dig1 - 1]:
            pi_str = '1'
        else:
            pi_str = 'p{}'.format(pi)
        for qi in xrange(num_dig2):
            if qi in [0, num_dig2 - 1]:
                qi_str = '1'
            else:
                qi_str = 'q{}'.format(qi)

            pq_str = '*'.join([pi_str, qi_str])
            eqns_lhs[pi + qi] += sympy.sympify(pq_str)

    ## Now loop over and add the carry variables
    for column_ind, sum_ in enumerate(eqns_lhs):
        if sum_ == 0:
            max_val = 1
        else:
            max_val = max_value(sum_)
        max_pow_2 = int(math.floor(math.log(max_val, 2)))
        for i in xrange(1, max_pow_2 + 1):
            z = sympy.Symbol('z{}{}'.format(column_ind, column_ind + i))
            eqns_rhs[column_ind] += (2 ** i) * z
            eqns_lhs[column_ind + i] += z

    eqns = [sympy.Eq(lhs, rhs) for lhs, rhs in zip(eqns_lhs, eqns_rhs)]
    eqns = filter(is_equation, eqns)

    return eqns

def generate_carry_equations_raw(num_dig1=None, num_dig2=None, product=None):
    ''' Generate the carry equations for a given factorisation

        >>> product = 25
        >>> eqns = generate_carry_equations_raw(product=product)
        >>> for e in eqns: print e
        p0*q0 == 1
        p0*q1 + p1*q0 == 2*z12
        p0*q2 + p1*q1 + p2*q0 + z12 == 2*z23 + 4*z24
        p1*q2 + p2*q1 + z23 == 2*z34 + 1
        p2*q2 + z24 + z34 == 2*z45 + 1
        z45 == 0

        >>> product = 143
        >>> eqns = generate_carry_equations_raw(product=product)
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
        pi_str = 'p{}'.format(pi)
        for qi in xrange(num_dig2):
            qi_str = 'q{}'.format(qi)

            pq_str = '*'.join([pi_str, qi_str])
            eqns_lhs[pi + qi] += sympy.sympify(pq_str)

    ## Now loop over and add the carry variables
    for column_ind, sum_ in enumerate(eqns_lhs):
        if sum_ == 0:
            max_val = 1
        else:
            max_val = max_value(sum_)
        max_pow_2 = int(math.floor(math.log(max_val, 2)))
        for i in xrange(1, max_pow_2 + 1):
            z = sympy.Symbol('z{}{}'.format(column_ind, column_ind + i))
            eqns_rhs[column_ind] += (2 ** i) * z
            eqns_lhs[column_ind + i] += z

    eqns = [sympy.Eq(lhs, rhs) for lhs, rhs in zip(eqns_lhs, eqns_rhs)]
    eqns = filter(is_equation, eqns)

    return eqns


def generate_factorisation_equation(num_dig1=None, num_dig2=None, product=None):
    ''' Generate the carry equations for a given factorisation

        >>> product = 9
        >>> eqn = generate_factorisation_equation(product=product)
        >>> print eqn
        p0*q0 + 2*p0*q1 + 2*p1*q0 + 4*p1*q1 == 9

        >>> product = 143
        >>> eqn = generate_factorisation_equation(product=product)
        >>> for i in xrange(0, len(str(eqn)), 80): print str(eqn)[i:i+80]
        p0*q0 + 2*p0*q1 + 4*p0*q2 + 8*p0*q3 + 2*p1*q0 + 4*p1*q1 + 8*p1*q2 + 16*p1*q3 + 4
        *p2*q0 + 8*p2*q1 + 16*p2*q2 + 32*p2*q3 + 8*p3*q0 + 16*p3*q1 + 32*p3*q2 + 64*p3*q
        3 == 143
    '''
    if product is None:
        raise ValueError('generate_carry_equations must be given a product')

    if num_dig1 is None:
        assert num_dig2 is None
        num_dig1, num_dig2 = num_to_factor_num_qubit(product)

    eqn_rhs = sympy.sympify(int(bin(product), 2))
    eqn_lhs = 0

    ## Now add the contributions from the actual factors
    for pi in xrange(num_dig1):
        for qi in xrange(num_dig2):
            pq_str = 'p{} * q{}'.format(pi, qi)
            eqn_lhs += sympy.sympify(pq_str) * 2 ** (pi + qi)

    return sympy.Eq(eqn_lhs, eqn_rhs)

def generate_carry_equations_auxiliary(num_dig1=None, num_dig2=None, product=None):
    ''' Given a product, generate the carry equations that express this
        using auxiliary variables for the p_i * q_j interactions

        >>> product = 25
        >>> eqns = generate_carry_equations_auxiliary(product=product)
        >>> for e in eqns: print e
        m0_0 == 1
        m0_1 + m1_0 == 2*z12
        m0_2 + m1_1 + m2_0 + z12 == 2*z23 + 4*z24
        m1_2 + m2_1 + z23 == 2*z34 + 1
        m2_2 + z24 + z34 == 2*z45 + 1
        z45 == 0
        p0*q0 == m0_0
        p0*q1 == m0_1
        p0*q2 == m0_2
        p1*q0 == m1_0
        p1*q1 == m1_1
        p1*q2 == m1_2
        p2*q0 == m2_0
        p2*q1 == m2_1
        p2*q2 == m2_2

        >>> product = 143
        >>> eqns = generate_carry_equations_auxiliary(product=product)
        >>> for e in eqns: print e
        m0_0 == 1
        m0_1 + m1_0 == 2*z12 + 1
        m0_2 + m1_1 + m2_0 + z12 == 2*z23 + 4*z24 + 1
        m0_3 + m1_2 + m2_1 + m3_0 + z23 == 2*z34 + 4*z35 + 1
        m1_3 + m2_2 + m3_1 + z24 + z34 == 2*z45 + 4*z46
        m2_3 + m3_2 + z35 + z45 == 2*z56 + 4*z57
        m3_3 + z46 + z56 == 2*z67
        z57 + z67 == 2*z78 + 1
        z78 == 0
        p0*q0 == m0_0
        p0*q1 == m0_1
        p0*q2 == m0_2
        p0*q3 == m0_3
        p1*q0 == m1_0
        p1*q1 == m1_1
        p1*q2 == m1_2
        p1*q3 == m1_3
        p2*q0 == m2_0
        p2*q1 == m2_1
        p2*q2 == m2_2
        p2*q3 == m2_3
        p3*q0 == m3_0
        p3*q1 == m3_1
        p3*q2 == m3_2
        p3*q3 == m3_3
    '''
    if product is None:
        raise ValueError('generate_carry_equations must be given a product')
    if num_dig1 is None:
        assert num_dig2 is None
        num_dig1, num_dig2 = num_to_factor_num_qubit(product)

    eqns_rhs = [int(digit) for digit in bin(product)[2:][::-1]]
    eqns_lhs = [0 for _ in eqns_rhs]

    # Create a holder for the equations that constrain pi.qj=mi_j
    constraints = []

    # Now pad them
    for i in xrange(5):
        eqns_lhs.append(0)
        eqns_rhs.append(0)

    ## Now add the contributions from the actual factors
    for pi in xrange(num_dig1):
        pi_str = 'p{}'.format(pi)
        for qi in xrange(num_dig2):
            qi_str = 'q{}'.format(qi)

            # Add the single interaction variable
            interaction_str = 'm{}_{}'.format(pi, qi)
            interaction_var = sympy.sympify(interaction_str)
            eqns_lhs[pi + qi] += interaction_var

            constraint = sympy.Eq(sympy.sympify('{}*{}'.format(pi_str, qi_str)),
                                  interaction_var)
            constraints.append(constraint)


    ## Now loop over and add the carry variables
    for column_ind, sum_ in enumerate(eqns_lhs):
        if sum_ == 0:
            max_val = 1
        else:
            max_val = max_value(sum_)
        max_pow_2 = int(math.floor(math.log(max_val, 2)))
        for i in xrange(1, max_pow_2 + 1):
            z = sympy.Symbol('z{}{}'.format(column_ind, column_ind + i))
            eqns_rhs[column_ind] += (2 ** i) * z
            eqns_lhs[column_ind + i] += z

    eqns = [sympy.Eq(lhs, rhs) for lhs, rhs in zip(eqns_lhs, eqns_rhs)]
    eqns = filter(is_equation, eqns)

    return eqns + constraints
if __name__ == "__main__":
    import doctest
    doctest.testmod()