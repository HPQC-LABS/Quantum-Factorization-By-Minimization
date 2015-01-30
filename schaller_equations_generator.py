# -*- coding: utf-8 -*-
"""
Created on Sun Jan 25 17:03:10 2015

@author: Richard
"""

import sympy

from sympy_helper_fns import remove_binary_squares, balance_terms, is_equation

def _sympify_dict(dict_):
    ''' Given a dictionary, build up a new dict that is sympified '''
    out = {}
    for k, v in dict_.iteritems():
        sym_k = sympy.sympify(k)
        sym_v = sympy.sympify(v)
        
        if out.get(sym_k) is not None:
            raise ValueError('Duplicate keys found upon sympification')
        out[sym_k] = sym_v
    return out

def schaller_transform(a, b, s):
    ''' Take a single 2-qubit interaction and a sum and return the new
        Hamiltonian

        >>> a, b, s, s1, s2 = sympy.symbols('a b s s1 s2')
        >>> schaller_transform(a, b, s)
        a*b + 2*a*s + 2*b*s + s

        >>> schaller_transform(a, b, s1 + s2)
        a*b + 2*a*s1 + 2*a*s2 + 2*b*s1 + 2*b*s2 + 4*s1*s2 + s1 + s2
    '''
    TWO = sympy.sympify('2')
    sqrd = (a + b - (1 / TWO)) / TWO + s
    sqrd = remove_binary_squares((sqrd ** 2).expand())
    return (2 * sqrd) - (1 / TWO ** 3)

def generate_schaller_equations(num_digits_multiplicand_1,
                                num_digits_multiplicand_2,
                                product):
    ''' Generate the Schaller equations by introducing auxiliary variables 
    
        >>> eqns = generate_schaller_equations(4, 4, 143)
        >>> for e in eqns: print e
        2*S1_0 + 4*S1_1*z1_1 + 5*S1_1 + 5*z1_1 + 1 == 8*S1_0*S1_1 + 8*S1_0*z1_1
        8*S1_1*S2_1 + 6*S1_1 + 2*S1_2*q2 + 4*S1_2*z1_2 + 3*S1_2 + S2_1 + 2*q2*z1_2 + q2 + 3*z1_2 == 8*S1_1*S1_2 + 4*S1_1*q2 + 8*S1_1*z1_2 + 4*S1_2*S2_1 + 2*S2_1*q2 + 4*S2_1*z1_2
        8*S1_2*S2_2 + 6*S1_2 + S2_2 + 2*q3*z1_3 + q3 + 3*z1_3 == 4*S1_2*q3 + 8*S1_2*z1_3 + 2*S2_2*q3 + 4*S2_2*z1_3
        8*S1_3*S2_3 + 2*S1_3 + 1 == S2_3
        2*S2_1*p2 + 4*S2_1*z2_1 + 3*S2_1 + 2*p2*z2_1 + p2 + 6*z1_1 + 3*z2_1 == 8*S2_1*z1_1 + 4*p2*z1_1 + 8*z1_1*z2_1
        2*S2_2*p2 + 2*S2_2*q2 + 4*S2_2*z2_2 + S2_2 + 8*S3_1*z1_2 + 3*S3_1 + p2*q2 + 2*p2*z2_2 + 2*q2*z2_2 + 10*z1_2 + z2_2 == 4*S2_2*S3_1 + 8*S2_2*z1_2 + 2*S3_1*p2 + 2*S3_1*q2 + 4*S3_1*z2_2 + 4*p2*z1_2 + 4*q2*z1_2 + 8*z1_2*z2_2
        2*S2_3*p2 + 2*S2_3*q3 + 4*S2_3*z2_3 + S2_3 + 8*S3_2*z1_3 + 3*S3_2 + p2*q3 + 2*p2*z2_3 + 2*q3*z2_3 + 10*z1_3 + z2_3 == 4*S2_3*S3_2 + 8*S2_3*z1_3 + 2*S3_2*p2 + 2*S3_2*q3 + 4*S3_2*z2_3 + 4*p2*z1_3 + 4*q3*z1_3 + 8*z1_3*z2_3
        S3_3 + p2 == 2*S3_3*p2
        2*S3_1*p3 + 4*S3_1*z3_1 + 3*S3_1 + 2*p3*z3_1 + p3 + 6*z2_1 + 3*z3_1 == 8*S3_1*z2_1 + 4*p3*z2_1 + 8*z2_1*z3_1
        2*S3_2*p3 + 2*S3_2*q2 + 4*S3_2*z3_2 + S3_2 + 8*S4_1*z2_2 + 3*S4_1 + p3*q2 + 2*p3*z3_2 + 2*q2*z3_2 + 10*z2_2 + z3_2 == 4*S3_2*S4_1 + 8*S3_2*z2_2 + 2*S4_1*p3 + 2*S4_1*q2 + 4*S4_1*z3_2 + 4*p3*z2_2 + 4*q2*z2_2 + 8*z2_2*z3_2
        2*S3_3*p3 + 2*S3_3*q3 + 4*S3_3*z3_3 + S3_3 + 8*S4_2*z2_3 + 3*S4_2 + p3*q3 + 2*p3*z3_3 + 2*q3*z3_3 + 10*z2_3 + z3_3 == 4*S3_3*S4_2 + 8*S3_3*z2_3 + 2*S4_2*p3 + 2*S4_2*q3 + 4*S4_2*z3_3 + 4*p3*z2_3 + 4*q3*z2_3 + 8*z2_3*z3_3
        S4_3 + p3 == 2*S4_3*p3
        S4_1 + 10*z3_1 == 8*S4_1*z3_1
        2*S4_2*q2 + 14*z3_2 + 1 == 8*S4_2*z3_2 + S4_2 + 4*q2*z3_2 + q2
        2*S4_3*q3 + 14*z3_3 + 1 == 8*S4_3*z3_3 + S4_3 + 4*q3*z3_3 + q3
    
    '''
    if num_digits_multiplicand_1 < num_digits_multiplicand_2:
        temp = num_digits_multiplicand_1
        num_digits_multiplicand_1 = num_digits_multiplicand_2
        num_digits_multiplicand_2 = temp
    
    target_digits = bin(product)[2:]
    
    n = num_digits_multiplicand_1 + num_digits_multiplicand_2
    k = num_digits_multiplicand_1
    
    if len(target_digits) < n:
        target_digits = target_digits.rjust(n, '0')
    
    equations = []
    edge_subs = {}
    
    # We know p1, q1 are equal to 1 as the numbers must be odd
    edge_subs['p1'] = 1
    edge_subs['q1'] = 1
    
    # Also the largest digit must be 1
    edge_subs['p{k}'.format(k=k)] = 1
    edge_subs['q{nmk}'.format(nmk=n - k)] = 1
        
    # Loop over and create the edge cases
    for i in xrange(1, k + 1):
        for j in xrange(1, n - k + 1):
            # Now add in the edge cases
            if i == 1:
                edge_subs['z0_{j}'.format(j=j)] = 'S1_{jm1}'.format(jm1=j - 1)
            
            if j == 1:
                # We want i - 1 as the digit indices start at 1
                edge_subs['S{i}_0'.format(i=i)] = target_digits[i - 1]
            
            if i == k:
                edge_subs['S{kp1}_{jm1}'.format(kp1=k + 1, jm1=j - 1)] = target_digits[k + j - 1]
            
            if j == n - k:
                edge_subs['S{i}_{nmk}'.format(i=i, nmk=n - k)] = 0
                edge_subs['z{i}_{nmk}'.format(i=i, nmk=n - k)] = 0
            
            if i == k:
                edge_subs['z{k}_{j}'.format(k=k, j=j)] = 0
            
            if i == 1:
                edge_subs['S1_{nmkm1}'.format(nmkm1=n - k - 1)] = 0
    
    for i in xrange(1, k + 1):
        for j in xrange(1, n - k + 1):
            
            formatter = {'i': i, 'j': j, 'ip1': i + 1, 'jm1': j - 1, 
                         'im1': i - 1}

            a = 'p{i}'.format(**formatter)
            a_sub = edge_subs.get(a)
            if a_sub is not None:
                a = a_sub
            a = sympy.sympify(a)
            
            b = 'q{j}'.format(**formatter)
            b_sub = edge_subs.get(b)
            if b_sub is not None:
                b = b_sub
            b = sympy.sympify(b)
            
            s = 0 
            lterms = ['S{i}_{j}', 'z{i}_{j}']
            for term in lterms:
                term = term.format(**formatter)
                sub = edge_subs.get(term)
                if sub is not None:
                    term = sub
                term = sympy.sympify(term)
                s += term
            

            rterms = ['S{ip1}_{jm1}', 'z{im1}_{j}']
            for coef, term in enumerate(rterms):
                term = term.format(**formatter)
                sub = edge_subs.get(term)
                if sub is not None:
                    term = sub
                term = sympy.sympify(term)
                s -= (coef + 1) * term
            
            equations.append(schaller_transform(a, b, s))
    
    #equations = map(remove_binary_squares, equations)
    equations = map(sympy.Eq, equations)
    equations = map(balance_terms, equations)
    equations = filter(is_equation, equations)
    
    return equations

def generate_unschallered_equations(num_digits_multiplicand_1,
                                num_digits_multiplicand_2,
                                product):
    ''' Generate the Schaller equations by introducing auxiliary variables '''
    if num_digits_multiplicand_1 < num_digits_multiplicand_2:
        temp = num_digits_multiplicand_1
        num_digits_multiplicand_1 = num_digits_multiplicand_2
        num_digits_multiplicand_2 = temp
    
    target_digits = bin(product)[2:]
    
    n = num_digits_multiplicand_1 + num_digits_multiplicand_2
    k = num_digits_multiplicand_1
    
    equations = []
    edge_subs = {}
    
    for i in xrange(1, k + 1):
        for j in xrange(1, n - k + 1):
            # Add the multiplicand term
            lhs = 'p{i}*q{j} + S{i}_{j} + z{i}_{j}'.format(i=i, j=j)
            
            rhs = 'S{ip1}_{jm1} + 2*z{im1}_{j}'.format(ip1=i + 1, jm1=j - 1, 
                                                     im1=i - 1, j=j)
            
            lhs = sympy.sympify(lhs)
            rhs = sympy.sympify(rhs)
            
            equations.append(sympy.Eq(lhs, rhs))
    
            # Now add in the edge cases
            if i == 1:
                edge_subs['z0_{j}'.format(j=j)] = 'S1_{jm1}'.format(jm1=j - 1)
            
            if j == 1:
                # We want i + 1 as the digit indices start at 1
                edge_subs['S{i}_0'.format(i=i)] = target_digits[i - 1]
            
            if i == k:
                edge_subs['S{kp1}_{jm1}'.format(kp1=k + 1, jm1=j - 1)] = target_digits[k + j - 1]
            
            if j == n - k:
                edge_subs['S{i}_{nmk}'.format(i=i, nmk=n - k)] = 0
                edge_subs['z{i}_{nmk}'.format(i=i, nmk=n - k)] = 0
            
            if i == k:
                edge_subs['z{k}_{j}'.format(k=k, j=j)] = 0
            
            if i == 1:
                edge_subs['S1_{nmkm1}'.format(nmkm1=n - k - 1)] = 0
    
    edge_subs = _sympify_dict(edge_subs)
    
    equations = [eqn.subs(edge_subs) for eqn in equations]
    
    return equations

if __name__ == "__main__":
    import doctest
    doctest.testmod()