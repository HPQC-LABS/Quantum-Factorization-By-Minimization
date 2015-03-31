# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 09:42:07 2015

@author: Richard
"""

import itertools


FACTOR_ROOTS = ('p', 'q')

def cmp_variables(x, y):
    ''' Given 2 variables, compare them according to root (first letter) and
        then index.

        >>> variables = 'p1, p2, p10, p11, p22, q1, q10, q9, z12'.split(', ')
        >>> sorted(variables, cmp=cmp_variables)
        ['p1', 'p2', 'p10', 'p11', 'p22', 'q1', 'q9', 'q10', 'z12']
    '''
    x, y = str(x), str(y)
    if x[0] != y[0]:
        return cmp(x[0], y[0])
    else:
        return cmp(int(x[1:]), int(y[1:]))

def cmp_variables_right(x, y):
    ''' Given 2 variables, compare them according to whether their root is a
        factor, and then according to the power of 2 the variable represents.
        Low powers of 2 are given priority

        >>> variables = 'p1, p2, p10, p11, p22, q1, q10, q9, z12'.split(', ')
        >>> sorted(variables, cmp=cmp_variables_right)
        ['p1', 'q1', 'p2', 'q9', 'p10', 'q10', 'p11', 'p22', 'z12']
    '''
    x, y = str(x), str(y)
    x_root = x[0]
    y_root = y[0]
    x_ind = int(x[1:])
    y_ind = int(y[1:])

    # Now force p and q to be the same
    if (x_root in FACTOR_ROOTS) and (y_root in FACTOR_ROOTS):
        x_root = y_root = 'p'

    if x_root != y_root:
        return cmp(x_root, y_root)
    else:
        return cmp(x_ind, y_ind)

def cmp_variables_left(x, y):
    ''' Given 2 variables, compare them according to whether their root is a
        factor, and then according to the power of 2 the variable represents.
        High powers of 2 are given priority

        >>> variables = 'p1, p2, p10, p11, p22, q1, q10, q9, z12'.split(', ')
        >>> sorted(variables, cmp=cmp_variables_left)
        ['p22', 'p11', 'p10', 'q10', 'q9', 'p2', 'p1', 'q1', 'z12']
    '''
    x, y = str(x), str(y)
    x_root = x[0]
    y_root = y[0]
    x_ind = int(x[1:])
    y_ind = int(y[1:])

    # Now force p and q to be the same
    if (x_root in FACTOR_ROOTS) and (y_root in FACTOR_ROOTS):
        x_root = y_root = 'p'

    if x_root != y_root:
        return cmp(x_root, y_root)
    else:
        return -1 * cmp(x_ind, y_ind)

def factorisation_summary(prod):
    ''' Print a summary of what's going on.
        >>> print factorisation_summary(143)
        Hamming distance: 2
        10001111 =
        1101 x
        1011

        >>> print factorisation_summary(536275525501)
        Hamming distance: 11
        111110011011100100000110001111101111101 =
        10000000000000110101 x
        11111001101100101001
    '''
    # Hack around the circular import
    from verification import get_target_factors

    summary = []
    p, q = get_target_factors(prod)
    summary.append('Hamming distance: {}'.format(factor_binary_differences(p, q)))
    summary.append('{} =\n{} x\n{}'.format(bin(prod)[2:], bin(p)[2:], bin(q)[2:]))
    return '\n'.join(summary)

def factor_binary_differences(p, q):
    ''' Given 2 factors, work out how many places they differ by when expressed
        in binary form

        >>> for m, n in itertools.combinations_with_replacement(range(1, 10), 2):
        ...     if len(bin(m)) != len(bin(n)):
        ...         continue
        ...     print m, n, factor_binary_differences(m, n)
        1 1 0
        2 2 0
        2 3 1
        3 3 0
        4 4 0
        4 5 1
        4 6 1
        4 7 2
        5 5 0
        5 6 2
        5 7 1
        6 6 0
        6 7 1
        7 7 0
        8 8 0
        8 9 1
        9 9 0

        >>> factor_binary_differences(524309, 534167)
        5
        >>> factor_binary_differences(1048573, 1048423)
        4
        >>> factor_binary_differences(1267650600228229401496703205361, 633825993891935921676532842551)
        78
    '''
    p_str = bin(p)[2:]
    q_str = bin(q)[2:]
    if len(p_str) != len(q_str):
        return None
    diffs = 0
    for pi, qi in itertools.izip(p_str, q_str):
        if pi != qi:
            diffs += 1
    return diffs

def num_to_factor_num_qubit(prod):
    ''' Given a number, work out how many qubits the factors should be.
        Return the largest first

        >>> num_to_factor_num_qubit(143)
        (4, 4)
        >>> num_to_factor_num_qubit(56153)
        (8, 8)
        >>> num_to_factor_num_qubit(1099551473989)
        (21, 21)
        >>> num_to_factor_num_qubit(309485009822787627980424653)
        (45, 45)
        >>> num_to_factor_num_qubit(1267650600228508624673600186743)
        (51, 51)
    '''
    bin_str = bin(prod)[2:]
    num_qub = len(bin_str)
    if num_qub % 2:
        return (num_qub + 1) / 2, (num_qub + 1) / 2
    else:
        return num_qub / 2, num_qub / 2


if __name__ == '__main__':
    import doctest
    doctest.testmod()