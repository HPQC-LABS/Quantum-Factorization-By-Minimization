# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 09:42:07 2015

@author: Richard
"""

import itertools


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