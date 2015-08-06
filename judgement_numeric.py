"""
Created on Wed Jul 22 17:10:16 2015

@author: richard
"""

import math

import sympy

from semiprime_tools import num_to_factor_num_qubit
from verification import get_target_factors

def get_pq_interval(product):
    ''' Given that pq = product, and p and q have equal given lengths, work
        out a interval they must fall in.
        Suppose also that p < q
        
        >>> get_pq_interval(143)
        ((9, 11), (11, 15))

        >>> product = EXPERIMENTS[2].product
        >>> get_pq_interval(product)
        ((220, 236), (236, 255))

        >>> product = EXPERIMENTS[10].product
        >>> get_pq_interval(product)
        ((65536, 65558), (65558, 65582L))
    '''
    p_len, q_len = num_to_factor_num_qubit(product)

    root_product = int(math.sqrt(product))

    lower_bound = 2 ** (p_len - 1)
    upper_bound = 2 ** q_len - 1

    factors = get_target_factors(product)
    p, q = min(factors), max(factors)

    p_interval = max(lower_bound, product / upper_bound), root_product
    q_interval = root_product, min(upper_bound, product / lower_bound)

    assert lower_bound <= p <= q <= upper_bound
    assert p_interval[0] <= p <= p_interval[1]
    assert q_interval[0] <= q <= q_interval[1]

    return p_interval, q_interval

def get_pq_deductions(product):
    ''' Given a product which is the product of p < q, deduce anything about
        the interval in which p and q lie and turn it into deductions

        >>> get_pq_deductions(143)
        {p2: 0}

        >>> product = EXPERIMENTS[2].product
        >>> get_pq_deductions(product)
        {q6: 1, p6: 1, q5: 1}

        >>> product = EXPERIMENTS[10].product
        >>> get_pq_deductions(product)
        {p5: 0, q14: 0, q10: 0, p9: 0, p10: 0, q6: 0, p14: 0, q13: 0, p13: 0, p7: 0, p8: 0, q12: 0, q8: 0, q9: 0, p12: 0, p6: 0, p11: 0, q7: 0, q11: 0, q15: 0, p15: 0}
    '''
    deductions = {}
    for stem, interval in zip(('p', 'q'), get_pq_interval(product)):
        lower, upper = [bin(num)[2:] for num in interval]

        assert len(lower) == len(upper)

        for i, (lo, hi) in enumerate(zip(lower, upper)):
            if i == 0:
                continue
            if lo == hi:
                symbol = stem + str(len(lower) - i - 1)
                symbol = sympy.Symbol(symbol)
                deductions[symbol] = int(lo)
            else:
                break
    return deductions


if __name__ == '__main__':
    from cfg_sympy_solver import EXPERIMENTS
    import doctest
    doctest.testmod()
