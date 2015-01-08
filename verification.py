# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 18:35:14 2014

@author: Richard
"""

import itertools
import sympy

from sympy_helper_fns import is_constant
from rsa_constants import RSA100, RSA100_F1, RSA100_F2

KNOWN_FACTORISATIONS = {
        RSA100: (RSA100_F1, RSA100_F2)

}
BRUTE_FORCE_FACTORISATION_LIMIT = 10**16

def _extract_solutions(solutions, max_digits):
    ''' Solutions is a dict of solutions we want to look in. num_digits is
        the total number of digits in a factor.
        Returns the solution according to solutions.
        NOTE only works for factors with same number of digits
    '''
    p = [1]
    q = [1]
    for i in xrange(max_digits - 2, 0, -1):
        p.append(solutions.get(sympy.Symbol('p{}'.format(i))))
        q.append(solutions.get(sympy.Symbol('q{}'.format(i))))
    p.append(1)
    q.append(1)
    return p, q

def get_target_factors(product):
    ''' Return tuple of target factors, or None if the can't be found '''

    known_factors = KNOWN_FACTORISATIONS.get(product)
    if known_factors is not None:
        return known_factors

    # Work out our target ps and qs
    if product < BRUTE_FORCE_FACTORISATION_LIMIT:
        return factorise(product)

    return None


def check_solutions(product, solutions, verbose=False):
    ''' Check that solutions are consistent with the binary factorisation.
        NOTE Only works with prime numbers with the same number of digits in
        their factors

        >>> p1, p2, p3, p10, p160, q1, q3, q164 = sympy.var('p1 p2 p3 p10 p160 q1 q3 q164')
        >>> soln = {p1: 0, p2: 1, p3: 1, p10: 1, p160: 1, q1: 1, q3: 0, q164: 1}
        >>> check_solutions(RSA100, soln, verbose=True)
        All assertions passed.
    '''

    # Get the target factors, turn them into binary and do some checks
    target_factors = get_target_factors(product)
    if target_factors is None:
        return
    target_factors = map(bin, target_factors)

    # Check we have 2 factors
    assert len(target_factors) == 2

    # Check the same length
    assert len(target_factors[0]) == len(target_factors[1])

    # Trim off the first '0b' and check we have a 1 at each end
    target_factors = [fact[2:] for fact in target_factors]
    target_factors = [map(int, fact) for fact in target_factors]
    for fact in target_factors:
        assert fact[0] == fact[-1] == 1
    target_p, target_q = target_factors

    # Now extract a similar list from the given solutions dict
    digits_in_multiplicand = len(target_p)
    soln_p, soln_q = _extract_solutions(solutions, digits_in_multiplicand)

    symbolic_pairs = []
    # Check fully determined ps and qs and extract symbolic matchings
    for sol, target in itertools.chain(itertools.izip(soln_p, target_p),
                                       itertools.izip(soln_q, target_q)):
        # No solutions found whatsoever
        if sol is None:
            continue
        # If constant, we can check that easily
        if is_constant(sol):
            assert sol == target
        else:
            symbolic_pairs.append((sol, target))

    # For symbolic matchings, just check that no symbolic expression is mapped to
    # two different values
    #TODO Do something cleverer here, like plug into a new EquationSolver.
    symbolic_map = {}
    for sol, tar in symbolic_pairs:
        prev_tar = symbolic_map.get(sol)
        if prev_tar is None:
            symbolic_map[sol] = tar
        else:
            assert tar == prev_tar

    if verbose:
        print 'All assertions passed.'
        if symbolic_map:
            print 'Check the below are consistent'
            print symbolic_map


def factorise(n):
    ''' Return list of factors '''
    factors = []
    i = 2
    while True:
        if n == 1:
            break
        while not n % i:
            factors.append(i)
            n /= i
        i += 1
    factors.reverse()
    return factors

def binary_factorisation(n):
    ''' Return the binary factorisation of a number '''
    factors = factorise(n)
    bin_fact = map(bin, factors)
    return bin_fact


def print_binary_factorisation(n):
    ''' Print the binary factorisation of a number

        >>> print_binary_factorisation(143)
        143
        =0b10001111
        Factors:
          1101
          1011
    '''
    bin_fact = binary_factorisation(n)
    max_len = len(max(bin_fact, key=len))
    fact_str = '\n'.join([b_f[2:].rjust(max_len) for b_f in bin_fact])
    print '{}\n={}\nFactors:\n{}'.format(n, bin(n), fact_str)


if __name__ == '__main__':
    import doctest
    doctest.testmod()