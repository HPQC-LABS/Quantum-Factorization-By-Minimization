# -*- coding: utf-8 -*-
"""
Created on Sat Apr 04 23:27:43 2015

@author: Richard
"""

import itertools
import re
import sympy
import timeit

from sympy.core.cache import clear_cache

from sympy_helper_fns import is_equation
from equivalence_dict import EquivalenceDict



### Plain sympy functions

def subs1(expr, to_sub):
    ''' Given an equation and dictionary of values to substitute, return the
        new equation. Added for completeness
    '''
    return subs1_many([expr], to_sub)[0]

def subs1_many(exprs, to_sub):
    ''' Substitute to_sub into many equations. Barebones wrapper to check we
        follow the original implementation.
    '''
    subbed = [expr.subs(to_sub, simultaneous=True) for expr in exprs]
    clear_cache()
    return subbed

### Our better sympy subs

def subs2(expr, to_sub):
    ''' Our own proper, recursive subs function '''
    if isinstance(expr, sympy.Integer):
        return expr

    elif isinstance(expr, sympy.Symbol):
        return to_sub.get(expr, expr)

    elif is_equation(expr, check_true=False):
        return sympy.Eq(subs2(expr.lhs, to_sub), subs2(expr.rhs, to_sub))
    
    elif isinstance(expr, sympy.Add):
        return sum([subs2(arg, to_sub) for arg in expr.args])

    elif isinstance(expr, sympy.Mul):
        out = 1
        for arg in expr.args:
            out *= subs2(arg, to_sub)
        return out

    elif isinstance(expr, sympy.Pow):
        base_, exp_ = expr.args
        return subs2(base_, to_sub) ** subs2(exp_, to_sub)
    
    else:
        raise ValueError('Unknown type of {}: {}'.format(expr, type(expr)))


def subs2_many(exprs, to_sub):
    ''' Wrapper for singular, since there is no benefit to considering them all
    '''
    return [subs2(expr, to_sub) for expr in exprs]

### Our nasty string subs

def subs3(expr, to_sub):
    ''' Use horrible string mangling to substitute variables
    '''
    return subs3_many([expr], to_sub)[0]

def subs3_many(exprs, to_sub):
    ''' Substitute using regular expressions and sympify
    '''
    if not exprs:
        return []

    exprs = map(str, exprs)
    exprs = ', '.join(exprs)

    allowed_neighbours = '(\D|\Z)'
    for k, v in to_sub.iteritems():
        pattern = str(k) + allowed_neighbours
        repl = '(' + str(v) + ')\\1'
        exprs = re.sub(pattern, repl, exprs)

    exprs = exprs.split(', ')

    out = []
    for eqn in exprs:
        if '==' in eqn:
            eqn = sympy.Eq(*map(sympy.sympify, eqn.split('==')))
        else:
            eqn = sympy.sympify(eqn)
        if isinstance(eqn, bool):
            if eqn:
                out.append(sympy.boolalg.BooleanTrue())
            else:
                out.append(sympy.boolalg.BooleanFalse())
        else:
            out.append(eqn)
    return out


### Main functions

## This determines which function subs and subs_many point to, and hence most
## of our work
DEFAULT_SUBS_MANY = subs2_many

def subs(expr, to_sub):
    ''' Function that most modules will call for substituting a to_sub into
        a single equation. Just point to our main subs_many function
    '''
    return subs_many([expr], to_sub)[0]

def subs_many(exprs, to_sub):
    ''' Function most scripts will call for substituting to_sub into multiple
        equations or expressions
    '''
    # First dictify our to_sub so it plays nicely with sympy
    if not isinstance(to_sub, dict):
        to_sub = dict(to_sub)
    
    unitary_subs = {}
    compound_subs = {}
    for k, v in to_sub.iteritems():
        if len(k.atoms()) == 1:
            unitary_subs[k] = v
        else:
            compound_subs[k] = v

    subbed = DEFAULT_SUBS_MANY(exprs, unitary_subs)
    # Revert back to sympy subs for anything complicated
    subbed = subs1_many(subbed, compound_subs)
    return subbed


### Debug and testing

def _are_equal(expr1, expr2):
    ''' Given 2 expressions, work out whether they are equal. Used only in
        tests below
    '''
    # Here we can use is as Python True is singular, as is sympy's True
    if expr1 is expr2:
        assert str(expr1) == str(expr2)
        return True
    # For tests, we don't care if the evaluation is true or not
    if is_equation(expr1, check_true=False):
        if not is_equation(expr2, check_true=False):
            print '{} != {}'.format(expr1, expr2)
            return False
        
        return _are_equal(expr1.lhs, expr2.lhs) and _are_equal(expr1.rhs, expr2.rhs)

    diff = (expr1 - expr2).expand()
    if diff == sympy.S.Zero:
        if str(expr1) != str(expr2):
            print 'Double check:\t{} == {}'.format(expr1, expr2)
        return True
    else:
        print '{} != {}'.format(expr1, expr2)
        return False

def _profile(sub_func=subs):
    ''' Profile a function against sympy's subs '''
    num_var = 6
    var = sympy.symbols(' '.join(['x{}'.format(i) for i in xrange(num_var)]))
    terms = itertools.product(var, repeat=2)
    expr = sum([a*b for a, b in terms])
    vals = itertools.product(range(2), repeat=num_var)
    for val in vals:
        to_sub = dict(zip(var, val))
        _expr = sub_func(expr, to_sub)
    return
    

if __name__ == "__main__":
    import doctest
    doctest.testmod()

    symbols = sympy.symbols('x1 x2 x3 x4')
    x1, x2, x3, x4 = symbols
    eqns = [sympy.Eq(x1, x2),
            sympy.Eq(x1*x2, 0),
            sympy.Eq(x1 + x2*x3),
            sympy.Eq(x1**2, x2 - 3*x3),
            x1*x2*x3 - 4*x4**3,
            sympy.Eq(x1 + x2 - 2*x3, x4)]
    to_subs = [{x1: 1},
               {x1: x2},
               {x2: x3, x4: x2},
               {x1: x2, x2: x3, x4: 0, x3: x1},
               {x1: x2, x2: 0, x4: 1},
               {x1: x2 + x4, x2: 2, x4: 1},
               {x1: 1 - x2, x2: -82, x4: 1},
               {x1*x2: 0, x2*x3: 1},
            ]
    
    for to_sub in to_subs:

        # Work it out the proper way
        sympy_sol = [eqn.subs(to_sub, simultaneous=True) for eqn in eqns]
        # Work it out with whatever our singular function is
        subs_sol = [subs(eqn, to_sub) for eqn in eqns]
        # Work it out with whatever our batch function is
        subs_many_sol = subs_many(eqns, to_sub)
        
        # Check we haven't done anything really crazy
        assert len(sympy_sol) == len(subs_sol) == len(subs_many_sol)
        
        # Now check they're all equal
        for orig, target, ssol, smsol in zip(eqns, sympy_sol, subs_sol, subs_many_sol):
            # Check we're doing what sympy is
            assert _are_equal(target, ssol)
            # Check we're doing what we think we're doing!
            assert _are_equal(ssol, smsol)


    ### Profile the subs methods
    setup_str = 'from __main__ import subs, subs1, subs2, subs3, _profile'
    num_trial = 10
    time0 = timeit.timeit("_profile(subs)", setup_str, number=num_trial)
    print 'subs: {:.2f}s'.format(time0)
    time1 = timeit.timeit("_profile(subs1)", setup_str, number=num_trial)
    print 'subs1: {:.2f}s'.format(time1)
    time2 = timeit.timeit("_profile(subs2)", setup_str, number=num_trial)
    print 'subs2: {:.2f}s'.format(time2)
    time3 = timeit.timeit("_profile(subs3)", setup_str, number=num_trial)
    print 'subs3: {:.2f}s'.format(time3)
