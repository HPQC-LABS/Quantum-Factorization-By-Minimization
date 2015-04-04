# -*- coding: utf-8 -*-
"""
Created on Sat Apr 04 23:27:43 2015

@author: Richard
"""

import re
import sympy
import timeit

from sympy_helper_fns import is_equation

### Plain sympy functions

def subs1(expr, to_sub):
    ''' Given an equation and dictionary of values to substitute, return the
        new equation. Added for completeness
    '''
    return subs1_many([expr], to_sub)[0]

def subs1_many(exprs, to_sub):
    ''' Substitute to_sub into many equations. Barebones wrapper to check we
        follow the original implementation.
    
        >>> symbols = sympy.symbols('x1 x2 x3 x4')
        >>> x1, x2, x3, x4 = symbols
        >>> eqns = [sympy.Eq(x1, x2),
        ...         sympy.Eq(x1*x2, 0),
        ...         sympy.Eq(x1 + x2*x3),
        ...         sympy.Eq(x1**2, x2 - 3*x3),
        ...         x1*x2*x3 - 4*x4**3]
        >>> to_subs = [{x1: 1},
        ...             {x1: x2},
        ...             {x2: x3, x4: x2},
        ...             {x1: x2, x2: x3, x4: 0},
        ...             {x1: x2, x2: x3, x4: 1, x3: x1}]
        
        >>> for to_sub in to_subs:
        ...     for new_eqn in subs1_many(eqns, to_sub):
        ...         print new_eqn
        1 == x2
        x2 == 0
        x2*x3 + 1 == 0
        1 == x2 - 3*x3
        x2*x3 - 4*x4**3
        True
        x2**2 == 0
        x2*x3 + x2 == 0
        x2**2 == x2 - 3*x3
        x2**2*x3 - 4*x4**3
        x1 == x3
        x1*x3 == 0
        x1 + x3**2 == 0
        x1**2 == -2*x3
        x1*x3**2 - 4*x2**3
        True
        x3**2 == 0
        x3**2 + x3 == 0
        x3**2 == -2*x3
        x3**3
        True
        x1**2 == 0
        x1**2 + x1 == 0
        x1**2 == -2*x1
        x1**3 - 4
    '''
    return [expr.subs(to_sub) for expr in exprs]

### Our better sympy subs

def subs2(expr, to_sub):
    ''' Our own proper, recursive subs function '''
    pass

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


def subs(expr, to_sub):
    ''' Function that most modules will call for substituting a to_sub into
        a single equation. Just point to our main subs_many function
    '''
    return subs_many([expr], to_sub)[0]

def subs_many(exprs, to_sub):
    ''' Function most scripts will call for substituting to_sub into multiple
        equations or expressions

        >>> symbols = sympy.symbols('x1 x2 x3 x4')
        >>> x1, x2, x3, x4 = symbols
        >>> eqns = [sympy.Eq(x1, x2),
        ...         sympy.Eq(x1*x2, 0),
        ...         sympy.Eq(x1 + x2*x3),
        ...         sympy.Eq(x1**2, x2 - 3*x3),
        ...         x1*x2*x3 - 4*x4**3]
        >>> to_subs = [{x1: 1},
        ...             {x1: x2},
        ...             {x2: x3, x4: x2},
        ...             {x1: x2, x2: x3, x4: 0},
        ...             {x1: x2, x2: 2, x4: 1}]
        
        >>> for to_sub in to_subs:
        ...     for new_eqn in subs_many(eqns, to_sub):
        ...         print new_eqn
        1 == x2
        x2 == 0
        x2*x3 + 1 == 0
        1 == x2 - 3*x3
        x2*x3 - 4*x4**3
        True
        x2**2 == 0
        x2*x3 + x2 == 0
        x2**2 == x2 - 3*x3
        x2**2*x3 - 4*x4**3
        x1 == x3
        x1*x3 == 0
        x1 + x3**2 == 0
        x1**2 == -2*x3
        x1*x3**2 - 4*x3**3
        True
        x3**2 == 0
        x3**2 + x3 == 0
        x3**2 == -2*x3
        x3**3
        True
        False
        2*x3 + 2 == 0
        4 == -3*x3 + 2
        4*x3 - 4
    '''
    return subs1_many(exprs, to_sub)

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
            return False
        
        diff1 = expr1.lhs - expr2.lhs
        diff2 = expr1.rhs - expr2.rhs
        if diff1 == diff2 == sympy.S.Zero:
            assert str(diff1) == str(diff2)
            return True
        return False

    if expr1 - expr2 == sympy.S.Zero:
        assert str(expr1) == str(expr2)
        return True
    else:
        return False

def _profile(func=subs_many):
    ''' Profile a function against sympy's subs '''
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
    symbols = sympy.symbols('x1 x2 x3 x4')
    x1, x2, x3, x4 = symbols
    eqns = [sympy.Eq(x1, x2),
            sympy.Eq(x1*x2, 0),
            sympy.Eq(x1 + x2*x3),
            sympy.Eq(x1**2, x2 - 3*x3),
            x1*x2*x3 - 4*x4**3]
    to_subs = [{x1: 1},
               {x1: x2},
               {x2: x3, x4: x2},
               {x1: x2, x2: x3, x4: 0, x3: x1},
               {x1: x2, x2: 2, x4: 1},
               {x1: x2 + x4, x2: 2, x4: 1},
               {x1: 1 - x2, x2: -82, x4: 1},
            ]
    
    for to_sub in to_subs:
        # Work it out the proper way
        sympy_sol = [eqn.subs(to_sub) for eqn in eqns]
        # Work it out with whatever our singular function is
        subs_sol = [subs(eqn, to_sub) for eqn in eqns]
        # Work it out with whatever our batch function is
        subs_many_sol = subs_many(eqns, to_sub)
        
        # Check we haven't done anything really crazy
        assert len(sympy_sol) == len(subs_sol) == len(subs_many_sol)
        
        # Now check they're all equal
        for target, ssol, smsol in zip(sympy_sol, subs_sol, subs_many_sol):
            # Check we're doing what sympy is
            _are_equal(target, ssol)
            # Check we're doing what we think we're doing!
            _are_equal(ssol, smsol)
