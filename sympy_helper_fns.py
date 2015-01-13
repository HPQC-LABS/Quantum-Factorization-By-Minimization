"""
Created on Fri Dec 26 12:35:16 2014

Helper functions to deal with sympy expressions and equations

@author: Richard Tanburn
"""

import fractions
import re
import sympy

from contradiction_exception import ContradictionException

## Helper functions

def is_one_or_zero(val):
    ''' Self explanatory '''
    return (val == 0) or (val == 1)

def num_add_terms(expr, check=False):
    ''' Return the number of additive terms in an expression.
        Note doesn't work for multiplicative terms!

        >>> expr = 'x + 2*y'
        >>> num_add_terms(sympy.sympify(expr))
        2
        >>> expr = 'x + 2*y + 3*z + 5'
        >>> num_add_terms(sympy.sympify(expr))
        4
        >>> expr = 'x * 2*y'
        >>> num_add_terms(sympy.sympify(expr))
        1
        >>> expr = '(x + 2*y) * z'
        >>> num_add_terms(sympy.sympify(expr), check=False)
        1
    '''
    if check:
        assert expr.func == sympy.Add
    return len(expr.as_ordered_terms())

def is_constant(expr):
    ''' Determine whether an expression is constant
        >>> expr = 'x + 2*y'
        >>> is_constant(sympy.sympify(expr))
        False
        >>> expr = 'x + 5'
        >>> is_constant(sympy.sympify(expr))
        False
        >>> expr = '3'
        >>> is_constant(sympy.sympify(expr))
        True
        >>> expr = '2*x - 4'
        >>> is_constant(sympy.sympify(expr))
        False
    '''
    if isinstance(expr, int):
        return True
    return len(expr.atoms(sympy.Symbol)) == 0

def is_monic(expr):
    ''' Determine whether an expression is monic
        >>> expr = 'x + 2*y'
        >>> is_monic(sympy.sympify(expr))
        False
        >>> expr = 'x + 5'
        >>> is_monic(sympy.sympify(expr))
        False
        >>> expr = '3'
        >>> is_monic(sympy.sympify(expr))
        False
        >>> expr = '2*x - 4'
        >>> is_monic(sympy.sympify(expr))
        False
        >>> expr = 'x'
        >>> is_monic(sympy.sympify(expr))
        True
        >>> expr = 'x*y + z'
        >>> is_monic(sympy.sympify(expr))
        True
        >>> expr = 'x*y + z**2 + 1'
        >>> is_monic(sympy.sympify(expr))
        True
    '''
    return all(coef == 1 for coef in expr.as_coefficients_dict().itervalues())

def is_equation(eqn):
    ''' Return True if it is an equation rather than a boolean value.
        If it is False, raise a ContradictionException. We never want anything
        that might be False

        >>> x, y = sympy.symbols('x y')
        >>> eq1 = sympy.Eq(x, y)
        >>> eq2 = sympy.Eq(x, x)
        >>> eq3 = sympy.Eq(x, y).subs(y, x)

        >>> is_equation(eq1)
        True
        >>> is_equation(eq2)
        False
        >>> is_equation(eq3)
        False
        
        Now check that it raises exceptions for the right things
        >>> is_equation(0)
        False
        >>> is_equation(False)
        Traceback (most recent call last):
            ...
        ContradictionException: False equation
    '''
    if isinstance(eqn, sympy.boolalg.BooleanFalse) or (eqn is False):
        raise ContradictionException('False equation')
    return isinstance(eqn, sympy.Equality)

def parity(expr):
    ''' Return parity:
        0 - even
        1 - odd
        None - undetermined

        >>> expr = 'x + 2*y'
        >>> parity(sympy.sympify(expr))

        >>> expr = '2*x + 2*y'
        >>> parity(sympy.sympify(expr))
        0
        >>> expr = 'x + 5'
        >>> parity(sympy.sympify(expr))

        >>> expr = '3'
        >>> parity(sympy.sympify(expr))
        1

        >>> expr = '2*x - 4'
        >>> parity(sympy.sympify(expr))
        0

        >>> expr = '2*x + 4*x*y + 1'
        >>> parity(sympy.sympify(expr))
        1
    '''
    parity = 0
    for term, coef in expr.as_coefficients_dict().iteritems():
        if is_constant(term):
            parity += coef % 2
        else:
            if coef % 2:
                return None
    return parity


def cancel_constant_factor(eqn, maintain_sign=False):
    ''' Divide the equation by the hcf of all the terms.
        If every term is negative, then also divide by -1.
        If maintain_sign is True, then the dividing factor cannot be negative

        >>> lhs = sympy.sympify('2*x + 2')
        >>> rhs = sympy.sympify('2*y + 3')
        >>> cancel_constant_factor(sympy.Eq(lhs, rhs))
        2*x + 2 == 2*y + 3

        >>> lhs = sympy.sympify('2*x - 2')
        >>> rhs = sympy.sympify('4*y + 6')
        >>> cancel_constant_factor(sympy.Eq(lhs, rhs))
        x - 1 == 2*y + 3

        >>> lhs = sympy.sympify('15*x + 3')
        >>> rhs = sympy.sympify('45*y')
        >>> cancel_constant_factor(sympy.Eq(lhs, rhs))
        5*x + 1 == 15*y

        >>> lhs = sympy.sympify('15*x + 3')
        >>> cancel_constant_factor(sympy.Eq(lhs))
        5*x + 1 == 0

        Negative equations
        >>> lhs = sympy.sympify('-3*x - 6')
        >>> rhs = sympy.sympify('-9*y')
        >>> cancel_constant_factor(sympy.Eq(lhs, rhs))
        x + 2 == 3*y
        >>> cancel_constant_factor(sympy.Eq(lhs, rhs), maintain_sign=True)
        -x - 2 == -3*y
    '''
    if not is_equation(eqn):
        return eqn

    coef = (eqn.lhs.as_coefficients_dict().values() +
            eqn.rhs.as_coefficients_dict().values())
    hcf = reduce(fractions.gcd, coef)
    if maintain_sign:
        hcf = abs(hcf)
    return sympy.Eq(eqn.lhs / hcf, eqn.rhs / hcf)


def balance_constant(eqn):
    ''' Take an equation and tidy up the constant part

        >>> lhs = sympy.sympify('x + 2')
        >>> rhs = sympy.sympify('y + 3')
        >>> balance_constant(sympy.Eq(lhs, rhs))
        x == y + 1

        >>> lhs = sympy.sympify('x - 2')
        >>> rhs = sympy.sympify('y + 3')
        >>> balance_constant(sympy.Eq(lhs, rhs))
        x == y + 5

        >>> lhs = sympy.sympify('x + 2')
        >>> rhs = sympy.sympify('y - 3')
        >>> balance_constant(sympy.Eq(lhs, rhs))
        x + 5 == y
    '''
    if not is_equation(eqn):
        return eqn

    lhs_c = eqn.lhs.as_coeff_add()[0]
    rhs_c = eqn.rhs.as_coeff_add()[0]

    if (lhs_c < 0) or (rhs_c < 0):
        to_add = abs(min(lhs_c, rhs_c))
    else:
        to_add = - min(lhs_c, rhs_c)

    return sympy.Eq(eqn.lhs + to_add,
                    eqn.rhs + to_add)

def balance_terms(eqn):
    ''' Take an equation and make sure all terms are positive. Also, if LHS is
        constant, flip.

        >>> lhs = sympy.sympify('x + 2')
        >>> rhs = sympy.sympify('- y + 3')
        >>> balance_terms(sympy.Eq(lhs, rhs))
        x + y == 1

        >>> lhs = sympy.sympify('x - 2 - x*y')
        >>> rhs = sympy.sympify('y + 3 + 10*x*y')
        >>> balance_terms(sympy.Eq(lhs, rhs))
        x == 11*x*y + y + 5
    '''
    if not is_equation(eqn):
        return eqn

    lhs_terms = eqn.lhs.as_coefficients_dict()
    rhs_terms = eqn.rhs.as_coefficients_dict()

    terms = set(lhs_terms.keys()).union(set(rhs_terms.keys()))

    for term in terms:
        lhs_c = lhs_terms.get(term)
        rhs_c = rhs_terms.get(term)

        if lhs_c is None:
            lhs_c = 0
        if rhs_c is None:
            rhs_c = 0

        if (lhs_c < 0) or (rhs_c < 0):
            to_add = abs(min(lhs_c, rhs_c))
        else:
            to_add = - min(lhs_c, rhs_c)

        eqn = sympy.Eq(eqn.lhs + to_add * term,
                       eqn.rhs + to_add * term)

    if is_constant(eqn.lhs):
        eqn = sympy.Eq(eqn.rhs, eqn.lhs)

    return eqn

def max_value(expr):
    ''' Return the max value of an expression

        >>> expr = sympy.sympify('x+1')
        >>> max_value(expr)
        2

        >>> expr = sympy.sympify('x-1')
        >>> max_value(expr)
        0

        >>> expr = sympy.sympify('x + y + 1')
        >>> max_value(expr)
        3

        >>> expr = sympy.sympify('2*x - y + 1')
        >>> max_value(expr)
        3

        >>> expr = sympy.sympify('x*y*z + a*b')
        >>> max_value(expr)
        2

        >>> expr = sympy.sympify('-2 * (x+1)')
        >>> max_value(expr)
        -2
    '''
    #expr = sympy.simplify(expr)
    coef = expr.as_coefficients_dict()
    max_ = 0
    for term, c in coef.iteritems():
        if c > 0 or term == 1:
            max_ += c
    return max_

def min_value(expr):
    ''' Return the min value of an expression

        >>> expr = sympy.sympify('x+1')
        >>> min_value(expr)
        1

        >>> expr = sympy.sympify('x-1')
        >>> min_value(expr)
        -1

        >>> expr = sympy.sympify('x + y +1')
        >>> min_value(expr)
        1

        >>> expr = sympy.sympify('2*x - y + 1')
        >>> min_value(expr)
        0

        >>> expr = sympy.sympify('x*y*z + a*b - 4*k + 2')
        >>> min_value(expr)
        -2

        >>> expr = sympy.sympify('-2 * (x+1)')
        >>> min_value(expr)
        -4
    '''
    #expr = sympy.simplify(expr)
    coef = expr.as_coefficients_dict()
    min_ = 0
    for term, c in coef.iteritems():
        if c < 0 or term == 1:
            min_ += c
    return min_

def remove_binary_squares_eqn(eqn):
    ''' Given an equation, remove all of the squares as any binary
        variable squared is itself.

        >>> lhs = sympy.sympify('x**2 * y + z**3 + 2*z - 4')
        >>> rhs = sympy.sympify('x + 1')
        >>> remove_binary_squares_eqn(sympy.Eq(lhs, rhs))
        x*y + 3*z - 4 == x + 1

        >>> expr1 = sympy.sympify('2*x**2 + 1')
        >>> expr2 = sympy.sympify('y**4 + z')
        >>> remove_binary_squares_eqn(sympy.Eq(expr1, expr2))
        2*x + 1 == y + z
    '''
    if not is_equation(eqn):
        return eqn
    #eqn = eqn.subs({var ** 2: var for var in eqn.atoms(sympy.Symbol)})

    return sympy.Eq(remove_binary_squares(eqn.lhs),
                    remove_binary_squares(eqn.rhs))

def remove_binary_squares(expr):
    ''' Given an equation, remove all of the squares as any binary
        variable squared is itself.

        >>> expr = 'x**2 * y + z**3 + 2*z - 4'
        >>> remove_binary_squares(sympy.sympify(expr))
        x*y + 3*z - 4

        >>> expr = '(x*y)**2 + z**3 + 1'
        >>> remove_binary_squares(sympy.sympify(expr))
        x*y + z + 1

        Because of the new implementation, we want to check the variables
        are exactly equivalent
        >>> x = sympy.symbols('x')
        >>> remove_binary_squares(x ** 3) == x
        True
        >>> remove_binary_squares(x ** 3) is x
        True
    '''
    exp_match = re.compile('[a-zA-Z][0-9]*\*\*[0-9]*')
    matches = re.findall(exp_match, str(expr))
    for match in matches:
        var, exp = match.split('**')
        var = sympy.sympify(var)
        exp = int(exp)
        expr = expr.subs(var ** exp, var)
    return expr

    # Old, memory hungry implementation
#    w = sympy.Wild('w')
#    p = sympy.Wild('p')
#    expr = expr.replace(w ** p, w, exact=True)
#    return expr

def expressions_to_variables(exprs):
    ''' Take a list of equations and return a set of variables 

        >>> eqn = sympy.Eq(sympy.sympify('x*a + 1'))
        >>> to_test = [sympy.sympify('x + y*z + 2*a^b'), eqn]
        >>> expressions_to_variables(to_test)
        set([x, z, a, b, y])
    '''
    assert all(map(lambda x: isinstance(x, sympy.Basic), exprs))
    return set.union(*[expr.atoms(sympy.Symbol) for expr in exprs])

if __name__ == "__main__":
    import doctest
    doctest.testmod()