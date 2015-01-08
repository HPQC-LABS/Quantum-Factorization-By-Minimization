# -*- coding: utf-8 -*-
"""
Created on Wed Jan 07 14:35:08 2015

@author: Richard
"""

from collections import defaultdict
import itertools

import sympy

#from sympy_helper_fns import


def _combine_terms((term1, term2)):
    term1, coef1 = term1
    term2, coef2 = term2

    atoms = term1.atoms().union(term2.atoms())
    return sympy.prod(atoms), coef1 * coef2

def _obj_func_term_dict_from_equations(equations):
    ''' Given a list of Sympy equations, calculate the term: coefficient dict
        of the objective function
    '''
    final_terms = defaultdict(int)
    for i, eqn in enumerate(equations):
#        print '{:.2f}% done'.format(100.0 * i / len(equations))
        terms = itertools.chain(eqn.lhs.as_coefficients_dict().iteritems(),
                                (-eqn.rhs).as_coefficients_dict().iteritems())
        products = itertools.product(terms, repeat=2)
        simplified = itertools.imap(_combine_terms, products)
        for t, c in simplified:
            final_terms[t] += c
    return final_terms

def objective_function_from_equations(equations):
    ''' Take a list of sympy equations and return an objective function

        >>> lhs = 'x'
        >>> rhs = '1'
        >>> eqns = [sympy.Eq(sympy.sympify(lhs), sympy.sympify(rhs))]
        >>> objective_function_from_equations(eqns)
        -x + 1

        >>> lhs = 'x + y'
        >>> rhs = 'x*y'
        >>> eqns = [sympy.Eq(sympy.sympify(lhs), sympy.sympify(rhs))]
        >>> objective_function_from_equations(eqns)
        -x*y + x + y
    '''
    # Watch out for edge cases
    if len(equations) == 0:
        return sympy.sympify('0')

    term_dict = _obj_func_term_dict_from_equations(equations)
    return sum([k*v for k, v in term_dict.iteritems()])

def _assign_atoms_to_index(atoms):
    ''' Take an iterable of atoms and return an {atom: integer id} map

        >>> atoms = sympy.sympify('x1*x2*y1*z*a').atoms()
        >>> _assign_atoms_to_index(atoms)
        {x2: 3, y1: 4, x1: 2, a: 1, z: 5}
    '''
    atoms = sorted(atoms, key=str)
    atom_map = {v : i + 1 for i, v in enumerate(atoms)}
    return atom_map

def _term_dict_to_coef_string(term_dict):
    ''' Given a dictionary of terms to coefficient, return the coefficient
        string for input into the optimisation code
    '''
    atoms = set.union(*[term.atoms(sympy.Symbol) for term in term_dict.iterkeys()])
    atom_map = _assign_atoms_to_index(atoms)

    lines = []
    for term, coef in term_dict.iteritems():
        var_num = sorted([atom_map[atom] for atom in term.atoms(sympy.Symbol)])
        line = ' '.join(map(str, var_num))
        if line:
            line += ' ' + str(coef)
        else:
            line = str(coef)
        lines.append(line)

    coef_str = '\n'.join(lines)
    atom_str = str(atom_map)
    return '\n\n'.join([coef_str, atom_str])

def expression_to_coef_string(expr):
    ''' Return the coefficient string for a sympy expression
        Also include the dictionary of variable number to original variable

        >>> inp = '40*s_1 + 30*s_1*s_2 + 100*s_1*s_2*s_3 - 15*s_2*s_3 - 20*s_3 + 4'
        >>> print expression_to_coef_string(inp)
        4
        1 2 30
        2 3 -15
        1 40
        3 -20
        1 2 3 100
        <BLANKLINE>
        {s_3: 3, s_2: 2, s_1: 1}
    '''
    if isinstance(expr, str):
        expr = sympy.sympify(expr)

    return _term_dict_to_coef_string(expr.as_coefficients_dict())

def equations_to_coef_string(equations):
    ''' Return the coefficient string of the objective function of some equations
        Also include the dictionary of variable number to original variable

        >>> lhs = 'x'
        >>> rhs = '1'
        >>> eqns = [sympy.Eq(sympy.sympify(lhs), sympy.sympify(rhs))]
        >>> print equations_to_coef_string(eqns)
        1
        1 -1
        <BLANKLINE>
        {x: 1}

        >>> lhs = 'x + y'
        >>> rhs = 'x*y'
        >>> eqns = [sympy.Eq(sympy.sympify(lhs), sympy.sympify(rhs))]
        >>> print equations_to_coef_string(eqns)
        1 2 -1
        1 1
        2 1
        <BLANKLINE>
        {x: 1, y: 2}
    '''
    term_dict = _obj_func_term_dict_from_equations(equations)
    return _term_dict_to_coef_string(term_dict)

if __name__ == "__main__":
    import doctest
    doctest.testmod()