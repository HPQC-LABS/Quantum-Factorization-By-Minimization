# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:16:25 2015

@author: Richard
"""

import itertools
import operator
import os

import numpy
import sympy

from carry_equations_generator import generate_carry_equations
from semiprime_tools import num_to_factor_num_qubit
from sympy_helper_fns import expressions_to_variables, is_constant, is_equation
from verification import check_solutions
from cfg_sympy_solver import EXPERIMENTS

TERMS_PER_LINE = 40

# Script from http://www.doc.ic.ac.uk/~sgc/teaching/pre2012/v231/lecture15.html

SCRIPT = '''?- use_module(library(clpfd)).

/* type go(A). */

go([{VAR_STR}]) :-

    /* Start the timer */
    statistics(runtime,_),

    /* Domains */
    {DOM_STR}


    /* Equation Constraints */
    {EQN_STR}

    /* Find the Answer */
    labeling([], [{VAR_STR}]),

    /* Write the Answer */
    {PRINT_STR}nl,

    /* Stop the timer */
    statistics(runtime,[_,TimeTaken]),
    write(TimeTaken),nl,nl.'''

def coef_matrix_to_prolog_script(filename='20x20_016_0_coef.txt', 
                                 directory=os.path.join('objective_functions')):
    ''' Given a coefficient matrix filepath, create a Prolog script that can
        be used to try and check the Hamiltonian it describes
    '''
    variables = set()
    terms = []
    with open(os.path.join(directory, filename), 'r') as f:
        for line in f.readlines():
            if (not line) or line.isspace() or line.startswith('{'):
                break
            data = line.strip().split(' ')
            coef = data.pop()
            vars_ = map(lambda x: 'X{}'.format(x), data)
            variables.update(vars_)
            term = '*'.join([coef] + vars_)
            terms.append('({})'.format(term))

    VAR_STR = var_str(variables)
    EQN_STR = hamiltonian_str(terms)
    DOM_STR = domain_str(variables)
    PRINT_STR = print_str(variables)
    script = SCRIPT.format(DOM_STR=DOM_STR, EQN_STR=EQN_STR, PRINT_STR=PRINT_STR,
                         VAR_STR=VAR_STR)
    
    with open(os.path.join(directory, filename.replace('.txt', '.pl')), 'w') as f:
        f.write(script)

            

def var_str(variables):
    ''' Return a list of variables '''
    str_var = map(str, variables)
    return ','.join(sorted(str_var))

def domain_str(variables):
    ''' Return the domain string of binary variables '''
    statements = ['{} in 0..1,'.format(v) for v in variables]
    return '\n    '.join(statements)

def hamiltonian_str(terms):
    ''' Given a list of terms, all of which are strings, add them together.
        If needed, break them up over several lines by adding new variables.    
    '''
#    return equation_str(['+'.join(terms) + ' #=< 0'])
    count = 0    
    equations, terms = terms_to_equations(terms, count)
    while len(terms) > TERMS_PER_LINE:
        count += len(terms)
        _equations, terms = terms_to_equations(terms, count)
        equations.extend(_equations)
    equations.append('+'.join(terms) + ' #=< 0')
    return equation_str(equations)
        

def terms_to_equations(terms, count):
    ''' Add terms to new lines and equate them to new line variables.
        Return the new lines, and a list of variables we've used to equate them    
    '''
    terms = itertools.izip_longest(*[iter(terms)] * TERMS_PER_LINE, 
                                    fillvalue=None)
    
    equations = []
    for l, line in enumerate(terms):
        line = '+'.join(filter(None, line))
        equations.append(line + ' #= L{}'.format(l + count))
        
    return equations, tuple('L{}'.format(_l + count) for _l in xrange(l + 1))

def equation_str(eqns):
    ''' Return a string representing the equations '''
    eqns = map(str, eqns)
    return '\n    '.join([eqn.replace('==', '#=').upper() + ',' for eqn in eqns])

def print_str(variables):
    ''' Return a string for the print statement '''
    variables = sorted(map(str, variables))
    writes = ['''write('{v} = '),write({v}),nl,'''.format(v=v) for v in variables]
    return '\n    '.join(writes)

def get_script(equations):
    ''' Return a script for prolog '''
    equations = filter(is_equation, equations)
    variables = expressions_to_variables(equations)
    variables = sorted([str(v).upper() for v in variables])
    VAR_STR = var_str(variables)
    EQN_STR = equation_str(equations)
    DOM_STR = domain_str(variables)
    PRINT_STR = print_str(variables)
    return SCRIPT.format(DOM_STR=DOM_STR, EQN_STR=EQN_STR, PRINT_STR=PRINT_STR,
                         VAR_STR=VAR_STR)

def write_script(script_str, filename='factor.pl'):
    ''' Write it down '''
    filename = os.path.join(os.path.pardir, 'prolog', filename)
    with open(filename, 'w') as f:
        f.write(script_str)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
