# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:16:25 2015

@author: Richard
"""

import itertools
import os

import numpy
import sympy

from carry_equations_generator import generate_carry_equations
from semiprime_tools import num_to_factor_num_qubit
from sympy_helper_fns import expressions_to_variables, is_constant, is_equation
from verification import check_solutions
from cfg_sympy_solver import EXPERIMENTS

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

def var_str(variables):
    ''' Return a list of variables '''
    str_var = map(str, variables)
    return ','.join(sorted(str_var))
    
def domain_str(variables):
    ''' Return the domain string of binary variables '''
    statements = ['{} in 0..1,'.format(v) for v in variables]
    return '\n    '.join(statements)

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
    PATH = '''C:\Users\Richard\Dropbox\Code\RSA\prolog'''
    filename = os.path.join(PATH, filename)
    with open(filename, 'w') as f:
        f.write(script_str)

if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    prod = 143
#    prod = 56153
#    prod = 4306239659
#    prod = EXPERIMENTS[100].product

    # 200x200    
#    prod = 645562469521727147413979793000752968582426448207305878207844816196118912415070998411753548531516944594279465790740639879

#    prod = 2945340432158418383223693624588738123559693482299075088767878449688292160397327779966295692450325070170031945807812908771881611572255401942922812303597161984656837427025571643491822766545899513818123285400690072855199653418293946527676736528882321172981935465021616458364488640088706931812060921466459680818892943

#    prod = 

#    fact1, fact2 = num_to_factor_num_qubit(prod)
#    equations = generate_carry_equations(fact1, fact2, prod)
#    equations = filter(is_equation, equations)


    from sympy_helper_fns import str_eqns_to_sympy_eqns
    equations = '''q2 + q7 + z8586 + z89 + 1 == 2*q2*q7 + q2*z8586 + 2*z910 + 4*z911
    q2*z8586 + 2*z910 + 4*z911 == z8586 + z89 + 2
    q7*z8586 == 0
    z8586*z910 == 0'''.split('\n')
    equations = str_eqns_to_sympy_eqns(equations)
    print equations

    write_script(get_script(equations), 'present2.pl')

#    variables = expressions_to_variables(equations)
#    print var_str(variables)
#    print domain_str(variables)
#    print equation_str(equations)
#    print print_str(variables)