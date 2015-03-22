# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 14:53:21 2015

@author: Richard
"""

#!/usr/bin/env python


import sys
from time import time

from cfg_sympy_solver import EXPERIMENTS, QUBIT_REDUCTION_ID
from objective_function_helper import coef_str_to_file
from sympy_assumptions import (make_simultaneous_assumptions, 
                               frequency_rank_variables,
                               weighted_frequency_rank_variables,
                               max_coef_rank_variables,
                               lexographical_rank_variable)
from sympy_solver import EquationSolver
from verification import check_solutions


__author__ = "Nathaniel Bryans"
__credits__ = ["Nathaniel Bryans", "Nikesh Dattani"]
__version__ = "0.0.6"
__status__ = "Prototype"

# Default file name
OutputFileName = "output.txt"

# A default experiment to run
exp = 100
params = EXPERIMENTS[exp]
digitsInMultiplicand1, digitsInMultiplicand2, product = params

# Default assumption parameters
num_assumptions = 0
limit_assumptions = 1
# Use standard carry equations by default
qubit_reduction_method = 0

equation_generator, coef_str_generator = QUBIT_REDUCTION_ID[qubit_reduction_method]

eqns = equation_generator(digitsInMultiplicand1, digitsInMultiplicand2, product)

s = time()
# None means the result will be printed to screen
output = None#OutputFileName


# We can use the handy state caching    
cache_name = '_state_rsa100'#'_state_{}'.format(str(product)[-6:])
log_deductions = False
if cache_name is not None:
    try:
        raise Exception()
        system = EquationSolver.from_disk(cache_name)
    except Exception as e:
        print e
        system = EquationSolver(eqns, output_filename=output, 
                                            log_deductions=log_deductions)
        system.solve_equations(verbose=True)
        system.to_disk(cache_name)

# Do it normally
else:
    system = EquationSolver(eqns, output_filename=output, 
                                        log_deductions=log_deductions)
    system.solve_equations(verbose=True)
    

system.print_summary()
print '\nSolved in {}s'.format(time() - s)

check_solutions(product, system.solutions.copy(), verbose=True)