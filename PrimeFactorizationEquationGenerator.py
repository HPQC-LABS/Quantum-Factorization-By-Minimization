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
exp = 2
params = EXPERIMENTS[exp]
digitsInMultiplicand1, digitsInMultiplicand2, product = params

# Default assumption parameters
num_assumptions = 0
limit_assumptions = 1
# Use standard carry equations by default
qubit_reduction_method = 2

#We can override the digit and product values above using arguments
if len(sys.argv) > 2:
    digitsInMultiplicand1 = int(sys.argv[1])
    digitsInMultiplicand2 = int(sys.argv[2])
    product = int(sys.argv[3])
    qubit_reduction_method = int(sys.argv[4])
    num_assumptions = int(sys.argv[5])
    limit_assumptions = int(sys.argv[6])
    OutputFileName = str(sys.argv[7])

equation_generator, coef_str_generator = QUBIT_REDUCTION_ID[qubit_reduction_method]

eqns = equation_generator(digitsInMultiplicand1, digitsInMultiplicand2, product)

s = time()
# None means the result will be printed to screen
output = None#OutputFileName


# We can use the handy state caching    
cache_name = None#'_state_{}'.format(str(product)[-6:])
log_deductions = False
if cache_name is not None:
    try:
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

try:
    coef_filename = OutputFileName.replace('.txt', '_coef.txt')
    coef_str = coef_str_generator(system.final_equations)
    coef_str_to_file(coef_str, coef_filename)
except Exception as e:
    print 'Failed to write the coefficient'
    print e


check_solutions(product, system.solutions.copy(), verbose=True)

print

## Now lets do the assumptions stuff
if len(system.unsolved_var) and num_assumptions:    

    solns = zip(*make_simultaneous_assumptions(system, 
                                          num_assumptions=num_assumptions,
                                          verbose=True,
                                          rank_func=max_coef_rank_variables,
                                          return_variables=True,
                                          limit_permutations=limit_assumptions))
    
    for i, sol in enumerate(solns):
        sol, sub = sol
        print '\n' + 'Case {}'.format(i + 1)
        print sub
        #sol.print_summary()
        print 'Num Qubits End: {}'.format(len(sol.unsolved_var))
        
        correct = check_solutions(product, sol.solutions.copy(), verbose=True)
