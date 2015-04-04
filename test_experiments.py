# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 19:48:44 2015

@author: Richard
"""

#!/usr/bin/env python
from semiprime_tools import num_to_factor_num_qubit

from time import time

from cfg_sympy_solver import (EXPERIMENTS, EXTENDED_EXPERIMENTS, EXPERIMENTS_20,
                              EXPERIMENTS_21, QUBIT_REDUCTION_ID)

from sympy_solver import EquationSolver
from solver_hybrid import SolverHybrid
from verification import check_substitutions, check_solutions


if __name__ == '__main__':
    solver = SolverHybrid  
    SKIP = [100]
    #for exp, params in sorted(EXTENDED_EXPERIMENTS.iteritems(), 
    #                          key=lambda x: x[1].num_qubits_expected,
    #                          reverse=True):
    for exp, params in EXTENDED_EXPERIMENTS.iteritems():
        try:
            if exp in SKIP:
                continue
    
            print '***\tExperiment {}\t***'.format(exp)
        
            # A default experiment to run
#            digitsInMultiplicand1, digitsInMultiplicand2, product, expected_qubits = params
            product = params.product
            expected_qubits = params.num_qubits_expected
            digitsInMultiplicand1, digitsInMultiplicand2 = num_to_factor_num_qubit(product)
            
            # Default assumption parameters
            # Use standard carry equations by default
            qubit_reduction_method = 0
            
            equation_generator, coef_str_generator = QUBIT_REDUCTION_ID[qubit_reduction_method]
            
            eqns = equation_generator(digitsInMultiplicand1, digitsInMultiplicand2, product)
            
            s = time()
            # None means the result will be printed to screen
            output = None#OutputFileName
            
            
            # We can use the handy state caching    
            cache_name = None#'_state_{}'.format(str(product)[-6:])
            log_deductions = False
            invariant_interactions_on_substitution = True
        
            system = solver(eqns, output_filename=output, 
                                                log_deductions=log_deductions,
                                                invariant_interactions_on_substitution=invariant_interactions_on_substitution,
                                                parallelise=False)
            system.solve_equations(verbose=False, max_iter=200)
    
            if len(system.unsolved_var) != expected_qubits:        
                print 'Num Qubits Start: {}'.format(system.num_qubits_start)
                print 'Num Qubits End: {}'.format(len(system.unsolved_var))
                print 'Num Qubits Expected: {}'.format(expected_qubits)
        
            print 'Solved in {:.3f}s'.format(time() - s)
                
#            check_solutions(product, system.solutions.copy(), verbose=True)
            check_substitutions(product, system.copy(), verbose=True)
    
        except Exception as e:
            print e
        
        print '\n' * 2