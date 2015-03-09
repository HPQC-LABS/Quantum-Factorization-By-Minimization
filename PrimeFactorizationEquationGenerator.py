#!/usr/bin/env python


import sys
from time import time

from cfg_sympy_solver import EXPERIMENTS, QUBIT_REDUCTION_ID, EXPERIMENTS_21
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

def num_to_factor_num_qubit(prod):
    ''' Given a number, work out how many qubits the factors should be.
        Return the largest first
    
        >>> num_to_factor_num_qubit(143)    
        (4, 4)
        >>> num_to_factor_num_qubit(56153)
        (8, 8)
        >>> num_to_factor_num_qubit(1099551473989)
        (21, 21)
        >>> num_to_factor_num_qubit(1267650600228508624673600186743)
        (51, 51)
    '''
    bin_str = bin(prod)[2:]
    num_qub = len(bin_str)
    if num_qub % 2:
        return (num_qub + 1) / 2, (num_qub + 1) / 2
    else:
        return num_qub / 2, num_qub / 2

def run_experiment(exp_num, *args, **kwargs):
    # A default experiment to run
    params = EXPERIMENTS[exp_num]
    digitsInMultiplicand1, digitsInMultiplicand2, product = params[:3]
    factorize(product, digitsInMultiplicand1=digitsInMultiplicand1,
              digitsInMultiplicand2=digitsInMultiplicand2, *args, **kwargs)

def factorize(product, digitsInMultiplicand1=None, digitsInMultiplicand2=None,
              num_assumptions=0, limit_assumptions=1, qubit_reduction_method=0,
              output=None, invariant_interactions_on_substitution=True):
    ''' Notes:
        output = None -> output is printed to screen
    '''
    log_deductions = False
    
    if digitsInMultiplicand1 is None:
        assert digitsInMultiplicand2 is None
        digitsInMultiplicand1, digitsInMultiplicand2 = num_to_factor_num_qubit(product)
    
    equation_generator, coef_str_generator = QUBIT_REDUCTION_ID[qubit_reduction_method]
    
    eqns = equation_generator(digitsInMultiplicand1, digitsInMultiplicand2, product)
    
    s = time()
    
#    # We can use the handy state caching    
#    cache_name = None#'_state_{}'.format(str(product)[-6:])
#    if cache_name is not None:
#        try:
#            system = EquationSolver.from_disk(cache_name)
#        except Exception as e:
#            print e
#            system = EquationSolver(eqns, output_filename=output, 
#                                                log_deductions=log_deductions,
#                                                invariant_interactions_on_substitution=invariant_interactions_on_substitution)
#            system.solve_equations(verbose=True)
#            system.to_disk(cache_name)
#    
#    # Do it normally
#    else:
    system = EquationSolver(eqns, output_filename=output, 
                                        log_deductions=log_deductions,
                                        invariant_interactions_on_substitution=invariant_interactions_on_substitution)
    system.solve_equations(verbose=True, max_iter=250)
        
    
    system.print_summary()
    print '\nSolved in {:.3f}s'.format(time() - s)
    
    #try:
    #    coef_filename = OutputFileName.replace('.txt', '_coef.txt')
    #    coef_str = coef_str_generator(system.final_equations)
    #    coef_str_to_file(coef_str, coef_filename)
    #except Exception as e:
    #    print 'Failed to write the coefficient'
    #    print e
    
    
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


if __name__ == '__main__':
    try:
        if len(sys.argv) == 2:
            product = int(sys.argv[1])
            factorize(product=product)
    
        #We can override the digit and product values above using arguments
        elif len(sys.argv) > 2:
            digitsInMultiplicand1 = int(sys.argv[1])
            digitsInMultiplicand2 = int(sys.argv[2])
            product = int(sys.argv[3])
            qubit_reduction_method = int(sys.argv[4])
            num_assumptions = int(sys.argv[5])
            limit_assumptions = int(sys.argv[6])
        #    OutputFileName = str(sys.argv[7])
            factorize(product=product, digitsInMultiplicand1=digitsInMultiplicand1,
                      digitsInMultiplicand2=digitsInMultiplicand2,
                      num_assumptions=num_assumptions,
                      limit_assumptions=limit_assumptions,
                      qubit_reduction_method=qubit_reduction_method)
        else:
            product = 143
            factorize(product)
    except Exception as e:
        print e