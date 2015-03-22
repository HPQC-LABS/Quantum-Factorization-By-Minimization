# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 19:25:24 2015

@author: Richard
"""

import itertools

import numpy
import sympy

from carry_equations_generator import generate_carry_equations
from semiprime_tools import num_to_factor_num_qubit
from sympy_helper_fns import (max_value, min_value, is_equation,
                              remove_binary_squares_eqn, balance_terms,
                              cancel_constant_factor, is_constant,
                              num_add_terms, parity, is_monic, is_one_or_zero,
                              remove_binary_squares, expressions_to_variables,
                              gather_monic_terms, square_equations,
                              str_eqns_to_sympy_eqns, standardise_equation,
                              is_simple_binary)
from verification import check_solutions
from sympy_solver import EquationSolver
from sympy.core.cache import clear_cache

class EquationSolver2(EquationSolver):
    ''' Stripped down version '''

    def solve_equations(self, max_iter=250, verbose=False):
        ''' Solve a system of equations
        '''
        state_summary = self._length_tuple
        # The number of iterations in which we've made no new deductions
        num_constant_iter = 0

        # Keep track of the last few number of solutions, so that if we get
        # too many repeats we can break the cycle
        prev_num_solns = []

        if verbose:        
            self.print_('Num variables: {}'.format(len(self.variables)))
            self.print_('Iter\tNum Eqn\tNum Ded\tNum Sol')
        for i in xrange(max_iter):
            # Check we're not going around in circles
            prev_num_solns.append(self._length_tuple[2])
            if len(prev_num_solns) > 10:
                comp = prev_num_solns.pop(0)
                if all([comp-1 <= v <= comp+1 for v in prev_num_solns]):
                    break            

            # Clear the cache so that we don't blow up memory when working with
            # large numbers
            clear_cache()

            if verbose:
                self.print_('\t'.join(['{}'] * 4).format(i, *state_summary))

            self.equations = self.clean_equations(self.equations)
            
            # Extract all the equations from the system
            all_equations = self.final_equations

            all_equations.extend(self.non_trivial_solns)

            self.apply_judgements(all_equations)
            self.apply_contradictions(all_equations)

            # Slightly mysterious clean that fixes judgement blow up.
            # Something to do with the way clean_solutions cleans cycling imports,
            # cleaning re-updates values in the deductions and such.
            self.clean_deductions()

            if self._length_tuple == state_summary:
                num_constant_iter += 1

                # Here lets apply some slower, complex judgements to try and
                # unstick ourselves
                self.apply_judgements_complex(all_equations, num_constant_iter,
                                              verbose=verbose)

                if num_constant_iter > 4 or (self._length_tuple[:2] == (0, 0)):
                    break
                
                self.clean_deductions()

            if self._length_tuple != state_summary:
                num_constant_iter = 0
                state_summary = self._length_tuple

        # Final clean again, for good luck
        self.equations = self.clean_equations(self.equations)
        # and clear the cache for future generations
        clear_cache()
        
        # Close the pool
        self.close_pool()

    def apply_judgements_complex(self, equations, num_constant_iter, 
                                 verbose=False):
        ''' Apply more complex or slow judgements if we get stuck.
            num_constant_iter is the number of iterations that we have been
            stuck for.
        '''
        equations = self.equations
        if num_constant_iter == 0:
            return
        state_summary = self._length_tuple
        
        if num_constant_iter > 0:
            # Use mini-assumptions
            # Do for every stuck iteration, but with more variables
            for eqn in equations:
                # Limit the substitutions at 2^6=64
                num_var = min(3*num_constant_iter + 2, 6) + 2
                # Rank by number of times each occurs
                self.judgement_mini_assumption(eqn, num_var=num_var, 
                                               coef_transform=lambda x: pow(x, 0.01),
                                                log_flag='1')
                # Rank by sum of coefficients
#                self.judgement_mini_assumption(eqn, num_var=num_var,
#                                               coef_transform=lambda x: pow(x, 1.01), log_flag='2')

        if num_constant_iter == 1:
            # Use the multi-equation mini-assumption to find simple
            # correspondences between simple equations
            num_var = 5
            max_num_eqn = 3
#            filter_func = lambda eq: eq.atoms(sympy.Symbol) < 5
            filter_func = lambda eq: (num_add_terms(eq.lhs) + num_add_terms(eq.rhs)) < 5
            short_equations = filter(filter_func, equations)
            # Make sure we do something with the short equations
            if len(short_equations) < max_num_eqn:
                eqn_comb = [short_equations]
            else:
                eqn_comb = itertools.combinations(short_equations, max_num_eqn)

            for comb in eqn_comb:
                self.judgement_mini_assumption_multi_eqn(comb, num_var=num_var,
                                                         coef_transform=lambda x: pow(x, 0.01),
                                                         cutoff=0.01, log_flag='3')
#                self.judgement_mini_assumption_multi_eqn(comb, num_var=num_var,
#                                                         coef_transform=lambda x: pow(x, 1.01),
#                                                         cutoff=0.01, log_flag='4')
            
            for eqn in equations:

                # Apply the judgements that may add complexity
                self.judgement_5(eqn, increase_complexity=True)
                self.judgement_6(eqn, increase_complexity=True)
                self.judgement_9(eqn, increase_complexity=True)

        if num_constant_iter > 1:
            
            ## Work on the whole thing
            # Unleash the multiple equation mini-assumption
            #TODO put in a config or improve filtering
            num_var = num_constant_iter * 2#min(3 * num_constant_iter + 2, 6)
            num_eqn = 1 + int(num_constant_iter / 2.0)#max(2, int(num_constant_iter / 2.0)) + 1
            num_eqn = min(len(equations), num_eqn)
            eqn_comb = itertools.combinations(equations, num_eqn)
            for comb in eqn_comb:
                self.judgement_mini_assumption_multi_eqn(comb, num_var=num_var,
                                                         coef_transform=lambda x: pow(x, 0.01),log_flag='5')
#                self.judgement_mini_assumption_multi_eqn(comb, num_var=num_var,
#                                                         coef_transform=lambda x: pow(x, 1.01),log_flag='6')


            ## Work on the tails
#            # Unleash the multiple equation mini-assumption
#            #TODO put in a config or improve filtering
##            num_var = min(3 * num_constant_iter + 2, 6) + 10
#            tail_size = 1 + num_constant_iter
#            taill = equations[:tail_size]
#            tailr = equations[-tail_size:]
##            num_varl = len(expressions_to_variables(taill)) - 5 + num_constant_iter
##            num_varr = len(expressions_to_variables(tailr)) - 5 + num_constant_iter
#            num_varl = 3 + num_constant_iter*2
#            num_varr = 3 + num_constant_iter*2
#            print num_varl
#            print num_varr
#            self.judgement_mini_assumption_multi_eqn(taill, num_var=num_varl,
#                                                     coef_transform=lambda x: pow(x, 0.01),log_flag='5')
#            self.judgement_mini_assumption_multi_eqn(tailr, num_var=num_varr,
#                                                     coef_transform=lambda x: pow(x, 0.01),log_flag='6')


#    def sequential_search(equations):
#        for eqn in equations:

    @property
    def non_trivial_solns(self):
        ''' A list of everything from self.solutions that might hold 
            information
            
            >>> a, b, c, u, v, x, y, z = sympy.symbols('a b c u v x y z')
            >>> system = EquationSolver()
            >>> solutions = {a: 1, b: c, u: 1 - v, x: y*z, y: x - 2*z}
            >>> system.solutions = solutions
            >>> system.non_trivial_solns
            [y + 2*z == x, x == y*z]
        '''
        # Now fetch non-trivial solutions
        non_trivial_soln = []
        for variable, soln in self.solutions.iteritems():
            # If we've found the solution, don't bother trying to apply
            # judgements to it                
            if is_simple_binary(soln):
                continue
            non_trivial_soln.append(sympy.Eq(variable, soln))
        non_trivial_soln = map(remove_binary_squares_eqn, non_trivial_soln)
        non_trivial_soln = map(balance_terms, non_trivial_soln)
        non_trivial_soln = map(cancel_constant_factor, non_trivial_soln)
        non_trivial_soln = filter(is_equation, non_trivial_soln)
        non_trivial_soln = itertools.ifilterfalse(is_simple_binary, non_trivial_soln)
        return list(non_trivial_soln)

    @property
    def deductions_as_equations(self):
        ''' Return deductions as a list of equations '''
        new_equations = EquationSolver._dict_as_equations(self.deductions)
        new_equations = filter(is_equation, new_equations)        
        new_equations = [eqn.expand() for eqn in new_equations]
        new_equations = map(remove_binary_squares_eqn, new_equations)
        new_equations = map(balance_terms, new_equations)
        new_equations = map(cancel_constant_factor, new_equations)
        new_equations = filter(is_equation, new_equations)
        new_equations = itertools.ifilterfalse(is_simple_binary, new_equations)
        return list(new_equations)

    @property
    def final_equations(self):
        ''' final_equations are the final filtered equations that also
            include deductions
        '''
        final_equations = self.equations + self.deductions_as_equations
        return final_equations

    def apply_judgements(self, equations):
        ''' Apply judgements to a list of sympy equations and directly update
            self.deductions
        '''
        for eqn in equations:
            self.judgement_0(eqn)
            self.judgement_prod(eqn)
            self.judgement_min_max(eqn)
            self.judgement_1(eqn)
            self.judgement_2(eqn)
            self.judgement_3(eqn)
            self.judgement_4(eqn)
            self.judgement_5(eqn, increase_complexity=False)
            self.judgement_6(eqn, increase_complexity=False)
            self.judgement_7(eqn)
            self.judgement_8(eqn)
            self.judgement_9(eqn)

    def print_summary(self):
        ''' Print a summary of the information held in the object '''
        unsolved_var = self.unsolved_var
        if self.log_deductions:
            self.print_deduction_log()

        self.print_('Num Qubits Start: {}'.format(self.num_qubits_start))
        self.print_('Num Qubits End: {}'.format(len(unsolved_var)), close=True)

        for s in self.solutions.iteritems(): print s


if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    prod = 143
    prod = 56153
#    prod = 4306239659
    fact1, fact2 = num_to_factor_num_qubit(prod)
    equations = generate_carry_equations(fact1, fact2, prod)
    equations = filter(is_equation, equations)
    
#    var_dict = {str(v): v for v in expressions_to_variables(equations)}
    system = EquationSolver2(equations, log_deductions=True)

    system.solve_equations(verbose=True)
    
    system.print_summary()

#    for e in equations:
#        system.equations.append(e)
#        system.solve_equations()
#        print system._length_tuple
    
#    system = Eq(equations)
#    system.run(stuck_iter=100, var_step=4)
    
#    check_solutions(prod, system.state_dict, verbose=True)
#    print system.state_dict
#    for e in system.equations: print e, is_equation(e), is_false(e)