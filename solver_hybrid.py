# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 19:25:24 2015

@author: Richard
"""

import itertools

import sympy

from carry_equations_generator import generate_carry_equations
from contradictions import apply_contradictions
from semiprime_tools import num_to_factor_num_qubit
from sympy_helper_fns import (is_equation, expressions_to_variables,
                              str_eqns_to_sympy_eqns, standardise_equation,
                              is_simple_binary, dict_as_eqns, degree,
                              eqns_with_variables, num_add_terms,
                              remove_binary_squares, square_equations)
from solver_base import unique_array_stable
from solver_judgement import SolverJudgement
from solver_sequential import SolverSequential
from sympy.core.cache import clear_cache

class SolverHybrid(SolverJudgement):
    ''' A solver that combines the best of the judgement and CSP worlds '''

    CONST_ITER_LIMIT = 3

    def apply_judgements_complex(self, equations, deductions, num_constant_iter,
                                 verbose=False):
        ''' Apply more complex or slow judgements if we get stuck.
            num_constant_iter is the number of iterations that we have been
            stuck for.
        '''
        super(SolverHybrid, self).apply_judgements_complex(equations=equations,
                                                           deductions=deductions,
                                                           num_constant_iter=num_constant_iter,
                                                           verbose=verbose)

        # Do a sequential search starting from the left or right of the equations
        # Limit the number of states so that we don't go crazy
        # Limit the number of equations so we don't spend forever substituting
        # into equations that don't contain the first few variables

        # For now don't substitute into thousands of deductions, while we wait
        # for interleaving equation adding
        if num_constant_iter in (2, 3):

            # Work out which way we are slicing
            if not hasattr(self, '_slice_dir'):
                self._slice_dir = -1
            else:
                self._slice_dir *= -1
            
            #num_eq = 20
            max_states = 2 ** 14 + 2

            # Reverse the equations depending on what we did last time            
            eqn_to_search = equations[::self._slice_dir]

            deduction_eqns = dict_as_eqns(deductions)

            # Prune the deductions
            atoms_of_interest = expressions_to_variables(eqn_to_search)
            if atoms_of_interest:
                deduction_eqns = eqns_with_variables(deduction_eqns,
                                                     atoms_of_interest,
                                                     strict=True)
            # Now interleave them for maximum effectiveness!
            all_eqn = SolverSequential.interleave_equations(eqn_to_search,
                                                            deduction_eqns,
                                                            priority=0)

            self.judgement_sequential_search(all_eqn,
                                             max_states=max_states)

    def judgement_sequential_search(self, equations, max_states=600):
        ''' Use a sequential solver to solve a bunch of equations, start from
            the left, adding any deductions that share atoms with the equations
        '''
        search = SolverSequential()
        variables = set()
        for e in equations:
            search.add_equation(e)
            subbed = search.sub_var(None, max_states=max_states)

#            for method in (1, 2):
#                sqrd = square_equations([e], method=method)
#                for _e in sqrd:
#                    search.add_equation(_e)

            variables.update(subbed)

#            print search._length_tuple

            if ((max_states is not None) and
                (len(search.valid_states) > max_states)):
                break

#        for s in search.valid_states: print s
        deductions = search.get_solutions()

        for k in variables:
            v = deductions[k]
            if k == v:
                continue
            if num_add_terms(v) == 2:
                self.judgement_two_term(sympy.Eq(k - v + 1, 1))
            else:
                self.update_value(k, v)

    def judgement_sequential_search2(self, equations, max_states=600, sorter=None):
        ''' Use a sequential solver to solve a bunch of equations, start from
            the left, adding any deductions that share atoms with the equations
        '''
        if not equations:
            return
#        max_states = 100
#        equations = equations[:3]
        equations = equations[::-1]
        search = SolverSequential()
        from semiprime_tools import cmp_variables_left
        sorter = cmp_variables_left
#        search.reorder_variables_to_sub(sorter=sorter)
#        print search.variables_to_sub

        ## Add 2 equations to start off with
        search.add_equation(equations.pop())
        search.add_equation(equations.pop())
        variables = set()
        for i in xrange(100):
            subbed = search.sub_var(4)
            search.add_equation(equations.pop())

            if subbed:
                print subbed
            else:
                break

            variables.update(subbed)

            print search._length_tuple

            if ((max_states is not None) and
                (len(search.valid_states) > max_states)):
                break

        deductions = search.get_solutions()

        for k in variables:
            v = deductions[k]
            if k == v:
                continue
            if num_add_terms(v) == 2:
                self.judgement_two_term(sympy.Eq(k - v + 1, 1))
            else:
                self.update_value(k, v)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
