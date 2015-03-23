# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 19:25:24 2015

@author: Richard
"""

import itertools

import sympy

from carry_equations_generator import generate_carry_equations
from contradictions import apply_contradictions
from judgement_mixin import JudgementMixin
from semiprime_tools import num_to_factor_num_qubit
from sympy_helper_fns import (is_equation, expressions_to_variables,
                              str_eqns_to_sympy_eqns, standardise_equation,
                              is_simple_binary, dict_as_eqns, degree,
                              eqns_with_variables, num_add_terms)
from solver_base import BinarySolutionSolverBase, unique_array_stable
from solver_sequential import SolverSequential
from sympy.core.cache import clear_cache

class SolverHybrid(BinarySolutionSolverBase, JudgementMixin):
    ''' A solver that combines the best of the judgement and CSP worlds '''

    def __init__(self, *args, **kwargs):
        super(SolverHybrid, self).__init__(*args, **kwargs)
        
        # Holding place for deductions that aren't quite simple enough to be
        # solutions
        self.deductions = {}
    
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
            self.print_('Iter\tNum Eqn\tNum Sol')
        for i in xrange(max_iter):
            # Check we're not going around in circles
            prev_num_solns.append(self._length_tuple[1])
            if len(prev_num_solns) > 10:
                comp = prev_num_solns.pop(0)
                if all([comp-1 <= v <= comp+1 for v in prev_num_solns]):
                    break            

            # Clear the cache so that we don't blow up memory when working with
            # large numbers
            clear_cache()

            if verbose:
                self.print_('\t'.join(['{}'] * 3).format(i, *state_summary))

            self.equations = self.clean_equations(self.equations)
            
            # Extract all the equations from the system
            all_equations = self.final_equations

            self.apply_judgements(all_equations)
            apply_contradictions(all_equations)

            # Slightly mysterious clean that fixes judgement blow up.
            # Something to do with the way clean_solutions cleans cycling imports,
            # cleaning re-updates values in the deductions and such.
            self.clean_deductions()

            if self._length_tuple == state_summary:
                num_constant_iter += 1

                # Here lets apply some slower, complex judgements to try and
                # unstick ourselves
                self.apply_judgements_complex(self.equations, num_constant_iter,
                                              verbose=verbose)

                if (num_constant_iter > 4) or (self._length_tuple[0] == 0):
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

    def clean_equations(self, equations):
        ''' Substitute in values, filter, normalise, unique-ify '''
        cleaned = equations[:]
        to_sub = self.deductions.copy()

        # Extract solutions using the magic __getitem__ attribute
        var_of_interest = expressions_to_variables(cleaned + dict_as_eqns(to_sub))
        for var in var_of_interest:
            sol = self.solutions[var]
            if sol == var:
                continue
            to_sub[var] = sol

        cleaned = self.batch_substitutions(equations, to_sub)
        cleaned = itertools.ifilter(is_equation, cleaned)
        cleaned = itertools.imap(standardise_equation, cleaned)
        cleaned = itertools.ifilter(is_equation, cleaned)
        return unique_array_stable(cleaned)        

    def clean_deductions(self):
        ''' Clean our deductions. Involves caching solved values and rearranging
            some equations, now we can have negative variables and substitutions

            >>> a, b, c, x, y, z = sympy.symbols('a b c x y z')
            >>> variables = [a, b, c, x, y, z]
            >>> system = SolverHybrid([], {str(v) : v for v in variables})
            >>> ZERO, ONE = sympy.sympify(0), sympy.sympify(1)
            >>> deductions = {a: ONE, b: ZERO, ONE: c, x: 1 - y, z*x: ONE}
            >>> system.deductions = deductions
            >>> system.clean_deductions()
            >>> system.solutions
            {c: 1, x: -y + 1, b: 0, a: 1}

            >>> system.deductions
            {-y*z + z: 1}

            >>> variables = [a, b, c, x, y, z]
            >>> system = SolverHybrid([], {str(v) : v for v in variables})
            >>> deductions = {a: a*b, b: a*b, a: b}
            >>> system.deductions = deductions
            >>> system.clean_deductions()
            >>> system.solutions
            {a: b}
            
            Sort out the x = xy case
            >>> a, b, c, x, y, z = sympy.symbols('a b c x y z')
            >>> variables = [a, b, c, x, y, z]
            >>> system = SolverHybrid([], {str(v) : v for v in variables})
            >>> ZERO, ONE = sympy.sympify(0), sympy.sympify(1)
            >>> deductions = {a: a*b, x:x*y + y*z}
            >>> system.deductions = deductions
            >>> system.clean_deductions()
            >>> system.solutions
            {}
            >>> system.deductions
            {a*b: a, x*y + y*z: x}
        '''
        # First trawl through the deductions for definite solutions
        for expr, val in self.deductions.copy().iteritems():
            if is_simple_binary(expr) and is_simple_binary(val):
                self.deductions.pop(expr)
                self.solutions[expr] = val

        # Clean the solutions so we don't spend so long in subs
        # Extract only the atoms we would like to try and find
        ded_as_eqn = self.deductions_as_equations
        atoms = expressions_to_variables(ded_as_eqn)
        cleaned_sol = {a: self.solutions[a] for a in atoms if (a != self.solutions[a])}

        old_deductions = self.deductions.copy()
        self.deductions = {}
        for expr, val in old_deductions.iteritems():
            new_expr = expr.subs(cleaned_sol).expand()
            new_val = val.subs(cleaned_sol).expand()
            
            # If equal, continue
            if new_expr == new_val:
                continue
            if degree(new_expr) < degree(new_val):
                new_expr, new_val = new_val, new_expr
            self.deductions[new_expr] = new_val

    def apply_judgements_complex(self, equations, num_constant_iter, 
                                 verbose=False):
        ''' Apply more complex or slow judgements if we get stuck.
            num_constant_iter is the number of iterations that we have been
            stuck for.
        '''
        if num_constant_iter == 0:
            return
        
        equations = unique_array_stable(equations)
        
        # Do a sequential search starting from the left of the equations
        #TODO We should limit on the number of states
#        num_eq = 100 * num_constant_iter
        max_states = 2 ** (num_constant_iter + 8)
        if not (num_constant_iter % 2):
            eqn_to_search = equations#[:num_eq]
        # Or the right
        else:
            # Reverse the equations again so the solver has a better chance
            eqn_to_search = equations[::-1]#[-num_eq:][::-1]

        deductions = dict_as_eqns(self.deductions)
        # Prune the deductions
        deduction_eqns = eqns_with_variables(deductions, 
                                             expressions_to_variables(eqn_to_search),
                                             strict=True)

        self.judgement_sequential_search(deduction_eqns + eqn_to_search, max_states=max_states)

        # Now unleash the mini-assumption        
        if num_constant_iter > 2:
            # Use mini-assumptions
            # Do for every stuck iteration, but with more variables
            for eqn in equations + dict_as_eqns(self.deductions):
                # Limit the substitutions at 2^6=64
                num_var = min(2*num_constant_iter + 2, 8)
                # Rank by number of times each occurs
                self.judgement_mini_assumption(eqn, num_var=num_var, 
                                               coef_transform=lambda x: pow(x, 1.01))

    @property
    def deductions_as_equations(self):
        ''' Return deductions as a list of equations '''
        eqns = dict_as_eqns(self.deductions)
        eqns = map(standardise_equation, eqns)
        eqns = filter(is_equation, eqns)
        return unique_array_stable(eqns)

    @property
    def final_equations(self):
        ''' final_equations are the final filtered equations that also
            include deductions
        '''
        final_equations = self.equations + self.deductions_as_equations
        return unique_array_stable(final_equations)

    def update_value(self, expr, val):
        ''' Simple function that allows us to leverage the judgements from
            the original equation solver
        '''
        if isinstance(expr, int): expr = sympy.sympify(expr)
        if isinstance(val, int): val = sympy.sympify(val)
        if degree(val) > degree(expr):
            expr, val = val, expr
        if is_simple_binary(expr) and is_simple_binary(val):
            self.solutions[expr] = val
        else:
            self.deductions[expr] = val

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
            self.judgement_2_extended(eqn)
            self.judgement_3(eqn)
            self.judgement_4(eqn)
            self.judgement_5(eqn, increase_complexity=False)
            self.judgement_6(eqn, increase_complexity=False)
            self.judgement_7(eqn)
            self.judgement_8(eqn)
            self.judgement_9(eqn)
            
            self.judgement_two_term(eqn)

    def judgement_sequential_search(self, equations, max_states=600):
        ''' Use a sequential solver to solve a bunch of equations, start from
            the left, adding any deductions that share atoms with the equations
        '''
        search = SolverSequential()
        variables = set()
        for e in equations:
            search.add_equation(e)
            search.sub_var(None)
            
            variables.update(e.atoms(sympy.Symbol))
            
            if ((max_states is not None) and 
                (len(search.valid_states) > max_states)):
                break
        
        deductions = search.get_deductions()
        
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
