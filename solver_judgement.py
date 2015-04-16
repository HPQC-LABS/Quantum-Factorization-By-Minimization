# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 19:25:24 2015

@author: Richard
"""

from copy import deepcopy
import itertools

import sympy

from carry_equations_generator import generate_carry_equations
from contradictions import apply_contradictions
from judgement_mixin import JudgementMixin
from semiprime_tools import num_to_factor_num_qubit
from sympy_helper_fns import (is_equation, expressions_to_variables,
                              str_eqns_to_sympy_eqns, standardise_equation,
                              is_simple_binary, dict_as_eqns, degree,
                              eqns_with_variables, num_add_terms,
                              remove_binary_squares, square_equations,
                              is_constant)
from solver_base import BinarySolutionSolverBase, unique_array_stable
from solver_sequential import SolverSequential
from sympy.core.cache import clear_cache

class SolverJudgement(BinarySolutionSolverBase, JudgementMixin):
    ''' A solver that uses the vanilla judgements to the best of it's ability
    '''

    CONST_ITER_LIMIT = 1

    def __init__(self, *args, **kwargs):
        super(SolverJudgement, self).__init__(*args, **kwargs)

        # Holding place for deductions that aren't quite simple enough to be
        # solutions
        self.deductions = {}

    @property    
    def _length_tuple(self):
        ''' Return a tuple of the lengths of equations, deductions, solutions 
        '''
        return len(self.equations), len(self.solutions), len(self.deductions)

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
            self.print_('Iter\tNum Eqn\tNum Sol\tNum Ded')
        for i in xrange(max_iter):
            # Check we're not going around in circles
            prev_num_solns.append(self._length_tuple[1])
            if len(prev_num_solns) > 20:
                comp = prev_num_solns.pop(0)
                if all([comp-1 <= v <= comp+1 for v in prev_num_solns]):
                    break

            # Clear the cache so that we don't blow up memory when working with
            # large numbers
            clear_cache()

            if verbose:
                self.print_('\t'.join(['{}'] * 4).format(i, *state_summary))

            # First do some spring cleaning
            self.clean_deductions()
            self.equations = self.clean_equations(self.equations)

            # Extract all the equations from the system
            all_equations = self.final_equations

            self.apply_judgements(all_equations)
            self.equations.extend(self.match_eqns(all_equations))
            apply_contradictions(all_equations)

            # Re clean the deductions after applying the judgements
            self.clean_deductions()

            if self._length_tuple == state_summary:
                num_constant_iter += 1

                # Here lets apply some slower, complex judgements to try and
                # unstick ourselves
                self.apply_judgements_complex(equations=self.equations[:],
                                              deductions=self.deductions.copy(),
                                              num_constant_iter=num_constant_iter,
                                              verbose=verbose)

                self.clean_deductions()

                if (num_constant_iter > self.CONST_ITER_LIMIT) or (len(self.unsolved_var) < 2):
                    break


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
        to_sub = {}

        # First copy the deductions
        to_sub = self.deductions.copy()

        # Extract solutions using the magic __getitem__ attribute
        var_of_interest = expressions_to_variables(cleaned + dict_as_eqns(self.deductions))
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
            >>> system = SolverJudgement([], {str(v) : v for v in variables})
            >>> ZERO, ONE = sympy.sympify(0), sympy.sympify(1)
            >>> deductions = {a: ONE, b: ZERO, ONE: c, x: 1 - y, z*x: ONE}
            >>> system.deductions = deductions
            >>> system.clean_deductions()
            >>> system.solutions
            {c: 1, x: -y + 1, b: 0, a: 1}
            >>> system.deductions
            {-y*z + z: 1}

            >>> variables = [a, b, c, x, y, z]
            >>> system = SolverJudgement([], {str(v) : v for v in variables})
            >>> deductions = {a: a*b, b: a*b, a: b}
            >>> system.deductions = deductions
            >>> system.clean_deductions()
            >>> system.solutions
            {a: b}

            Sort out the x = xy case
            >>> a, b, c, x, y, z = sympy.symbols('a b c x y z')
            >>> variables = [a, b, c, x, y, z]
            >>> system = SolverJudgement([], {str(v) : v for v in variables})
            >>> ZERO, ONE = sympy.sympify(0), sympy.sympify(1)
            >>> deductions = {a: a*b, x:x*y + y*z}
            >>> system.deductions = deductions
            >>> system.clean_deductions()
            >>> system.solutions
            {}
            >>> system.deductions
            {a*b: a, x*y + y*z: x}
        '''
        # First extract solutions so that we have maximum substituting power
        for expr, val in self.deductions.copy().iteritems():
            # If it's simple binary, then add it to the solutions
            if is_simple_binary(expr) and is_simple_binary(val):
                self.solutions[expr] = val
                self.deductions.pop(expr)

        # Clean the solutions so we don't spend so long in subs
        # Extract only the atoms we would like to try and find
        ded_as_eqn = self.deductions_as_equations
        atoms = expressions_to_variables(ded_as_eqn)
        cleaned_sol = {a: self.solutions[a] for a in atoms if (a != self.solutions[a])}

        old_deductions = self.deductions.copy()
        self.deductions = {}
        for expr, val in old_deductions.iteritems():
            # Sub in all of the solutions we know
            new_expr = remove_binary_squares(expr.subs(cleaned_sol).expand())
            new_val = remove_binary_squares(val.subs(cleaned_sol).expand())

            # If equal, continue
            if new_expr == new_val:
                continue

            # Make sure we're always reducing the qubit interactions in our
            # quest for simplicity
            if degree(new_expr) < degree(new_val):
                new_expr, new_val = new_val, new_expr
            # Finally, add it back in
            self.deductions[new_expr] = new_val



    def apply_judgements_complex(self, equations, deductions, num_constant_iter,
                                 verbose=False):
        ''' Apply more complex or slow judgements if we get stuck.
            num_constant_iter is the number of iterations that we have been
            stuck for.
        '''
        if num_constant_iter == 0:
            return

        # First make sure we don't have any duplicates, since these will make
        # searching much longer
        equations = unique_array_stable(equations)

        # First unleash the mini-assumption. This is checking single-equation
        # consistency
        if num_constant_iter == 1:
            # Use mini-assumptions
            all_eqns = equations + dict_as_eqns(deductions)            
            for eqn in all_eqns:
                # Limit the substitutions at 2^8=256
                num_var = 8
                # Rank by number of times each occurs
                self.judgement_mini_assumption(eqn, num_var=num_var,
                                               coef_transform=lambda x: pow(x, 1.01))        
        
            ## Apply square judgements
            for method in (1, 2):
                sqrd = square_equations(all_eqns, method=method, term_limit=10)
                self.apply_judgements(sqrd)

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
        final_equations = self.equations[:] + self.deductions_as_equations
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
            self.judgement_3(eqn)
            self.judgement_4(eqn)
            self.judgement_5(eqn, increase_complexity=True)#False)
            self.judgement_6(eqn, increase_complexity=True)#False)
            self.judgement_7(eqn)
            self.judgement_8(eqn)
            self.judgement_9(eqn)

            self.judgement_two_term(eqn)
            
#            self.judgement_2_extended(eqn)
#            self.judgement_8_extended(eqn)

    def match_eqns(self, equations):
        ''' Given a list of equations, see if any of the LHS == RHS so we can
            make new deductions
            
            >>> eqns = ['p4 + q4 + z34 == 2*z9394 + 1',
            ...         'p4 + q4 + z34 == 2*z45 + 1']
            >>> eqns = str_eqns_to_sympy_eqns(eqns)
            >>> solver = SolverJudgement()
            >>> solver.match_eqns(eqns)
            [z9394 == z45]
        '''
        EQUATION_QUADRATIC_LIMIT = 300
        new_eqns = []
        if len(equations) < EQUATION_QUADRATIC_LIMIT:
            to_add = []
            # Now add any equations where LHS = RHS1, LHS = RHS2 and permutations
            def _helper(eqn1, eqn2, to_add):
                if ((eqn1.lhs == eqn2.lhs) and
                    (eqn1.rhs != eqn2.rhs) and
                    (not is_constant(eqn1.lhs))):
                    new_eq = sympy.Eq(eqn1.rhs, eqn2.rhs)
                    new_eq = standardise_equation(new_eq)
                    
                    to_add.append(new_eq)
#                    self.print_('Equation added! {}, {}\t=>\t{}'.format(eqn1, eqn2, new_eq))
    
            for eqn1, eqn2 in itertools.combinations(equations, 2):
                _helper(eqn1, eqn2, to_add)
                _helper(sympy.Eq(eqn1.rhs, eqn1.lhs), eqn2, to_add)
                _helper(eqn1, sympy.Eq(eqn2.rhs, eqn2.lhs), to_add)
                _helper(sympy.Eq(eqn1.rhs, eqn1.lhs), sympy.Eq(eqn2.rhs, eqn2.lhs),
                        to_add)
            to_add = filter(is_equation, to_add)
            new_eqns.extend(to_add)
        
        # Now filter out new_eqns that are already in eqns
        seen = set(equations[:])
        seen_add = seen.add
        new_eqns = [x for x in new_eqns if not (x in seen or seen_add(x))]
        
        return new_eqns

    def copy(self):
        ''' Return a new instance of itself '''
        cls = type(self)
        copy = cls(deepcopy(self.equations),
                   deepcopy(self.variables),
                   output_filename=self.output_filename,
                   parallelise=self.parallelise)

        # Now use deepcopy to copy everything else
        copy.num_qubits_start = self.num_qubits_start
        copy.deductions = deepcopy(self.deductions)
        copy.solutions = deepcopy(self.solutions)
        copy._fix_pq_soln_used = self._fix_pq_soln_used
        return copy

    # Pickling
    # This isn't mature/finished, but manages to write equations, deductions and
    # solutions to disk
    def __getstate__(self):
        return (self.equations, self.solutions, self.variables, self.deductions,
                self.num_qubits_start, self.output_filename, self.parallelise, 
                self._fix_pq_soln_used)
        
    def __setstate__(self, state):
        (self.equations, self.solutions, 
         self.variables, self.deductions, self.num_qubits_start,
         self.output_filename, self.parallelise, 
         self._fix_pq_soln_used) = state
        
        self._pool = None

if __name__ == '__main__':
    import doctest
    doctest.testmod()
