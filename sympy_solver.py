"""
Created on Fri Dec 26 12:35:16 2014

Solve a system of equations with binary variables

@author: Richard Tanburn
"""
from copy import deepcopy
from collections import defaultdict
import inspect
import itertools
from operator import itemgetter
from random import shuffle
import sympy
from sympy.core.cache import clear_cache

import ReHandler
from contradiction_exception import ContradictionException
from sympy_helper_fns import (max_value, min_value, is_equation,
                              remove_binary_squares_eqn, balance_terms,
                              cancel_constant_factor, is_constant,
                              num_add_terms, parity, is_monic, is_one_or_zero,
                              remove_binary_squares, expressions_to_variables,
                              gather_monic_terms, square_equations,
                              str_eqns_to_sympy_eqns, standardise_equation,
                              is_simple_binary)
from objective_function_helper import (equations_to_vanilla_coef_str, 
                                       equations_to_vanilla_objective_function,
                                       equations_to_auxillary_coef_str)

from sympy_paralleliser import paralellised_subs, get_pool, DEFAULT_NUM_PROCESSES


__author__ = "Richard Tanburn"
__credits__ = ["Richard Tanburn", "Nathaniel Bryans", "Nikesh Dattani"]
__version__ = "0.0.1"
__status__ = "Prototype"

# Maximum number of equations before quadratic equality checking kicks in
EQUATION_EQUAL_CHECK_LIMIT = 350


class EquationSolver(object):
    ''' Solver of equations '''

    def __init__(self, equations=None, variables=None, log_deductions=False,
                 output_filename=None, invariant_interactions_on_substitution=True,
                 parallelise=False):
        if variables is None:
            if equations is None:
                variables = {}
            else:
                equations = filter(is_equation, equations)
                variables = {str(v): v for v in expressions_to_variables(equations)}
        
        if equations is None:
            equations = []

        # Number of variables at the start
        self.num_qubits_start = len(variables)

        # Dict of string tag to sympy variable instance
        self.variables = variables
        # List of Equation instances to be solved
        self.equations = equations
        # Set of deductions we have made
        self.deductions = {}

        # Solutions. Subset of deductions, where the key is a single variable.
        self.solutions = {}

        # File to print to
        self.output_filename = output_filename
        self._file = None

        # And keep a nested dictionary of who made them, if we want
        self.log_deductions = log_deductions
        self.deduction_record = defaultdict(lambda : defaultdict(list))
        
        # if invariant_interactions_on_substitution is True, then only
        # substitute x = (1-y + z1) type deductions, where each term on the
        # RHS has 1 atom
        self.invariant_interactions_on_substitution = invariant_interactions_on_substitution
        
        # Allow parallelisation
        self.parallelise = parallelise
        self._pool = None

    def print_(self, output, close=False):
        ''' Print either to screen or a file if given '''
        if self.output_filename is None:
            print output
        else:
            if (self._file is None) or self._file.closed:
                self._file = open(self.output_filename, 'a')
            output = str(output)
            self._file.write(output + '\n')
            if close:
                self._file.close()

    def copy(self):
        ''' Return a new instance of itself '''
        copy = EquationSolver(deepcopy(self.equations), 
                              deepcopy(self.variables),
                              log_deductions=self.log_deductions, 
                              output_filename=self.output_filename,
                              parallelise=self.parallelise)
                              
        # Now use deepcopy to copy everything else
        copy.num_qubits_start = self.num_qubits_start
        copy.deductions = deepcopy(self.deductions)
        copy.solutions = deepcopy(self.solutions)
        copy.deduction_record = deepcopy(self.deduction_record)
        copy.invariant_interactions_on_substitution = self.invariant_interactions_on_substitution
        return copy

    # Pickling
    # This isn't mature/finished, but manages to write equations, deductions and
    # solutions to disk
    def __getstate__(self):
        return (self.equations, self.deductions, self.solutions, 
                self.invariant_interactions_on_substitution, self.log_deductions,
                # We can't pickle defaultdicts apparently
                dict(self.deduction_record), self.variables, self.num_qubits_start,
                self.output_filename, self._file, self.parallelise)
        
    def __setstate__(self, state):
        (self.equations, self.deductions, self.solutions, 
         self.invariant_interactions_on_substitution, self.log_deductions,
         deduction_record, self.variables, self.num_qubits_start,
         self.output_filename, self._file, self.parallelise) = state
         
        # Re cast to a defaultdict
        self.deduction_record = defaultdict(lambda : defaultdict(list))
        for k, v in deduction_record.iteritems():
            self.deduction_record[k] = v

    def to_disk(self, filename):
        ''' Write a state to disk '''
        if filename is None:
            return
        import pickle
        pickle.dump(self, open(filename, 'w'))
    
    @staticmethod
    def from_disk(filename, **kwargs):
        ''' Load from disk '''
        if filename is None:
            raise ValueError('Cannot load from filename None')
        import pickle
        data = pickle.load(open(filename, 'r'))
        return data

    @staticmethod
    def _dict_as_equations(dict_):
        ''' Return deductions as a list of equations '''
        new_equations = []
        for lhs, rhs in dict_.iteritems():
            new_equations.append(sympy.Eq(lhs, rhs))
        
#        new_equations = filter(is_equation, new_equations)        
#        new_equations = [eqn.expand() for eqn in new_equations]
#        new_equations = map(remove_binary_squares_eqn, new_equations)
#        new_equations = map(balance_terms, new_equations)
#        new_equations = map(cancel_constant_factor, new_equations)
#        new_equations = filter(is_equation, new_equations)

        return sorted(new_equations, key=lambda x: str(x))
    
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
        return new_equations
    
    @property
    def solutions_as_equations(self):
        ''' Return solutions as a list of equations '''
        return EquationSolver._dict_as_equations(self.solutions)

    def print_deduction_log(self):
        ''' Print the judgements and the deductions they have made '''
        to_skip = ['clean_deductions', 'clean_solutions']
        for judgement in sorted(self.deduction_record.keys()):
            if judgement in to_skip:
                continue
            ded_info = self.deduction_record[judgement]
            self.print_('\n' + judgement)
            for eqn, deds in ded_info.iteritems():
                eqn_str = str(eqn).ljust(25)
                ded_str = map(lambda (x, y) : '{}={}'.format(x, y), deds)
                ded_str = ', '.join(ded_str)
                self.print_('{}\t=>\t{}'.format(eqn_str, ded_str))
        self.print_('\n')

    @property    
    def _length_tuple(self):
        ''' Return a tuple of the lengths of equations, deductions, solutions 
        '''
        return len(self.equations), len(self.deductions), len(self.solutions)

    def solve_equations(self, max_iter=250, verbose=False):
        ''' Solve a system of equations
        '''
        state_summary = self._length_tuple
        # The number of iterations in which we've made no new deductions
        num_constant_iter = 0

        if verbose:        
            self.print_('Num variables: {}'.format(len(self.variables)))
            self.print_('Iter\tNum Eqn\tNum Ded\tNum Sol')
        for i in xrange(max_iter):
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
                self.apply_judgements_complex(all_equations, num_constant_iter)

                # Now apply judgements to the squares of the equations
                # Since applying judgements to the square of the equations 
                # doesn't change behaviour as we get stuck, just apply it once
                if num_constant_iter == 2:
                    self.apply_judgements_square(all_equations, verbose=verbose)

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

    @property
    def final_equations(self):
        ''' final_equations are the final filtered equations that also
            include deductions
        '''
        final_equations = self.equations + self.deductions_as_equations
        final_equations = sorted(set(final_equations), key=str)
        return final_equations

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
        return non_trivial_soln

    def print_summary(self):
        ''' Print a summary of the information held in the object '''
        unsolved_var = self.unsolved_var
#
#        self.final_variables = unsolved_var
#
#        self.final_equations = final_equations
#        self.final_solutions = self.solutions.copy()

        if self.log_deductions:
            self.print_deduction_log()

#        self.print_('Unsimplified equations')
#        for e in self.equations:
#            self.print_(e)
#        self.print_('Deductions')
#        for e in self.deductions_as_equations:
#            self.print_(e)

#        self.print_('Solns')
#        for k in sorted(self.solutions.keys(), key=str):
#            self.print_('{} = {}'.format(k, self.solutions[k]))

        self.print_('Final Variables')
        self.print_(unsolved_var)


        self.print_('Final Equations')
        for e in self.final_equations:
            self.print_(e)

        self.print_('Num Qubits Start: {}'.format(self.num_qubits_start))
        self.print_('Num Qubits End: {}'.format(len(unsolved_var)), close=True)

        # Print the p and q solution
        pqs = {}
        zs = {}
        self.print_('p, q Solutions')
        for var, sol in self.solutions.iteritems():
            svar = str(var)
            if svar.startswith('p') or svar.startswith('q'):
                pqs[var] = sol
            else:
                zs[var] = sol
        for k in sorted(pqs.keys(), key=str):
            self.print_('{} = {}'.format(k, pqs[k].subs(zs)))

#        self.print_('Final coefficients')
#        self.print_(equations_to_coef_string(self.final_equations), close=True)

    @property
    def unsolved_var(self):
        ''' Return a set of variables we haven't managed to eliminate '''
        return set(self.variables.values()).difference(self.solutions.keys())

    @property
    def objective_function(self):
        ''' Return the final objective function, using self.final_equations '''
        return equations_to_vanilla_objective_function(self.equations)

    def objective_function_to_file(self, filename=None):
        ''' Write the objective function to a file, or printing if None.
            Also include the dictionary of variable number to original variable
        '''
        out = equations_to_vanilla_coef_str(self.equations)
        #out = equations_to_auxillary_coef_str(self.equations)

        if filename is None:
            self.print_(out)
        else:
            f = open(filename, 'a')
            f.write(out)
            f.close()

    def close_pool(self):
        ''' Close the pool if it's open and re-assign to None '''
        if self._pool is not None:
            self._pool.close()
            self._pool.join()
            self._pool = None

    def batch_substitutions(self, equations, substitutions):
        ''' Helper method that substitutes large dicts into large sets of
            equations. Deals with memory management and parallelisation.
            
            >>> x, y, z = sympy.symbols('x y z')
            >>> system = EquationSolver()
            >>> eqns = [x + y - 1,
            ...         x*z - 1,
            ...         x - 1]
            >>> eqns = map(balance_terms, map(sympy.Eq, eqns))
            >>> subs = {x: 1, z: 2}
            >>> subbed = system.batch_substitutions(eqns, subs)
            >>> for eqn in subbed: print eqn
            y + 1 == 1
            False
            True
        '''
        if len(equations) == 0:
            return []

        #TODO Move into a cfg file
        batch_size = 30
        min_batches = 6
        fill = None
        batch_equations = list(_batcher(equations, batch_size=batch_size, 
                                        fill_value=fill))
        # Reduce the final set so we only have the original equations
        last = batch_equations.pop()
        last = filter(lambda x: x is not None, last)
        batch_equations.append(last)        
        
        # substituted will be the holder for our new equations. It is None
        # while no method has worked.
        # Later it will be a nested list of batched results, which we need to
        # flatten later on
        substituted = None        
        
        # Try to parallelise the slow substitution
        if self.parallelise and (len(batch_equations) >= min_batches):
            try:
                if self._pool is None:
                    self._pool = get_pool()
                substituted = paralellised_subs(batch_equations, substitutions, 
                                                 pool=self._pool)
            except Exception as e:
                print e
                self.parallelise = False
                self.close_pool()

        if substituted is None:
            substituted = []
            for batch in batch_equations:
                substituted.append([eqn.subs(substitutions) for eqn in batch])
                clear_cache()

        # Now flatten substituted
        substituted = list(itertools.chain(*substituted))
        
        return substituted

    def clean_equations(self, eqns):
        ''' Remove True equations and simplify '''
        # First clean up the deductions so we can use them
        self.clean_deductions()

        cleaned = filter(is_equation, eqns[:])

        # Extract only the atoms we would like to try and find
        if len(cleaned):
            cleaned_atoms = expressions_to_variables(cleaned)
            cleaned_sol = ((var, self.solutions.get(var)) for var in cleaned_atoms)
            cleaned_sol = filter(lambda x: x[1] is not None, cleaned_sol)
            cleaned_sol = {x[0]: x[1] for x in cleaned_sol}
        else:
            cleaned_sol = {}

        # Combine all combinations into one dict, giving priority to cleaned_sol        
        combined_subs = self.deductions.copy()
        combined_subs.update(cleaned_sol)

        cleaned = self.batch_substitutions(cleaned, combined_subs)

        cleaned = filter(is_equation, cleaned)
        cleaned = [eqn.expand() for eqn in cleaned]
        cleaned = map(remove_binary_squares_eqn, cleaned)
        cleaned = map(balance_terms, cleaned)
        cleaned = map(cancel_constant_factor, cleaned)
        cleaned = filter(is_equation, cleaned)

        if len(cleaned) < EQUATION_EQUAL_CHECK_LIMIT:
            to_add = []
            # Now add any equations where LHS = RHS1, LHS = RHS2 and permutations
            def _helper(eqn1, eqn2, to_add):
                if ((eqn1.lhs == eqn2.lhs) and
                    (eqn1.rhs != eqn2.rhs) and
                    (not is_constant(eqn1.lhs))):
                    new_eq = sympy.Eq(eqn1.rhs, eqn2.rhs)
                    new_eq = balance_terms(new_eq)
                    
                    # Try only adding stuff with more than one additive term
                    if num_add_terms(new_eq.lhs) == num_add_terms(new_eq.rhs) == 1:
                        return
                    
                    to_add.append(new_eq)
#                    self.print_('Equation added! {}, {}\t=>\t{}'.format(eqn1, eqn2, new_eq))
    
            all_equations = itertools.chain(cleaned, self.deductions_as_equations)
            for eqn1, eqn2 in itertools.combinations(all_equations, 2):
                _helper(eqn1, eqn2, to_add)
                _helper(sympy.Eq(eqn1.rhs, eqn1.lhs), eqn2, to_add)
                _helper(eqn1, sympy.Eq(eqn2.rhs, eqn2.lhs), to_add)
                _helper(sympy.Eq(eqn1.rhs, eqn1.lhs), sympy.Eq(eqn2.rhs, eqn2.lhs),
                        to_add)
            to_add = filter(is_equation, to_add)
            cleaned.extend(to_add)

        return list(set(cleaned))

    def clean_deductions(self):
        ''' Clean our deductions. Involves caching solved values and rearranging
            some equations, now we can have negative variables and substitutions

            >>> a, b, c, x, y, z = sympy.symbols('a b c x y z')
            >>> variables = [a, b, c, x, y, z]
            >>> system = EquationSolver([], {str(v) : v for v in variables})
            >>> ZERO, ONE = sympy.sympify(0), sympy.sympify(1)
            >>> deductions = {a: ONE, b: ZERO, ONE: c, x: 1 - y, z*x: ONE}
            >>> system.deductions = deductions
            >>> system.clean_deductions()
            >>> system.solutions
            {c: 1, x: -y + 1, b: 0, a: 1}

            >>> system.deductions
            {z: y*z + 1}

            >>> variables = [a, b, c, x, y, z]
            >>> system = EquationSolver([], {str(v) : v for v in variables})
            >>> deductions = {a: a*b, b: a*b, a: b}
            >>> system.deductions = deductions
            >>> system.clean_deductions()
            >>> system.solutions
            {a: b}
            
            Sort out the x = xy case
            >>> a, b, c, x, y, z = sympy.symbols('a b c x y z')
            >>> variables = [a, b, c, x, y, z]
            >>> system = EquationSolver([], {str(v) : v for v in variables})
            >>> ZERO, ONE = sympy.sympify(0), sympy.sympify(1)
            >>> deductions = {a: a*b, x:x*y + y*z}
            >>> system.deductions = deductions
            >>> system.clean_deductions()
            >>> system.solutions
            {}
            >>> system.deductions
            {a*b: a, x: x*y + y*z}
        '''
        # First trawl through the deductions for definite solutions
        for expr, val in self.deductions.copy().iteritems():
            latoms = expr.atoms(sympy.Symbol)
            ratoms = val.atoms(sympy.Symbol)

            # Hack around the dodgy edge case xy = y
            if len(ratoms.intersection(latoms)):
                if len(latoms) == 1:
                    possible_other_value = self.deductions.get(val)
                    if possible_other_value is not None:
                        self.update_value(expr, possible_other_value)
                continue

            if (len(latoms) == 1) and is_monic(expr):
                self.deductions.pop(expr)
                curr_sol = self.solutions.get(expr)

                if (curr_sol is not None) and (curr_sol != val):
                    # We have different things. Better be careful!!
                    if is_constant(curr_sol):
                        # Both are constant and unequal
                        if is_constant(val):
                            err_str = 'clean_deductions: {} = {} != {}'.format(expr, curr_sol, val)
                            raise ContradictionException(err_str)
                        else:
                            # We have a variable and constant
                            self.update_value(val, curr_sol)
                    else:
                        # Once again, we have a constant and a value
                        if is_constant(val):
                            self.update_value(curr_sol, val)
                        # Both are symbolic
                        else:
                            self.update_value(curr_sol, _simplest(curr_sol, val))
                else:
                    if is_monic(expr):
                        self.solutions[expr] = val

            # If the RHS of a deduction is monic, then go again!
            elif (len(ratoms) == 1) and is_monic(val):
                self.deductions.pop(expr)
                curr_sol = self.solutions.get(val)

                # We might have some disagreement
                if (curr_sol is not None) and (curr_sol != expr):
                    # But if they're both symbolic that is ok!
                    if is_constant(curr_sol) and is_constant(expr):
                        err_str = 'clean_deductions: {} = {} != {}'.format(expr, curr_sol, val)
                        raise ContradictionException(err_str)
                    
                    elif is_constant(curr_sol):
                        self.update_value(val, curr_sol)
                    else:
                        # The new val is constant
                        self.solutions[val] = _simplest(curr_sol, expr)
                        self.update_value(curr_sol, val)
                else:
                    self.solutions[val] = expr

        # Now clean up the solutions before we plug them in
        self.clean_solutions()

        # Now go over the remaining deductions and pick out equations which
        # include an unsolved variables
        unsolved_var = self.unsolved_var

        # Clean the solutions so we don't spend so long in subs
        # Extract only the atoms we would like to try and find
        ded_as_eqn = self.deductions_as_equations
        if len(ded_as_eqn):
            cleaned_atoms = expressions_to_variables(ded_as_eqn)
            cleaned_sol = ((var, self.solutions.get(var)) for var in cleaned_atoms)
            cleaned_sol = filter(lambda x: x[1] is not None, cleaned_sol)
            cleaned_sol = {x[0]: x[1] for x in cleaned_sol}
        else:
            cleaned_sol = self.solutions.copy()

        old_deductions = self.deductions.copy()
        self.deductions = {}
        for expr, val in old_deductions.iteritems():
            
            # If we have something like x = xy + xz, we don't want to expand
            # up to something harder to deal with
            # Note here we want to break the convention of simplest, as in
            # the event of a tie, the first value is returned. In this case,
            # we want to err on the side of caution and not throw away
            # deductions
            # 
            # Actually, don't get rid of anything we don't *need* to, so that
            # we don't miss possible contradictions in the assumption stages.
#            if latoms.intersection(ratoms) and (_simplest(val, expr) == expr):
#                continue

            # Substitute all of the solved variables
            expr = expr.subs(cleaned_sol).expand()
            val = val.subs(cleaned_sol).expand()
            val = remove_binary_squares(val)
            latoms = expr.atoms()
            ratoms = val.atoms()
            if (len(latoms.intersection(unsolved_var)) or
                len(ratoms.intersection(unsolved_var))):

                eqn = sympy.Eq(expr, val)

                if eqn == True:
                    continue
                eqn = eqn.expand()
                eqn = remove_binary_squares_eqn(eqn)
                eqn = balance_terms(eqn)
                if eqn == True:
                    continue

                self.update_value(eqn.lhs, eqn.rhs)

            # Else we want to check consistency
            else:
                eqn = sympy.Eq(expr, val)

                if is_equation(eqn) and (not eqn):
                    raise ContradictionException('Subbing solutions raised contradiction in deductions')

                if not is_constant(expr):
                    assert not is_constant(val)
                    if expr != val:
                        self.print_('Dropping deduction {} = {}'.format(expr, val))

    def clean_solutions(self, _prev_changed=None):
        ''' Remove cycles and chains in the solutions. Make sure every value is
            set to equal an expression involving constants and unsolved variables
            
            _prev_changed allows the recursive calling to keep track of cycles

            Simple solutions and cleaning
            >>> system = EquationSolver()
            >>> a, b, c, x, y, z = sympy.symbols('a b c x y z')
            >>> ZERO, ONE = sympy.sympify(0), sympy.sympify(1)
            >>> soln = {a: ONE, b: a, c: a + b - ONE, x: ONE, z: ONE - x + y}
            >>> system.solutions = soln
            >>> system.clean_solutions()
            >>> system.solutions
            {c: 1, x: 1, b: 1, a: 1, z: y}

            Dealing with cyclic solutions
            >>> system = EquationSolver()
            >>> a, b, c, x, y = sympy.symbols('a b c x y')
            >>> soln = {a: b, b: c, c: a, x: y, y: x}
            >>> system.solutions = soln
            >>> system.clean_solutions()
            >>> system.solutions
            {c: b, a: b, y: x}

            Cyclic solutions that have a tail
            >>> system = EquationSolver()
            >>> a, b, c, x, y = sympy.symbols('a b c x y')
            >>> soln = {a: b, b: c, c: x, x: y, y: x}
            >>> system.solutions = soln
            >>> system.clean_solutions()
            >>> system.solutions
            {c: y, x: y, b: y, a: y}

            Non-trival cyclic solutions. Little bit rough around the edges.
            Keep an eye on it
            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> soln = {x: 1 - y, y: 1 - x}
            >>> system.solutions = soln
            >>> system.clean_solutions()
            >>> system.equations
            [x + y == 1]
            >>> system.deductions
            {}
            >>> system.solutions
            {y: -x + 1}

            >>> system = EquationSolver()
            >>> soln = {z: - x + 2, x: - y + 1, y: z - 1}
            >>> system.solutions = soln
            >>> system.clean_solutions()
            >>> system.equations
            [x + y == 1]
            >>> system.deductions
            {}
            >>> system.solutions
            {z: -x + 2, y: -x + 1}

            Incorrect keys
            >>> system = EquationSolver()
            >>> x, y = sympy.symbols('x y')
            >>> soln = {1 - x: 0, y: x}
            >>> system.solutions = soln
            >>> system.clean_solutions()
            >>> system.solutions
            {x: 1, y: 1}
            >>> system = EquationSolver()
            >>> x, y = sympy.symbols('x y')
            >>> soln = {1 - x: y}
            >>> system.solutions = soln
            >>> system.clean_solutions()
            >>> system.solutions
            {x: -y + 1}
        '''
        #TODO Make the solver handle cycles properly!!!

        # First make sure every key is a single value
        for expr, val in self.solutions.copy().iteritems():
            add_coef = expr.as_coeff_Add()[0]
            if add_coef:
                assert len(expr.atoms(sympy.Symbol)) == 1
                variable = expr - add_coef
                rhs = val - add_coef
                if variable.as_coeff_Mul()[0] < 0:
                    variable *= -1
                    rhs *= -1
                self.solutions.pop(expr)
                self.solutions[variable] = rhs
        
        # Now go through and do some standard checks and cleaning
        for variable, value in self.solutions.copy().iteritems():
            # Remove binary squares
            self.solutions[variable] = remove_binary_squares(value)            
            
            # Now make sure every value in the dict can be binary, throwing if
            # not
            if (max_value(value) < 0) or (min_value(value) > 1):
                err_str = 'clean_solutions: {} != {}'.format(variable, value)
                raise ContradictionException(err_str)
            # Now add an equation if it is written in terms of more than 3
            # other variables
#            if len(value.atoms(sympy.Symbol)) >= 3:
#                print 'Adding {} = {}'.format(variable, value)
#                self.update_value(variable, value)
#                self.solutions.pop(variable)

        changed = set()
        new_solutions = {}
        to_skip = []  # Skip for the infinite loops
        for variable, value in self.solutions.iteritems():
            # If value is a constant, skip
            if is_constant(value):
                assert len(variable.atoms()) == 1
                if not is_one_or_zero(value):
                    err_str = '{} must be binary, not {}'.format(variable, value)
                    raise ContradictionException(err_str)
                new_solutions[variable] = value
                continue
            # If x == x, remove it as it as there is no new information in there
            if variable == value:
                continue
            # Break the infinite chain!
            if variable in to_skip:
                continue

            init_value = value
            # Keep a track of the variables we've visited
            seen_before = [variable, value]
            for i in xrange(1000):
                old_value = value
                value = self.solutions.get(value)

                # Watch out for the infinite loops!
                if value in seen_before:
                    value = old_value
                    to_skip.append(value)
                    break
                else:
                    seen_before.append(value)

                if value is None:
                    value = old_value.subs(self.solutions, simultaneous=True).expand()
                    if value == variable:
                        value = old_value
                        changed.add(variable)
                    break
                elif isinstance(value, int):
                    break
                else:
                    continue#value = value.subs(self.solutions, simultaneous=True)

            if i > 990:
                raise ValueError('Invalid solutions, check it out!')

            # If we have x = xy, then remove this from solutions and put it in
            # deductions
            if len(variable.atoms(sympy.Symbol).intersection(value.atoms(sympy.Symbol))):
                self.update_value(variable, value)
                continue

            if value != init_value:
                changed.add(variable)

            new_solutions[variable] = value

        self.solutions = new_solutions

        if len(changed):
            #TODO Find a smarter way of choosing which one to pop?
            if changed == _prev_changed:
                var = changed.pop()
                self.equations.append(standardise_equation(sympy.Eq(var, self.solutions[var])))
                self.solutions.pop(var)
            self.clean_solutions(_prev_changed=changed)

    def _update_log(self, expr, value):
        ''' Log an update under a judgement '''
        if not self.log_deductions:
            return

        # Update the deductions process dictionary
        judgement, eqn = _get_judgement()
        self.deduction_record[judgement][eqn].append((expr, value))

    def update_value(self, expr, value):
        ''' Update the global dictionary and check for contradictions.
            Make sure expr is always 'positive'.
            NOTE Only accepts single terms for expr

            >>> system = EquationSolver()
            >>> x = sympy.symbols('x')
            >>> system.update_value(x, 0)
            >>> system.deductions
            {x: 0}
            >>> system.update_value(x, 1)
            Traceback (most recent call last):
                ...
            ContradictionException: x is already set to 0 != 1

            >>> system = EquationSolver()
            >>> x = sympy.symbols('x')
            >>> system.update_value(-x, 1)
            >>> system.deductions
            {x: -1}

            >>> system = EquationSolver()
            >>> x = sympy.symbols('x')
            >>> system.update_value(2*x, 0)
            >>> system.deductions
            {x: 0}

            x = x*y case
            >>> x, y, z = sympy.symbols('x y z')
            >>> system = EquationSolver()
            >>> system.update_value(x, x*y)
            >>> system.deductions
            {x*y: x}
            >>> system = EquationSolver()
            >>> system.update_value(y*z, x*y*z)
            >>> system.deductions
            {x*y*z: y*z}
        '''
        # First do some preprocessing to make sure the deduction is in a
        # reasonably nice form
        # If expr = 2*x and value == 0, then we can get rid of the 2
        if value == 0:
            expr = expr.as_coeff_Mul()[1]

        # Make sure the left is positive
        if expr.as_coeff_Mul()[0] < 0:
            expr = - expr
            value = - value

        # If value is an int, sympify it
        if isinstance(value, int):
            value = sympy.sympify(value)
        
        # Remember, we live in binary land
        expr = remove_binary_squares(expr)
        value = remove_binary_squares(value)

        current_val = self.deductions.get(expr)

        # If value already maps to expr, avoid the cycle!
        if self.deductions.get(value) == expr:
            return

        # If we have the nasty case where x = x*y, then we want to flip the
        # values round
        # Actually we do this whenever the lhs is a factor of the rhs
        expr_atoms = expr.atoms(sympy.Symbol)
        lhs_simpler = all(expr_atoms.issubset(term.atoms(sympy.Symbol)) 
                      for term in value.as_coefficients_dict().iterkeys())
        if lhs_simpler and len(expr_atoms) < len(value.atoms(sympy.Symbol)):
            self.update_value(value, expr)
            return

        # No possible conflict
        if current_val is None:
            self.deductions[expr] = value
            self._update_log(expr, value)

        # If we already know a value of this family
        elif is_one_or_zero(current_val):
            # If we've found a numeric value
            if is_constant(value):
                if current_val != value:
                    raise ContradictionException('{} is already set to {} != {}'.format(expr,
                                                 current_val, value))
                else:
                    return
            # We know what we're trying to update to is not a numeric, so update it too
            else:
                self.update_value(value, current_val)
#                self.deductions[expr] = current_val
                self._update_log(expr, value)

        # Current_val is symbolic
        else:
            if is_constant(value):
                # Perform another error check
                cc_val = self.deductions.get(current_val)
                if is_one_or_zero(cc_val) and (cc_val != value):
                    raise ContradictionException(
                            '{} is already set to {} != {}'.format(
                            current_val,
                            cc_val,
                            value))
                self.deductions[current_val] = value
                self.deductions[expr] = value
                self._update_log(current_val, value)
                self._update_log(expr, value)
            # Both values are symbolic!
            else:
                #TODO Clean up the hack around this silly edge case
                # Right now, if the RHS is written in terms of the LHS, then
                # we'd prefer to use the new value
                if (expr.atoms(sympy.Symbol).issubset(current_val.atoms(sympy.Symbol)) and
                    not expr.atoms(sympy.Symbol).issubset(value.atoms(sympy.Symbol))):
                    self.deductions[expr] = value
                    self._update_log(expr, value)

                else:
                    simple = _simplest(current_val, value)
                    self.deductions[expr] = simple
                    self._update_log(expr, simple)
                    if value != current_val:
                        self.deductions[current_val] = value
                        self._update_log(current_val, value)

    def get_var(self, var):
        ''' Return symbolic variable

            >>> system = EquationSolver()
            >>> x = system.get_var('x')
            >>> x
            x

            >>> x1 = system.get_var('x')
            >>> x1
            x

            >>> x1 is x
            True
        '''
        res = self.variables.get(var)
        if res is None:
            res = sympy.symbols(var, integer=True)
            self.variables[var] = res
        return res

    def set_to_max(self, expr):
        ''' Given an expression, update all terms so that it achieves it's maximum

            >>> system = EquationSolver()
            >>> expr = sympy.sympify('-2*x*y - 3 + 3*z23')
            >>> system.set_to_max(expr)
            >>> system.deductions
            {x*y: 0, z23: 1}
        '''
        coef_dict = expr.as_coefficients_dict()
        for var, coef in coef_dict.iteritems():
            if var == 1:
                continue
            if coef > 0:
                self.update_value(var, 1)
            elif coef < 0:
                self.update_value(var, 0)

    def set_to_min(self, expr):
        ''' Given an expression, update all terms so that it achieves it's minumum

            >>> system = EquationSolver()
            >>> expr = sympy.sympify('-2*x*y - 3 + 3*z23')
            >>> system.set_to_min(expr)
            >>> system.deductions
            {x*y: 1, z23: 0}

        '''
        coef_dict = expr.as_coefficients_dict()
        for var, coef in coef_dict.iteritems():
            if var == 1:
                continue
            if coef < 0:
                self.update_value(var, 1)
            elif coef > 0:
                self.update_value(var, 0)

    def apply_judgements_square(self, equations, verbose=False):
        ''' Pick out equations that we can square in a reasonable amount of
            time and apply the judgements to them
        '''
        pre = self._length_tuple
        eqn_sq1 = square_equations(equations, term_limit=20, method=1)
        self.apply_judgements(eqn_sq1)

        eqn_sq2 = square_equations(equations, term_limit=20, method=2)
        self.apply_judgements(eqn_sq2)
        post = self._length_tuple

        # If we didn't find anything, try the first round of complex judgements        
        if pre == post:
            self.apply_judgements_complex(eqn_sq1, num_constant_iter=1)
            self.apply_judgements_complex(eqn_sq2, num_constant_iter=1)

            post = self._length_tuple
            if verbose:
                num_ded = post[1] - pre[1]
                print '{} deductions made from squaring and complex judgements'.format(num_ded)

        elif verbose:
            num_ded = post[1] - pre[1]
            print '{} deductions made from squaring'.format(num_ded)

    def apply_judgements_complex(self, equations, num_constant_iter):
        ''' Apply more complex or slow judgements if we get stuck.
            num_constant_iter is the number of iterations that we have been
            stuck for.
        '''
        if num_constant_iter == 0:
            return
        state_summary = self._length_tuple
        
        if num_constant_iter > 0:
            # Use mini-assumptions
            for eqn in equations:
                # Limit the substitutions at 2^6=64
                num_var = min(3*num_constant_iter + 2, 6)
                # Rank by number of times each occurs
                self.judgement_mini_assumption(eqn, num_var=num_var, 
                                               coef_transform=lambda x: pow(x, 0.01))
                # Rank by sum of coefficients
                self.judgement_mini_assumption(eqn, num_var=num_var,
                                               coef_transform=lambda x: pow(x, 1.01))

        if num_constant_iter > 1:
            for eqn in equations:
                # Apply the slow judgement 8 and 2
                self.judgement_2_slow(eqn)
                self.judgement_8_slow(eqn)

                # Apply the judgements that may add complexity
                self.judgement_5(eqn, increase_complexity=True)
                self.judgement_6(eqn, increase_complexity=True)
                self.judgement_9(eqn, increase_complexity=True)

        if num_constant_iter > 2:

            for eqn in equations:
                # Only do 1 at a time, so if we have a new deduction
                # go round again
#                if self._length_tuple != state_summary:
#                    return

                self.judgement_n_term(eqn, num_constant_iter + 2)
        
        if num_constant_iter > 3:
#            if (not self.invariant_interactions_on_substitution):
            for eqn in equations:

#                if self._length_tuple != state_summary:
#                    return

                self.judgement_5(eqn, increase_complexity=True, 
                                 invariant_interactions_on_substitution=False)
                self.judgement_6(eqn, increase_complexity=True, 
                                 invariant_interactions_on_substitution=False)
                self.judgement_9(eqn, increase_complexity=True, 
                                 invariant_interactions_on_substitution=False)


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
            
    def apply_contradictions(self, equations):
        ''' Now look for contradictions in the equations '''
        for eqn in equations:
            self.contradiction_1(eqn)
            self.contradiction_2(eqn)

    def judgement_0(self, eqn):
        ''' Add x=y to deductions. This shouldn't be needed, but it's nice to
            make sure we're not missing anything obvious
            
            Also, if a*x = b*y, where a!=b are constants, x, y are variables,
            then x and y both have to be 0
            
            >>> eqns = ['6*p5*q7 == 5*q5*q7',
            ...         'x == 0',
            ...         'p6 == 2*z2627']
            >>> eqns = str_eqns_to_sympy_eqns(eqns)
            >>> system = EquationSolver()
            >>> for eqn in eqns: system.judgement_0(eqn)
            >>> system.deductions
            {p6: 0, x: 0, z2627: 0, p5*q7: 0, q5*q7: 0}
        '''
        lhs, rhs = eqn.lhs, eqn.rhs
        if len(lhs.atoms()) == len(rhs.atoms()) == 1:
            if lhs.is_constant():
                self.update_value(rhs, lhs)
            else:
                self.update_value(lhs, rhs)
        
        if num_add_terms(lhs) == num_add_terms(rhs) == 1:
            l_coef, l_var = lhs.as_coeff_Mul()
            r_coef, r_var = rhs.as_coeff_Mul()
            if (l_coef != 0) and (r_coef != 0) and (l_coef != r_coef):
                self.update_value(l_var, 0)
                self.update_value(r_var, 0)

    def judgement_prod(self, eqn):
        ''' If RHS is a non-zero constant and the LHS is a product of variables,
            then the variables must all be 1
            x*y*z=1 => x = y = 1

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x * y * z, 1)
            >>> system.judgement_prod(eqn)
            >>> system.deductions
            {x: 1, z: 1, y: 1}

            >>> eqn = sympy.Eq(2*x*y*z, 2)
            >>> system = EquationSolver([eqn])
            >>> system.solve_equations()
            >>> system.solutions
            {x: 1, z: 1, y: 1}
            
            >>> eqn = sympy.Eq(2*x*y*z, 1)
            >>> system = EquationSolver()
            >>> system.judgement_prod(eqn)
            Traceback (most recent call last):
                ...
            ContradictionException: judgement_sq: 2*x*y*z == 1
        '''
        if eqn.rhs == 1:
            if len(eqn.lhs.as_ordered_terms()) == 1:
                if eqn.lhs.as_coeff_mul()[0] != 1:
                    err_str = 'judgement_sq: {}'.format(str(eqn))
                    raise ContradictionException(err_str)
                for var in eqn.lhs.atoms(sympy.Symbol):
                    self.update_value(var, 1)

    def judgement_two_term(self, eqn):
        ''' If an expression has 2 variable terms, sub it in!
            This adds lots of complexity so the infrastructure, which is
            why we don't do it with 3 or 4 term sums.
            
            Note this isn't applied to the equations directly, but is called
            via the parity judgement
            
            >>> system = EquationSolver(invariant_interactions_on_substitution=False)
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y*z, 1)
            >>> system.judgement_two_term(eqn)
            >>> system.deductions
            {x: -y*z + 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 1)
            >>> system.judgement_two_term(eqn)
            >>> system.deductions
            {x: -y + 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(1 + y, 1)
            >>> system.judgement_two_term(eqn)
            >>> system.deductions
            {y: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x - y, 1)
            >>> system.judgement_two_term(eqn)
            >>> system.deductions
            {x: y + 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(-x + y, 1)
            >>> system.judgement_two_term(eqn)
            >>> system.deductions
            {x: y - 1}
        '''
        lhs, rhs = eqn.lhs, eqn.rhs
        if (num_add_terms(lhs) == 2) and is_constant(rhs):

            term1, term2 = lhs.as_ordered_terms()
            
            term1_atoms = term1.atoms(sympy.Symbol)
            term2_atoms = term2.atoms(sympy.Symbol)            
            
#            # max with 1 to avoid ignoring constants
#            if (self.invariant_interactions_on_substitution and 
#                (max((len(term1_atoms), 1)) != max((len(term2_atoms), 1)))):
#                return

            if ((0 < len(term2_atoms) < len(term1_atoms))
                or is_constant(term1)):
                self.update_value(term2, rhs - term1)

            else:
                self.update_value(term1, rhs - term2)

    def judgement_n_term(self, eqn, max_num_terms=3):
        ''' If an expression has n or fewer variable terms, sub it in!
            Only substitute single atoms terms for the moment

            >>> system = EquationSolver(invariant_interactions_on_substitution=False)
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y*z, 1)
            >>> system.judgement_n_term(eqn)
            >>> system.deductions
            {x: -y*z + 1}


            Get rid of the invariant_interactions_on_substitution test while
            it's turned off
#            >>> system = EquationSolver(invariant_interactions_on_substitution=True)
#            >>> x, y, z = sympy.symbols('x y z')
#            >>> eqn = sympy.Eq(x + y*z, 1)
#            >>> system.judgement_n_term(eqn)
#            >>> system.deductions
#            {}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, z)
            >>> system.judgement_n_term(eqn)
            >>> system.deductions
            {x: -y + z}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(1 + y, 1)
            >>> system.judgement_n_term(eqn)
            >>> system.deductions
            {y: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x - y, 1)
            >>> system.judgement_n_term(eqn)
            >>> system.deductions
            {x: y + 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(-x + y, 1)
            >>> system.judgement_n_term(eqn)
            >>> system.deductions
            {y: x + 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y + z, 1)
            >>> system.judgement_n_term(eqn)
            >>> system.deductions
            {x: -y - z + 1}

            >>> system = EquationSolver(invariant_interactions_on_substitution=False)
            >>> x, y, z, z2, u, v = sympy.symbols('x y z z2 u v')
            >>> eqn = sympy.Eq(2*x + y + z*z2, u + v)
            >>> system.judgement_n_term(eqn)
            >>> system.deductions
            {y: u + v - 2*x - z*z2}

            >>> system = EquationSolver()
            >>> lhs = sympy.sympify('q2 + q3')
            >>> rhs = sympy.sympify('2*q2*q3 + 2*z2021')
            >>> eqn = sympy.Eq(lhs, rhs)
            >>> system.judgement_n_term(eqn)
            >>> system.deductions
            {q2: 2*q2*q3 - q3 + 2*z2021}

            >>> system = EquationSolver()
            >>> lhs = sympy.sympify('q2 + 4*q3')
            >>> rhs = sympy.sympify('2*q2*q3 + 2*z2021')
            >>> eqn = sympy.Eq(lhs, rhs)
            >>> system.judgement_n_term(eqn)
            >>> system.deductions
            {q2: 0}
        '''
        lhs, rhs = eqn.lhs, eqn.rhs
        if (num_add_terms(lhs) <= max_num_terms):
            term_to_sub = None
            for term in lhs.as_ordered_terms():
                term_atoms = term.atoms(sympy.Symbol)
                # We want a monic term of 1 variable if we haven't found one yet
                if ((len(term_atoms) == 1) and
                    (term.as_coeff_mul()[0] == 1) and
                    (term_to_sub is None)):
                    term_to_sub = term

#                if self.invariant_interactions_on_substitution and (len(term_atoms) > 1):
#                    return

            if term_to_sub is not None:
                value = rhs - lhs + term_to_sub
                if parity(value) is not None:
                    value = parity(value)
                self.update_value(term_to_sub, value)


    def judgement_min_max(self, eqn):
        ''' If min(rhs) == max(lhs), then we know what to do

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y + z, 3)
            >>> system.judgement_min_max(eqn)
            >>> system.deductions
            {x: 1, z: 1, y: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 2*y, 5 - 2*z)
            >>> system.judgement_min_max(eqn)
            >>> system.deductions
            {x: 1, z: 1, y: 1}
        '''

        def _helper(self, lhs, rhs):
            if min_value(lhs) == max_value(rhs):
                self.set_to_min(lhs)
                self.set_to_max(rhs)

        _helper(self, eqn.lhs, eqn.rhs)
        _helper(self, eqn.rhs, eqn.lhs)

    def judgement_mini_assumption(self, eqn, num_var=4, 
                                  coef_transform=lambda x: pow(x, 0.01)):
        ''' Given an equation, assume the most common num_var are 0/1 and see
            if we can get any contradictions.
            coef_transform is a function used to rank variables by taking 
            coef_transform(abs(coef))
            
            Now also check to see if we can get any values that have to be
            equal or unequal - again by checking each case
            
            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 4*y, 5 + 2*z)
            >>> system.judgement_mini_assumption(eqn, num_var=1)
            >>> system.deductions
            {y: 1}
            
            >>> system.judgement_mini_assumption(eqn, num_var=3)
            >>> system.deductions
            {x: 1, z: 0, y: 1}
            

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2*x*y + 2*z)
            >>> system.judgement_mini_assumption(eqn, num_var=2)
            >>> system.deductions
            {x: y}
            
            >>> system.judgement_mini_assumption(eqn, num_var=3)
            >>> system.deductions
            {x: y, z: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 10*y, 5 + 2*z)
            >>> system.judgement_mini_assumption(eqn, num_var=4)
            Traceback (most recent call last):
                ...
            ContradictionException: Assumption judgement contradiction
            
            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(2*x + y, 1 + x + 4*z)
            >>> system.judgement_mini_assumption(eqn, num_var=3)
            >>> system.deductions
            {x: -y + 1, z: 0}
            
            >>> equations = ['2*z67 + 5*z68 + z78 == z56 + 4*z810 + 2*z89 + 1', 
            ...              'p4 + q4 + 2*z66 + 4*z69 == 2*z2627 + 3']
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system = EquationSolver()
            >>> for eqn in equations: system.judgement_mini_assumption(eqn, 
            ...                       num_var=4)
            >>> system.deductions
            {z69: -z66 + 1, z66: -z69 + 1, z2627: -z66 + 1, z810: z68}
            >>> system.clean_deductions()
            >>> system.solutions
            {z66: -z69 + 1, z2627: z69, z810: z68}
        '''
        # Create a dictionary of scores that we're going to use to rank variables
        var_score = defaultdict(int)
        for term, coef in (eqn.lhs + eqn.rhs).as_coefficients_dict().iteritems():
            for atom in term.atoms(sympy.Symbol):
                var_score[atom] += coef_transform(abs(coef))

        # Choose our variables
        variables = sorted(var_score.items(), key=lambda x: x[1])[-num_var:]
        variables = [v[0] for v in variables]

        # Now generate the potential values
        values = list(itertools.product((sympy.sympify(0), sympy.sympify(1)), 
                                        repeat=len(variables)))
        shuffle(values)

        # Intersection will be a set of (variable, value) tuples, where value
        # will be 0 or 1. We can then intersect with other non-contradictory
        # solutions so we are left with deductions that must be true
        intersection = None
        
        # difference_grid is a dict of (var1, var2): difference, where var1 and
        # var2 are variables and difference is the deduced difference - None
        # initially when we don't know anything about the relationship.
        # When we get contradictory relations, the tuple will be popped since
        # we can't make any deduction.
        difference_grid = dict()
        for vars_ in itertools.combinations(variables, 2):
            difference_grid[vars_] = None
        

        for vals in values:
            to_sub = dict(zip(variables, vals))
            _eqn = eqn.subs(to_sub)
            try:
                if is_equation(_eqn):
                    _eqn = standardise_equation(_eqn)
                    self.apply_contradictions([_eqn])
                
                # Process the simple 0/1 deductions
                if intersection is None:
                    intersection = set(to_sub.items())
                else:
                    intersection.intersection_update(set(to_sub.items()))
                    
                # Process the difference relations
                for key, diff in difference_grid.copy().iteritems():
                    var1, var2 = key
                    # We know they can be equal                    
                    if to_sub[var1] == to_sub[var2]:
                        # If they can also be unequal, bin it
                        if diff == 1:
                            difference_grid.pop(key)
                        else:
                            difference_grid[key] = 0
                    else:
                        if diff == 0:
                            difference_grid.pop(key)
                        else:
                            difference_grid[key] = 1
                        
                # If our intersection or difference_grid is empty, then we 
                # already know we can't deduce anything
                if (len(intersection) == 0) and (len(difference_grid) == 0):
                    return

            except ContradictionException:
                continue
        
        # If we haven't found a single solution, we must have gone somewhere
        if intersection is None:
            raise ContradictionException('Assumption judgement contradiction')
        
        # Update definite solutions
        for var, val in intersection:
            self.update_value(var, val)

        # Keep track of updated variables so we can update more gracefully
        # later
        updated = dict(intersection)
        # Sort the items so we update direct equivalences (x=y) first, and 
        # (x=1-y) relations later
        sorted_items = sorted(difference_grid.iteritems(), key=itemgetter(1))
        for (var1, var2), diff in sorted_items:
            val1 = updated.get(var1)
            val2 = updated.get(var2)
            # If we've got a relationship and a definite value, we should damn
            # well have one for the other variable
            assert (((val1 is None) and (val2 is None)) or 
                    ((val1 is not None) and (val2 is not None)))
            if val1 is not None:
                continue
            # Now update the relationship
            if diff == 1:
                # Use judgement_two_term to gracefully handle the difference
                self.judgement_two_term(sympy.Eq(var1 + var2, 1))
            elif diff == 0:
                self.update_value(var1, var2)
            else:
                # This absolutely should not happen as it should be caught
                # by checking intersection
                raise ContradictionException('This should not happen!')


    def judgement_1(self, eqn):
        ''' If x + y + z = 1 then xy = yz = zx = 0
            Generally true for any number of terms that = 1

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y + z, 1)
            >>> system.judgement_1(eqn)
            >>> system.deductions
            {x*z: 0, x*y: 0, y*z: 0}

            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + y + z + z2, 1)
            >>> system.judgement_1(eqn)
            >>> system.deductions
            {x*z: 0, z*z2: 0, y*z2: 0, x*z2: 0, y*z: 0, x*y: 0}
            
            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + 2*y + z*z2, 1)
            >>> system.judgement_1(eqn)
            >>> system.deductions
            {x*y: 0, x*z*z2: 0, y*z*z2: 0}
        '''
        if eqn.rhs != 1:
            return

        terms = eqn.lhs.as_ordered_terms()
        # We need more than 2 terms
        if len(terms) < 2:
            return

        # We also need every term to be positive
        for term in terms:
            if min_value(term) < 0:
                return

        pairs = itertools.combinations(terms, 2)
        for t1, t2 in pairs:
            self.update_value(t1*t2, 0)


#        if len(terms) == 3:
#            variables, coef = zip(*terms)
#            if all([c == 1 for c in coef]):
#                for v in itertools.permutations(variables, 2):
#                    self.update_value(v[0] * v[1], 0)

    def judgement_2_slow(self, eqn):
        ''' If a term being 1 would tip max(rhs) > max(lhs), then it must be 0

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2 + z)
            >>> system.judgement_2_slow(eqn)
            >>> system.deductions
            {z: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2 + z)
            >>> system.judgement_2_slow(eqn)
            >>> system.deductions
            {z: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 3 * y, 2 + z)
            >>> system.judgement_2_slow(eqn)
            >>> system.deductions
            {}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 5 * y, 2 + z)
            >>> system.judgement_2_slow(eqn)
            >>> system.deductions
            {y: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 1 + z)
            >>> system.judgement_2_slow(eqn)
            >>> system.deductions
            {}

            >>> system = EquationSolver()
            >>> lhs = sympy.sympify('q3 + q4')
            >>> rhs = sympy.sympify('2*q3*q4 + 2*z78 + 4*z79')
            >>> eqn = sympy.Eq(lhs, rhs)
            >>> system.judgement_2_slow(eqn)
            >>> system.deductions
            {z79: 0}
            
            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x+y, 2*z + 1)
            >>> system.judgement_2_slow(eqn)
            >>> system.deductions
            {z: 0}
        '''
        def _helper(lhs, rhs):
            rhs_term_dict = rhs.as_coefficients_dict()
#            top_coef = max(map(abs, rhs_term_dict.values()))
            for term, coef in rhs_term_dict.iteritems():
#                if abs(coef) != top_coef:
#                    continue
                to_subs = {term: 1}                
                _rhs = rhs.subs(to_subs)
                _lhs = lhs.subs(to_subs)
                if max_value(_lhs) < min_value(_rhs):
                    self.update_value(term, 0)

        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)

    def judgement_2(self, eqn):
        ''' If a term being 1 would tip max(rhs) > max(lhs), then it must be 0

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2 + z)
            >>> system.judgement_2(eqn)
            >>> system.deductions
            {z: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2 + z)
            >>> system.judgement_2(eqn)
            >>> system.deductions
            {z: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 3 * y, 2 + z)
            >>> system.judgement_2(eqn)
            >>> system.deductions
            {}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 5 * y, 2 + z)
            >>> system.judgement_2(eqn)
            >>> system.deductions
            {y: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 1 + z)
            >>> system.judgement_2(eqn)
            >>> system.deductions
            {}

            >>> system = EquationSolver()
            >>> lhs = sympy.sympify('q3 + q4')
            >>> rhs = sympy.sympify('2*q3*q4 + 2*z78 + 4*z79')
            >>> eqn = sympy.Eq(lhs, rhs)
            >>> system.judgement_2(eqn)
            >>> system.deductions
            {z79: 0}
            
            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x+y, 2*z + 1)
            >>> system.judgement_2(eqn)
            >>> system.deductions
            {z: 0}
        '''
        def _helper(lhs, rhs):
            lhs_max = max_value(lhs)

            rhs_terms = rhs.as_ordered_terms()
            if is_constant(rhs_terms[-1]):
                rhs_const = rhs_terms.pop()
            else:
                rhs_const = 0

            for term in rhs_terms:
                if max_value(term) + rhs_const > lhs_max:
                    self.set_to_min(term)

        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)


    def _judgement_2_old(self, eqn):
        ''' If max(lhs) < max(rhs) and there is only 1 variable term
            in the RHS, then this term must be 0

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x+y, 2*z + 1)
            >>> system._judgement_2_old(eqn)
            >>> system.deductions
            {z: 0}
        '''
        lhs, rhs = eqn.lhs, eqn.rhs
        if (max_value(lhs) < max_value(rhs)):
            terms = rhs.as_ordered_terms()
            if (len(terms) == 2) and terms[-1].is_constant():
                self.set_to_min(rhs)

    def judgement_3(self, eqn):
        ''' If max(lhs) = max(rhs) and rhs is constant, then every term on the
            left is 1.
            Similarly for minimum

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2)
            >>> system.judgement_3(eqn)
            >>> system.deductions
            {x: 1, y: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y + z, 0)
            >>> system.judgement_3(eqn)
            >>> system.deductions
            {x: 0, z: 0, y: 0}
        '''
        def _helper(lhs, rhs):
            if not is_constant(rhs):
                return
            elif (max_value(lhs) == max_value(rhs)):
                self.set_to_max(lhs)
            elif (min_value(lhs) == min_value(rhs)):
                self.set_to_min(lhs)

        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)


    def judgement_4(self, eqn):
        ''' If min(LHS) > min(RHS) and we only have one variable term on the
            RHS, this must be 1

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 1, 2*z)
            >>> system.judgement_4(eqn)
            >>> system.deductions
            {z: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 1, z)
            >>> system.judgement_4(eqn)
            >>> system.deductions
            {z: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x * y + 1, 2 * z)
            >>> system.judgement_4(eqn)
            >>> system.deductions
            {z: 1}

            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + 3 * y + 1, 2 * z + z2)
            >>> system.judgement_4(eqn)
            >>> system.deductions
            {}
            
            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 3 * y + 4, 5 * z)
            >>> system.judgement_4(eqn)
            >>> system.deductions
            {z: 1}
            
            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 3 * y + 4, 5)
            >>> system.judgement_4(eqn)
            >>> system.deductions
            {}
        '''
        def _helper(lhs, rhs):
            if (min_value(lhs) > min_value(rhs)):
                if (num_add_terms(rhs) == 1) and (not is_constant(rhs)):
                        self.update_value(rhs.as_coeff_Mul()[1], 1)
        
        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)

    def judgement_5(self, eqn, increase_complexity=False, 
                    invariant_interactions_on_substitution=True):
        ''' Parity argument used when the RHS is always even.
        
            If we have 1 odd term in the LHS, it must be 0.            
            If we have 2 odd terms in the LHS, they must be equal.            
            x + 2y + z = 4z2 -> x = z

            If we have 3 odd terms and we are allowed to increase complexity,
            then we can say something more:
            x + y + z = 0 mod 2
            =>
            x = (y-z)^2 = y + z - 2*y*z
            
            Also works with any even RHS and any number of even variables on
            the LHS.

            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + 2*y + z, 2*z2)
            >>> system.judgement_5(eqn)
            >>> system.deductions
            {x: z}

            >>> system = EquationSolver()
            >>> x, z1, z2, z3, z4 = sympy.symbols('x z1 z2 z3 z4')
            >>> eqn = sympy.Eq(x + 2*z1 + 2*z2, 4*z3*z4 + 2)
            >>> system.judgement_5(eqn)
            >>> system.deductions
            {x: 0}

            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + 2*y + 3*z, 2 * z2)
            >>> system.judgement_5(eqn)
            >>> system.deductions
            {x: z}

            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + 2*y + 4*z, 2*z2)
            >>> system.judgement_5(eqn)
            >>> system.deductions
            {x: 0}

            >>> system = EquationSolver()
            >>> x, y, z, u, v = sympy.symbols('x y z u v')
            >>> eqn = sympy.Eq(x + y + z + 2*u, 4*v)
            >>> system.judgement_5(eqn, increase_complexity=True, 
            ... invariant_interactions_on_substitution=False)
            >>> system.deductions
            {x: -2*y*z + y + z}

            >>> eqns = ['q2 + q3 + 2*z4950 + 1 == 2*q2*q3 + 2*z56 + 4*z57']
            >>> eqn = str_eqns_to_sympy_eqns(eqns)[0]
            >>> system = EquationSolver()
            >>> system.judgement_5(eqn, increase_complexity=True,
            ... invariant_interactions_on_substitution=True)
            >>> system.deductions
            {q2: -q3 + 1}
        '''

        def _helper(lhs, rhs):
            if parity(rhs) != 0:
                return

            odd_terms = []
            const = 0
            for term, term_coef in lhs.as_coefficients_dict().iteritems():
                if (term_coef % 2):
                    if term == 1:
                        const = 1
                    odd_terms.append(term)

            if len(odd_terms) == 1:
                self.update_value(odd_terms.pop(), 0)

            elif len(odd_terms) == 2:
                if is_constant(odd_terms[0]):
                    self.update_value(odd_terms[1], odd_terms[0])
                else:
                    self.update_value(*odd_terms)
            
            elif (increase_complexity and (len(odd_terms) == 3)):
                if const or (not invariant_interactions_on_substitution):
                    x, y, z = odd_terms
                    # If we get a 1, permute the variables so it's nicer to play with
                    if x == 1:
                        y, z, x = x, y, z
                    self.update_value(x, y + z - 2*y*z)

        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)

    def judgement_6(self, eqn, increase_complexity=False,
                    invariant_interactions_on_substitution=True):
        ''' Parity argument used if RHS is always odd.
                    
            If we have 1 odd term in the LHS, it must be 1.
            If we have 2 odd terms in the LHS, they must sum to 1.            
            x + 2y + z = 4z2 + 1 -> x = 1 - z
            
            If we have 3 odd terms and we are allowed to increase complexity,
            then we can say something more:
            x + y + z = 1 mod 2
            =>
            x = 1 - (y-z)^2 = 1 - y - z + 2yz

            If we have 
            and the LHS has 2 monic terms
            then the sum must be 1. If it has 1 monic term, it must be 1
            Also make sure we don't replicate judgement_two_term

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y + 2*x*y, 2 * z + 1)
            >>> system.judgement_6(eqn, increase_complexity=True)
            >>> system.deductions
            {x: -y + 1}

            >>> system = EquationSolver()
            >>> x, y = sympy.symbols('x y')
            >>> eqn = sympy.Eq(x + y, 1)
            >>> system.judgement_6(eqn, increase_complexity=False)
            >>> system.deductions
            {}

            >>> system = EquationSolver()
            >>> x, y = sympy.symbols('x y')
            >>> eqn = sympy.Eq(x + y, 1)
            >>> system.judgement_6(eqn, increase_complexity=True)
            >>> system.deductions
            {x: -y + 1}
            
            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y + 2*z, 3)
            >>> system.judgement_6(eqn, increase_complexity=True)
            >>> system.deductions
            {x: -y + 1}
            
            >>> system = EquationSolver()
            >>> eqns = ['2*q5*z1213 + 2*q6*q7 + 1 == z1213 + 2*z89']
            >>> eqn = str_eqns_to_sympy_eqns(eqns)[0]
            >>> system.judgement_6(eqn)
            >>> system.deductions
            {z1213: 1}

            >>> system = EquationSolver()
            >>> x, y, z, u, v = sympy.symbols('x y z u v')
            >>> eqn = sympy.Eq(x + y + z + 2*u, 4*v + 1)
            >>> system.judgement_6(eqn, increase_complexity=True,
            ... invariant_interactions_on_substitution=False)
            >>> system.deductions
            {x: 2*y*z - y - z + 1}
        '''
        def _helper(lhs, rhs):        
            if parity(rhs) != 1:
                return
                
            odd_terms = []
            const = 0
            for term, term_coef in lhs.as_coefficients_dict().iteritems():
                if (term_coef % 2):
                    if term == 1:
                        const = 1
                    odd_terms.append(term)

            if len(odd_terms) == 1:
                self.update_value(odd_terms.pop(), 1)
            elif (len(odd_terms) == 2):
                if const:
                    self.update_value(odd_terms.pop(), 0)
                elif increase_complexity:
                    self.judgement_two_term(sympy.Eq(sum(odd_terms), 1))
            elif (increase_complexity and 
                 ((not invariant_interactions_on_substitution) or const) and
             (len(odd_terms) == 3)):
                x, y, z = odd_terms
                # Make another variable 1, so that it plays nicer
                if x == 1:
                    y, z, x = x, y, z
                self.update_value(x, 1 - y - z + 2*y*z)
        
        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)


    def judgement_7(self, eqn):
        ''' Special case of judgement_5
            x + y = 2z -> x = y = z

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2*z)
            >>> system.judgement_7(eqn)
            >>> system.deductions
            {z: x, y: x}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x * y + x, 2 * z)
            >>> system.judgement_7(eqn)
            >>> system.deductions
            {x*y: x, z: x*y}
        '''
        def _helper(lhs, rhs):
            if len(lhs.as_ordered_terms()) == 2:
                t1, t2 = lhs.as_ordered_terms()
                if ((t1.as_coeff_mul()[0] == 1) and
                    (not t1.is_constant()) and
                    (t2.as_coeff_mul()[0] == 1) and
                    (not t2.is_constant()) and
                    (rhs.as_coeff_mul()[0] == 2)):
                    self.update_value(t2, t1)
                    self.update_value(rhs / 2, t1)
        
        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)

    def judgement_8(self, eqn):
        ''' If a term being 0 would tip max(rhs) < min(lhs), then it must be 1

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 3, 2 + z)
            >>> system.judgement_8(eqn)
            >>> system.deductions
            {z: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(10, 9*x + y)
            >>> system.judgement_8(eqn)
            >>> system.deductions
            {x: 1, y: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(9, 9*x + y)
            >>> system.judgement_8(eqn)
            >>> system.deductions
            {x: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(4*x*y + z, 3)
            >>> system.judgement_8(eqn)
            >>> system.deductions
            {x*y: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(4*x, 3*z + y)
            >>> system.judgement_8(eqn)
            >>> system.deductions
            {}

            >>> system = EquationSolver()
            >>> lhs = sympy.sympify('x + y + 2')
            >>> rhs = sympy.sympify('10*x*y + 2*x')
            >>> eqn = sympy.Eq(lhs, rhs)
            >>> system.judgement_8(eqn)
            >>> system.deductions
            {}
        '''
        def _helper(lhs, rhs):
            lhs_min = min_value(lhs)
            rhs_max = max_value(rhs)

            for term in rhs.as_ordered_terms():
                if rhs_max - max_value(term) < lhs_min:
                    self.set_to_max(term)

        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)

    def judgement_8_slow(self, eqn):
        ''' If a term being 0 would tip max(rhs) < min(lhs), then it must be 1

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 3, 2 + z)
            >>> system.judgement_8_slow(eqn)
            >>> system.deductions
            {z: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(10, 9*x + y)
            >>> system.judgement_8_slow(eqn)
            >>> system.deductions
            {x: 1, y: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(9, 9*x + y)
            >>> system.judgement_8_slow(eqn)
            >>> system.deductions
            {x: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(4*x*y + z, 3)
            >>> system.judgement_8_slow(eqn)
            >>> system.deductions
            {x*y: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(4*x, 3*z + y)
            >>> system.judgement_2(eqn)
            >>> system.deductions
            {}

            >>> system = EquationSolver()
            >>> lhs = sympy.sympify('x + y + 2')
            >>> rhs = sympy.sympify('10*x*y + 2*x')
            >>> eqn = sympy.Eq(lhs, rhs)
            >>> system.judgement_8_slow(eqn)
            >>> system.deductions
            {x: 1}
        '''
        def _helper(lhs, rhs):
            rhs_term_dict = rhs.as_coefficients_dict()
#            top_coef = max(map(abs, rhs_term_dict.values()))
            for term, coef in rhs_term_dict.iteritems():
#                if abs(coef) != top_coef:
#                    continue
                to_subs = {term: 0}                
                _rhs = rhs.subs(to_subs)
                _lhs = lhs.subs(to_subs)
                if max_value(_rhs) < min_value(_lhs):
                    self.update_value(term, 1)
                
                
        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)


    def judgement_9(self, eqn, increase_complexity=False, 
                    invariant_interactions_on_substitution=True):
        ''' Parity argument for when we have 1 variable that determines parity
            on the LHS.
            
            If we have 1 parity-determining variable on the RHS, then preserve
            parity.
            
            If we have 2 parity-determining variables on the RHS, then also
            preserve parity
            
            >>> eqns = ['3*q5 + 2*z89 == 4*q5*z89 + 5*q7 + 1',
            ...         '3*x + 2*a1 == 4*a2*a3 + y + z + 1']
            >>> eqns = str_eqns_to_sympy_eqns(eqns)
            >>> system = EquationSolver()
            >>> for eqn in eqns: system.judgement_9(eqn, increase_complexity=True)
            >>> system.deductions
            {q5: -q7 + 1}

            >>> for eqn in eqns: system.judgement_9(eqn, increase_complexity=True,
            ... invariant_interactions_on_substitution=False)
            >>> system.deductions
            {q5: -q7 + 1, x: 2*y*z - y - z + 1}
        '''
        if num_add_terms(eqn.lhs) == num_add_terms(eqn.rhs) == 1:
            return

        def _helper(lhs, rhs):
            # First find the 1 odd variable on the LHS
            odd_terms = []
            for term, term_coef in lhs.as_coefficients_dict().iteritems():
                if (term_coef % 2):
                    if term == 1:
                        return
                    odd_terms.append(term)
            if len(odd_terms) != 1:
                return
            
            lhs_determining_var = odd_terms.pop()
            
            # Now if we can say something about parity, do it
            odd_terms = []
            for term, term_coef in rhs.as_coefficients_dict().iteritems():
                if (term_coef % 2) and (term != 1):
                    odd_terms.append(term_coef * term)
            
            if len(odd_terms) == 1:
                rhs_determining_var = odd_terms.pop()
                rhs_parity = parity(rhs - rhs_determining_var)
                rhs_determining_var = rhs_determining_var.as_coeff_Mul()[1]
                if rhs_parity == 0:
                    self.update_value(lhs_determining_var, rhs_determining_var)
                elif (rhs_parity == 1) and (increase_complexity):
                    self.update_value(lhs_determining_var, 1 - rhs_determining_var)

            elif len(odd_terms) == 2:
                # These judgements only increace complexity
                if not increase_complexity:
                    return

                x, y = odd_terms
                rhs_parity = parity(rhs - x - y)
                x = x.as_coeff_Mul()[1]
                y = y.as_coeff_Mul()[1]
                if rhs_parity == 0:
                    self.judgement_two_term(sympy.Eq(lhs_determining_var, x + y))
                elif (rhs_parity == 1) and (not invariant_interactions_on_substitution):
                    self.update_value(lhs_determining_var, 1 - x - y + 2*x*y)
        
        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)

        
    def _judgement_10i(self, eqn):
        ''' x + y + z + 2a = 2b -> xy + xz = yx + yz = zx + zy = 0
            Also works with any even RHS and any number of even variables on
            the LHS.

            >>> system = EquationSolver()
            >>> x, y, z, z2, z3 = sympy.symbols('x y z z2, z3')
            >>> eqn = sympy.Eq(x + 2*y + z + z2, 2*z3)
            >>> system._judgement_10i(eqn)
            >>> system.deductions
            {x*z2 + z*z2: 0, x*z + x*z2: 0, x*z + z*z2: 0}

            >>> system = EquationSolver()
            >>> eqn = sympy.Eq(x + 2*y + 3*z + 7 * z2, 4*z3)
            >>> system._judgement_10i(eqn)
            >>> system.deductions
            {x*z2 + z*z2: 0, x*z + x*z2: 0, x*z + z*z2: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y + 1, 2 * z)
            >>> system._judgement_10i(eqn)
            >>> system.deductions
            {}
        '''

        def _helper(lhs, rhs):
            for term in rhs.as_ordered_terms():
                if term.as_coeff_mul()[0] % 2:
                    return

            odd_terms = []
            for term in lhs.as_ordered_terms():
                term_coef, term = term.as_coeff_Mul()
                if (term_coef % 2):
                    if len(odd_terms) > 3:
                        return
                    else:
                        if term.is_constant():
                            return
                        odd_terms.append(term)

            if len(odd_terms) == 3:
                pairs = [a * b for a, b in itertools.combinations(odd_terms, 2)]
                expressions = [ab + bc for ab, bc in itertools.combinations(pairs, 2)]
                for expr in expressions:
                    self.update_value(expr, 0)


        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)




    ## Look for contradictions
    def contradiction_1(self, eqn):
        ''' Check the values could be equal 
        
        >>> x, y, z = sympy.symbols('x y z')

        >>> eqn = sympy.Eq(x*y*z, 2)
        >>> system = EquationSolver()
        >>> system.contradiction_1(eqn)
        Traceback (most recent call last):
            ...
        ContradictionException: contradiction_1: x*y*z == 2

        >>> eqn = sympy.Eq(x*y*z)
        >>> system = EquationSolver(equations=[eqn])
        >>> system.solve_equations()
        >>> system.solutions
        {}
        

        >>> eqn = sympy.Eq(x*y*z + 1)
        >>> system = EquationSolver()
        >>> system.contradiction_1(eqn)
        Traceback (most recent call last):
            ...
        ContradictionException: contradiction_1: x*y*z + 1 == 0
        '''
        def _helper(self, lhs, rhs):
            if min_value(lhs) > max_value(rhs):
                raise ContradictionException('contradiction_1: {}'.format(eqn))

        _helper(self, eqn.lhs, eqn.rhs)
        _helper(self, eqn.rhs, eqn.lhs)
    
    def contradiction_2(self, eqn):
        ''' Check the parity 
        
        >>> x, y, z = sympy.symbols('x y z')
        
        >>> eqn = sympy.Eq(2*x*y + 4*z + 1)
        >>> system = EquationSolver()
        >>> system.contradiction_2(eqn)
        Traceback (most recent call last):
            ...
        ContradictionException: contradiction_2: 2*x*y + 4*z + 1 == 0
        
        >>> eqn = sympy.Eq(2*x*y * 4*z - 2)
        >>> system = EquationSolver()
        >>> system.contradiction_2(eqn)

        >>> eqn = sympy.Eq(2*x*y + z - 1)
        >>> system = EquationSolver()
        >>> system.contradiction_2(eqn)
        '''
        l_parity = parity(eqn.lhs)
        if l_parity is not None:
            r_parity = parity(eqn.rhs)
            if (r_parity is not None) and (l_parity != r_parity):
                raise ContradictionException('contradiction_2: {}'.format(eqn))


## Simple chunker for partitioning lists
def _batcher(iterable, batch_size, fill_value=None):
    ''' Given a list of things, return a partition
        >>> list(_batcher(range(20), 5))
        [(0, 1, 2, 3, 4), (5, 6, 7, 8, 9), (10, 11, 12, 13, 14), (15, 16, 17, 18, 19)]
        >>> list(_batcher(range(20), 4))
        [(0, 1, 2, 3), (4, 5, 6, 7), (8, 9, 10, 11), (12, 13, 14, 15), (16, 17, 18, 19)]
        >>> list(_batcher(range(20), 3))
        [(0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11), (12, 13, 14), (15, 16, 17), (18, 19, None)]
        >>> list(_batcher(range(10), 3, fill_value=100))
        [(0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 100, 100)]
    '''
    return itertools.izip_longest(*[iter(iterable)] * batch_size, 
                                    fillvalue=fill_value)

## Conflict resolution
def _simplest(expr1, expr2):
    ''' Return the simplest of expr1 and expr2, giving precedence to expr1.
        Used by the solver when we have 2 different symbolic values and we
        want to determine which to assign.
        ASSUMES THE FIRST ARGUMENT IS THE OLD VALUE
    '''
    if len(expr2.atoms()) > len(expr1.atoms()):
        return expr2
    else:
        return expr1

# Inspection
def _get_judgement():
    ''' Find the judgement calling update_value. Horrible hackery, but at least
        it's in 1 place
    '''
    up_again = ['set_to_min', 'set_to_max', 'update_value', '_helper',
                'judgement_two_term']

    ind = 3
    frame = inspect.stack()[ind]
    caller_name = frame[-3]
    while caller_name in up_again:
        ind += 1
        frame = inspect.stack()[ind]
        caller_name = frame[-3]
        if ind > 100:
            break

    eqn = frame[0].f_locals.get('eqn')
    return caller_name, eqn


if __name__ == "__main__":
    import doctest
    doctest.testmod()
