# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 22:14:25 2015

@author: Richard
"""

"""
Created on Fri Dec 26 12:35:16 2014

Solve a system of equations with binary variables

@author: Richard Tanburn
"""
from copy import deepcopy
import itertools
import pickle
import sympy
from sympy.core.cache import clear_cache

from contradiction_exception import ContradictionException
from equivalence_dict import BinaryEquivalenceDict
from sympy_helper_fns import (max_value, min_value, is_equation,
                              remove_binary_squares_eqn, balance_terms,
                              cancel_constant_factor, is_constant,
                              num_add_terms, parity, is_monic, is_one_or_zero,
                              remove_binary_squares, expressions_to_variables,
                              gather_monic_terms, square_equations,
                              str_eqns_to_sympy_eqns, standardise_equation,
                              is_simple_binary, dict_as_eqns)
from objective_function_helper import (equations_to_vanilla_coef_str, 
                                       equations_to_vanilla_objective_function)

from sympy_paralleliser import paralellised_subs, get_pool


__author__ = "Richard Tanburn"
__credits__ = ["Richard Tanburn", "Nikesh Dattani"]
__version__ = "0.0.1"
__status__ = "Prototype"

class SolverBase(object):
    ''' Solver of equations. Enforce certain methods to allow reusability '''

    SOLUTIONS_TYPE = dict

    def __init__(self, equations=None, variables=None,
                 output_filename=None, parallelise=False, *args, **kwargs):
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

        # Solutions. Subset of deductions, where the key is a single variable.
        self.solutions = self.SOLUTIONS_TYPE()

        # File to print to
        self.output_filename = output_filename
        self._file = None

        # Allow parallelisation
        self.parallelise = parallelise
        self._pool = None
        
        # Determine whether we have fixed which way around p and q will go
        self._fix_pq_soln_used = False

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
        cls = type(self)
        copy = cls(deepcopy(self.equations), 
                   deepcopy(self.variables),
                   output_filename=self.output_filename,
                   parallelise=self.parallelise)
                              
        # Now use deepcopy to copy everything else
        copy.num_qubits_start = self.num_qubits_start
        copy.solutions = deepcopy(self.solutions)
        copy._fix_pq_soln_used = self._fix_pq_soln_used
        return copy

    # Pickling
    # This isn't mature/finished, but manages to write equations, deductions and
    # solutions to disk
    def __getstate__(self):
        return (self.equations, self.solutions, self.variables, self.num_qubits_start,
                self.output_filename, self.parallelise, self._fix_pq_soln_used)
        
    def __setstate__(self, state):
        (self.equations, self.solutions, 
         self.variables, self.num_qubits_start,
         self.output_filename, self.parallelise, 
         self._fix_pq_soln_used) = state
         
    def to_disk(self, filename):
        ''' Write a state to disk '''
        if filename is None:
            return
        pickle.dump(self, open(filename, 'w'))
    
    @staticmethod
    def from_disk(filename, **kwargs):
        ''' Load from disk '''
        if filename is None:
            raise ValueError('Cannot load from filename None')
        data = pickle.load(open(filename, 'r'))
        return data

    @property
    def solutions_as_equations(self):
        ''' Return solutions as a list of equations '''
        return dict_as_eqns(self.solutions)

    @property    
    def _length_tuple(self):
        ''' Return a tuple of the lengths of equations, deductions, solutions 
        '''
        return len(self.equations), len(self.solutions)

    def solve_equations(self, max_iter=250, verbose=False):
        ''' Solve a system of equations
        '''
        raise NotImplementedError('solve_equations not implemented')

    def add_solution(self, variable, value):
        ''' Add a solution to the solution dict '''
        assert is_simple_binary(value)
        assert isinstance(variable, sympy.Symbol)

        current_sol = self.solutions.get(variable)
        if ((current_sol is not None) and is_constant(value) and 
            is_constant(current_sol) and value != current_sol):
            raise ContradictionException('Contradiction in add_solution()')
        self.solutions[variable] = value

    @property
    def final_equations(self):
        ''' final_equations are the final equations that the solver just
            couldn't handle, and would like to pass on to someone else
        '''
        raise NotImplementedError('final_equations not implemented')

    def print_summary(self):
        ''' Print a summary of the information held in the object '''
        unsolved_var = self.unsolved_var

        self.print_('Final Variables')
        self.print_(unsolved_var)


        self.print_('Final Equations')
        for e in self.final_equations:
            self.print_(e)

        self.print_('Num Qubits Start: {}'.format(self.num_qubits_start))
        self.print_('Num Qubits End: {}'.format(len(unsolved_var)), close=True)

    @property
    def unsolved_var(self):
        ''' Return a set of variables we haven't managed to eliminate '''
        #TODO Think of a cleverer way to do this with solutions that still
        # hold some information, so aren't 'simple' solutions
        return set(self.variables.values()).difference(self.solutions.keys())

    @property
    def objective_function(self):
        ''' Return the final objective function, using self.final_equations '''
        return equations_to_vanilla_objective_function(self.final_equations)

    def objective_function_to_file(self, filename=None):
        ''' Write the objective function to a file, or printing if None.
            Also include the dictionary of variable number to original variable
        '''
        out = equations_to_vanilla_coef_str(self.final_equations)
        #out = equations_to_auxillary_coef_str(self.equations)

        if filename is None:
            self.print_(out)
        else:
            f = open(filename, 'a')
            f.write(out)
            f.close()

    def open_pool(self):
        ''' Open the pool, if it isn't already '''
        if self._pool is not None:
            self._pool = get_pool()

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
            >>> system = SolverBase()
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
                self.open_pool()
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

    def get_var(self, var):
        ''' Return symbolic variable

            >>> system = SolverBase()
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

    def fix_pq_soln(self, solutions=None):
        ''' Look through a dictionary of solutions for p_i, q_i such that 
            p_i = 1 - p_i, and both are symbolic. Then we can fix them by
            symmetry of factors!

            >>> p1, p2, p3, q1, q2, q3, z1, z2, z3 = sympy.symbols('p1, p2, p3, q1, q2, q3, z1, z2, z3')            
            >>> system = SolverBase()
            
            >>> system.fix_pq_soln({p1: 0, q2: 1})
            >>> system.solutions
            {}

            >>> system.fix_pq_soln({p1: z1, q2: z2})
            >>> system.solutions
            {}

            >>> system.fix_pq_soln({p1: p2, q2: 1 - q3})
            >>> system.solutions
            {}
            
            >>> system.fix_pq_soln({p1: 1 - z2, q2: 1 - z3})
            >>> system.solutions
            {}

            >>> system.fix_pq_soln({p1: 1 - z1, q2: z1})
            >>> system.solutions
            {}

            >>> system.fix_pq_soln({p2: 1 - q3})
            >>> system.solutions
            {}

            >>> system.fix_pq_soln({p2: q3})
            >>> system.solutions
            {}

            >>> system.fix_pq_soln({p1: q1})
            >>> system.solutions
            {}

            >>> system = SolverBase()
            >>> system.fix_pq_soln({p1: 1 - q1})
            >>> system.solutions
            {p1: 1, q1: 0}

            This needs fixing!
            >>> system = SolverBase()
            >>> system.fix_pq_soln({p1: z2, q1: 1 - z2})
            >>> system.solutions
            {}
        '''
        # We can only use this judgement once
        if self._fix_pq_soln_used:
            return

        if solutions is None:
            solutions = self.solutions
        
        for k, v in solutions.iteritems():
            strk = str(k)
            strv = str(v)
            if strk[0] == 'p':
                num = strk[1:]
                if strv == '-q{} + 1'.format(num):
                    self._fix_pq_soln_used = True
                    self.add_solution(k, 1)
                    self.add_solution(1 - v, 0)
                    return
            if strk[0] == 'q':
                num = strk[1:]
                if strv == '-p{} + 1'.format(num):
                    self._fix_pq_soln_used = True
                    self.add_solution(k, 0)
                    self.add_solution(1 - v, 1)
                    return

class BinarySolutionSolverBase(SolverBase):
    ''' A class that is to be used if the funky BinaryEquivalenceDict is wanted
        to hold solutions
    '''
    SOLUTIONS_TYPE = BinaryEquivalenceDict

    @property
    def unsolved_var(self):
        ''' Because of the funny way we store variables in a 
            BinaryEquivalenceDict, count the number of variables that have
            a non-constant root
            
            >>> variables = sympy.symbols('x, y, z')
            >>> x, y, z = variables            
            >>> solver = BinarySolutionSolverBase(variables={str(v): v for v in variables})
            >>> solver.unsolved_var
            set([x, z, y])
            >>> solver.add_solution(x, 1)
            >>> solver.unsolved_var
            set([z, y])
            >>> solver.add_solution(y, 1 - z)
            >>> solver.unsolved_var
            set([z])
            >>> solver.add_solution(y, 1)
            >>> solver.unsolved_var
            set([])
        '''
        sols = [(v, self.solutions[v]) for v in self.variables.values()]
        unsolved = set([var for var, sol in sols if (var.atoms(sympy.Symbol) == 
                                                     sol.atoms(sympy.Symbol))])
        return unsolved

def unique_array_stable(array):
    ''' Given a list of things, return a new list with unique elements with
        original order preserved (by first occurence)
        
        >>> print unique_array_stable([1, 3, 5, 4, 7, 4, 2, 1, 9])
        [1, 3, 5, 4, 7, 2, 9]
    '''
    seen = set()
    seen_add = seen.add
    return [x for x in array if not (x in seen or seen_add(x))]

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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
