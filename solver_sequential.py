# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 00:05:28 2015

@author: Richard
"""
import sys
sys.setrecursionlimit(100)

import itertools
from operator import itemgetter

import sympy

from contradictions import apply_contradictions
from contradiction_exception import ContradictionException
from equivalence_dict import BinaryEquivalenceDict
from solver_base import BinarySolutionSolverBase
from sympy_helper_fns import is_equation, balance_terms, is_constant, expressions_to_variables
from sympy.core.cache import clear_cache


class SolverSequential(BinarySolutionSolverBase):
    ''' Class that represents a search through the solution space '''

    def __init__(self, *args, **kwargs):
        super(SolverSequential, self).__init__(*args, **kwargs)

        # List of (state_dict, evaluated_eqns). None if no equations have been
        # added
        self.valid_states = None
        
        # Queue of variables to substitute
        self.variables_to_sub = []
        
        # Track of variables already subbed
        self.variables_subbed = set()    

        # Now we've initialised everything, add any equations we might have
        # been passed in the initialisation
        for eqn in self.equations:
            self.add_equation(eqn, _init=True)
    
    def reorder_variables_to_sub(self):
        ''' Reorder the queue 
        
            >>> variables = sympy.symbols('x y z p1 p2 p10')            
            >>> search = SolverSequential()
            >>> search.variables_to_sub = variables
            >>> search.reorder_variables_to_sub()
            >>> search.variables_to_sub
            [p1, p10, p2, x, y, z]
        '''
        self.variables_to_sub = sorted(self.variables_to_sub, key=str)
    
    def check_states(self, check_state_vars_subbed=False):
        ''' Check we have somewhere to keep going! 
        
            >>> search = SolverSequential()
            >>> search.check_states()
            
            >>> search.valid_states = []
            >>> search.check_states()
            Traceback (most recent call last):
                ...
            ContradictionException: No more valid states!
            
            >>> x, y = sympy.symbols('x y')            
            >>> search.valid_states = [{x: 1}, []]
            >>> search.check_states()
            
            >>> search = SolverSequential()
            >>> search.add_equation(sympy.Eq(x, 1))
            >>> search.variables_subbed.add(y)
            >>> search.check_states()
            >>> search.check_states(check_state_vars_subbed=True)
            Traceback (most recent call last):
                ...
            AssertionError
        '''
        if self.valid_states is None:
            return
        elif len(self.valid_states) == 0:
            raise ContradictionException('No more valid states!')

        if check_state_vars_subbed:
            for state in self.valid_states:
                subbed = state[0]                
                assert set(subbed.keys()) == self.variables_subbed

    def _pop_variables_from_queue(self, vars_to_sub):
        ''' Given a list of variables, add them to the set of subbed variables
            and remove them from the queue
            
            >>> u, v, x, y, z = sympy.symbols('u v x y z')
            >>> eqns = ['x + y == 1', 'x*y + u + v == 2', 'u + z == 1']
            >>> eqns = str_eqns_to_sympy_eqns(eqns)
            >>> solver = SolverSequential(eqns)
            >>> solver.variables_to_sub
            [x, y, u, v, z]
            
            >>> solver.sub_var(1)
            >>> solver.variables_to_sub
            [y, u, v, z]
            
            >>> solver.sub_var(v)
            >>> solver.variables_to_sub
            [y, u, z]
            
            >>> solver.sub_var([y, z])
            >>> solver.variables_to_sub
            [u]
            
            >>> for s in solver.valid_states: print s
            ({v: 1, x: 0, z: 0, y: 1}, [u == 1, u == 1])
            ({v: 1, x: 0, z: 1, y: 1}, [u == 1, u == 0])
            ({v: 1, x: 1, z: 0, y: 0}, [u == 1, u == 1])
            ({v: 1, x: 1, z: 1, y: 0}, [u == 1, u == 0])
            >>> solver.sub_var(v)
            >>> for s in solver.valid_states: print s
            ({v: 1, x: 0, z: 0, y: 1}, [u == 1, u == 1])
            ({v: 1, x: 0, z: 1, y: 1}, [u == 1, u == 0])
            ({v: 1, x: 1, z: 0, y: 0}, [u == 1, u == 1])
            ({v: 1, x: 1, z: 1, y: 0}, [u == 1, u == 0])
            >>> solver.variables_to_sub
            [u]
            
        '''
        set_vars_to_sub = set(vars_to_sub)
        
        self.variables_subbed.update(set_vars_to_sub)
        
        # Now remove them from the queue
        filter_ = lambda x: x not in set_vars_to_sub
        self.variables_to_sub = filter(filter_, self.variables_to_sub)

    def sub_var(self, vars_to_sub=5, max_states=10000):
        ''' Substitute more variables into each system. This can be:
            None - substitute all variables in the queue
            int - Take the first vars_to_sub from the queue
            iterable of variables - do these explicitly
            a variable - do just this one
            
            Break the process if we get more than max_states valid states

            Walk slowly through an example
            >>> eqns = ['x + y == 1', 'x*y + u + v == 2', 'u + z == 1']
            >>> eqns = str_eqns_to_sympy_eqns(eqns)
            >>> solver = SolverSequential(eqns)
            >>> solver.variables_to_sub
            [x, y, u, v, z]
            >>> solver.valid_states
            [({}, [x + y == 1, u + v + x*y == 2, u + z == 1])]
            
            Sub in 1 variable
            >>> solver.sub_var(1)
            >>> solver.variables_subbed
            set([x])
            >>> solver.variables_to_sub
            [y, u, v, z]
            >>> for s in solver.valid_states: print s
            ({x: 0}, [y == 1, u + v == 2, u + z == 1])
            ({x: 1}, [y == 0, u + v + y == 2, u + z == 1])

            Sub in 1 more
            >>> solver.sub_var(1)
            >>> solver.variables_subbed
            set([x, y])
            >>> solver.variables_to_sub
            [u, v, z]
            >>> for s in solver.valid_states: print s
            ({x: 0, y: 1}, [u + v == 2, u + z == 1])
            ({x: 1, y: 0}, [u + v == 2, u + z == 1])
            
            Now sub in 2 at one time
            >>> solver.sub_var(2)
            >>> solver.variables_subbed
            set([v, u, x, y])
            >>> solver.variables_to_sub
            [z]
            >>> for s in solver.valid_states: print s
            ({v: 1, u: 1, x: 0, y: 1}, [z == 0])
            ({v: 1, u: 1, x: 1, y: 0}, [z == 0])

            But we could have just done it all in one go!
            >>> eqns = ['x + y == 1', 'x*y + u + v == 2', 'u + z == 1']
            >>> eqns = str_eqns_to_sympy_eqns(eqns)
            >>> solver = SolverSequential(eqns)
            >>> solver.sub_var(None)
            >>> solver.variables_subbed
            set([v, u, x, z, y])
            >>> for s in solver.valid_states: print s
            ({v: 1, u: 1, x: 0, z: 0, y: 1}, [])
            ({v: 1, u: 1, x: 1, z: 0, y: 0}, [])
            >>> solver.get_solutions()
            {v: 1, u: 1, x: -y + 1, z: 0}
            
            So the system is under-determined. Add another equation!
            >>> x, y = sympy.symbols('x y')
            >>> solver.add_equation(sympy.Eq(x + 2*y, 2))
            >>> solver.valid_states
            [({v: 1, u: 1, x: 0, z: 0, y: 1}, [])]
            >>> solver.get_solutions()
            {v: 1, u: 1, x: 0, z: 0, y: 1}
        '''
        # Check that we haven't gone too far over the limit by putting it
        # in the outer loop. This way we can still have the post-processing
        # with a continue statement and not check too often
        if len(self.valid_states) > max_states:
            return

        if vars_to_sub is None:
            vars_to_sub = self.variables_to_sub
        elif isinstance(vars_to_sub, int):
            vars_to_sub = self.variables_to_sub[:vars_to_sub]
        elif isinstance(vars_to_sub, sympy.Symbol):
            vars_to_sub = [vars_to_sub]
        
        self._pop_variables_from_queue(vars_to_sub)

        if len(vars_to_sub) == 0:
            return

        possible_vals = list(itertools.product([sympy.S.Zero, sympy.S.One], 
                                               repeat=len(vars_to_sub)))

        old_states = self.valid_states
        self.valid_states = []


        for state_dict, eqns in old_states:
            for possible_val in possible_vals:
                to_sub = dict(zip(vars_to_sub, possible_val))
                try:
                    new_eqns = [eqn.subs(to_sub) for eqn in eqns]
                    new_eqns = map(balance_terms, new_eqns)
                    new_eqns = filter(is_equation, new_eqns)
                    apply_contradictions(new_eqns)
                    new_state = state_dict.copy()
                    for k, v in to_sub.iteritems():
                        # Check we've not already tried a different solution
                        if new_state.get(k) not in [None, v]:
                            raise ContradictionException('{} already set to {}'.format(k, new_state.get(k)))
                        new_state[k] = v
                    self.valid_states.append((new_state, new_eqns))
                    

                except ContradictionException:
                    continue
        
        self.check_states()
        
        # Reset the sympy cache when done
        clear_cache()

    @property
    def _length_tuple(self):
        ''' Return a summary of what's going on in the solver '''
        return (len(self.variables_subbed), len(self.variables_to_sub), 
                len(self.valid_states), len(self.solutions))


    def solve_equations(self, max_iter=200, verbose=False,
                        vars_at_a_time=3, max_states=2046):
        ''' Nice wrapper around sub_var that can be used by external test to
            say 'try and get as far as you can with sensible parameters'.
            Any specialised uses of this class should use sub_var
        '''
        if verbose:        
            self.print_('Num variables: {}'.format(len(self.variables)))
            self.print_('Iterations\tNum Subbed\tNum to Sub\tNum State\tNum Sol')

        for i in xrange(max_iter):
            self.sub_var(vars_to_sub=vars_at_a_time, max_states=max_states)
            if verbose:
                self.print_('\t\t'.join(['{}'] * 5).format(i, *self._length_tuple))
            
            if len(self.variables_to_sub) == 0:
                break

    def add_equation(self, eqn, reorder_variables_to_sub=False, _init=False):
        ''' Add an equation to each of the states, substituting in known values
            and checking for contradictions
       
            >>> search = SolverSequential()
            >>> eqns = ['x + y == z', 'x*y + u + v == 3']
            >>> eqns = str_eqns_to_sympy_eqns(eqns)
            >>> for eqn in eqns:
            ...     search.add_equation(eqn, reorder_variables_to_sub=False)
            ...     print search.variables_to_sub
            [x, y, z]
            [x, y, z, u, v]
       
            >>> search = SolverSequential()
            >>> eqns = ['x + y == z', 'x*y + u + v == 3']
            >>> eqns = str_eqns_to_sympy_eqns(eqns)
            >>> for eqn in eqns:
            ...     search.add_equation(eqn, reorder_variables_to_sub=True)
            ...     print search.variables_to_sub
            [x, y, z]
            [u, v, x, y, z]
        '''
        if not is_equation(eqn):
            raise ValueError('{} tried to be added as an equation')

        # While we're at it, add the equation to self.equations provided by
        # the base class
        if not _init:
            self.equations.append(eqn)

        exclude_var = self.variables_subbed.union(set(self.variables_to_sub))
        for variable in sorted(eqn.atoms(sympy.Symbol), key=str):
            if variable in exclude_var:
                continue
            self.variables_to_sub.append(variable)
        
        # And re-order
        if reorder_variables_to_sub:
            self.reorder_variables_to_sub()

        # Initialise if we have the first equation
        if self.valid_states is None:
            self.valid_states = [(dict(), [eqn])]
            return

        old_states = self.valid_states
        self.valid_states = []
        for state_dict, eqns in old_states:
            try:
                new_eqn = eqn.subs(state_dict)
                if is_equation(new_eqn, check_true=True):
                    apply_contradictions([new_eqn])
                    self.valid_states.append((state_dict, eqns + [new_eqn]))
                else:
                    self.valid_states.append((state_dict, eqns))
            except ContradictionException:
                continue
        
        self.check_states()

    def add_solution(self, variable, value):
        ''' Override add_variable, since we don't have a solutions dict 
        
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqns = ['x + y == 1', 'u + v == 1']
            >>> eqns = str_eqns_to_sympy_eqns(eqns)
            >>> solver = SolverSequential(eqns)
            >>> solver.sub_var(None)
            >>> for s in solver.valid_states: print s
            ({v: 1, u: 0, x: 0, y: 1}, [])
            ({v: 0, u: 1, x: 0, y: 1}, [])
            ({v: 1, u: 0, x: 1, y: 0}, [])
            ({v: 0, u: 1, x: 1, y: 0}, [])
            >>> solver.add_solution(x, 1)
            >>> for s in solver.valid_states: print s
            ({v: 1, u: 0, x: 1, y: 0}, [])
            ({v: 0, u: 1, x: 1, y: 0}, [])
            
            Make sure we inject anything entirely new into the states

            >>> solver.add_solution(z, 1)
            >>> for s in solver.valid_states: print s
            ({v: 1, u: 0, x: 1, z: 1, y: 0}, [])
            ({v: 0, u: 1, x: 1, z: 1, y: 0}, [])

            Make sure if we add a conflicting solution, we throw
            >>> solver.add_solution(z, 0)
            Traceback (most recent call last):
                ...
            ContradictionException: No more valid states!
            
            Check again with a variable we know from the equations
            >>> solver = SolverSequential(eqns)
            >>> solver.sub_var(None)
            >>> solver.add_solution(x, 1)
            >>> solver.add_solution(y, 1)
            Traceback (most recent call last):
                ...
            ContradictionException: No more valid states!
        '''
        assert len(variable.atoms(sympy.Symbol)) == 1
        assert is_constant(value)
        self.add_equation(sympy.Eq(variable, value))
        self.sub_var([variable])

    @property
    def solutions(self):
        ''' Override the solutions dict '''
        return self.get_solutions()

    @solutions.setter
    def solutions(self, val):
        ''' Use add solution '''
        pass

    def add_interleaving_equations(self, deductions):
        ''' Given equations, and a dictionary of deductions, return a new list of
            equations with the deductions interleaved so the search space of
            final system is reduced as much as possible
        '''
        pass
    
    def get_solutions(self):
        ''' Work out what we know from the valid states '''
        # Check we have something to begin with
        if self.valid_states is None:
            return BinaryEquivalenceDict()
        elif len(self.valid_states) == 0:
            raise ContradictionException('No valid states')
            
        # Extract the substitution dictionaries
        valid_states = map(itemgetter(0), self.valid_states)
    
        # Intersection will be a set of (variable, value) tuples, where value
        # will be 0 or 1. We can then intersect with other non-contradictory
        # solutions so we are left with deductions that must be true
        intersection = set.intersection(*[set(state.items()) for state in valid_states])
        
        # Store our deductions somewhere
        deductions = BinaryEquivalenceDict(intersection)

        # difference_grid is a dict of (var1, var2): difference, where var1 and
        # var2 are variables and difference is the deduced difference - None
        # initially when we don't know anything about the relationship.
        # When we get contradictory relations, the tuple will be popped since
        # we can't make any deduction.
        difference_grid = {}
        
        # Check all the valid states have the same set of variables
        variables = set(valid_states[0].keys())
        
        for vars_ in itertools.combinations(variables, 2):
            difference_grid[vars_] = None
        
        for state in valid_states:
            # Process the difference relations
            for key, diff in difference_grid.copy().iteritems():
                var1, var2 = key
                # We know they can be equal                    
                if state[var1] == state[var2]:
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
                    
        # Now update the deductions with differences
        for (var1, var2), diff in difference_grid.iteritems():
            if diff == 1:
                deductions[var1] = sympy.S.One - var2
            elif diff == 0:
                deductions[var1] = var2
            else:
                # This absolutely should not happen as it should be caught
                # by checking intersection
                raise ContradictionException('This should not happen!')

        return deductions

    

if __name__ == '__main__':
    import doctest
    from sympy_helper_fns import str_eqns_to_sympy_eqns
    doctest.testmod()
