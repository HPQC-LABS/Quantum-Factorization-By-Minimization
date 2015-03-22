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
from sympy_helper_fns import is_equation, balance_terms
from sympy.core.cache import clear_cache


class SolverSequential(object):
    ''' Class that represents a search through the solution space '''
    def __init__(self):
        # List of (state_dict, evaluated_eqns)
        self.valid_states = None
        
        # Queue of variables to substitute
        self.variables_to_sub = []
        
        # Track of variables already subbed
        self.variables_subbed = set()
    
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
    
    def check_states(self):
        ''' Check we have somewhere to keep going! 
        
            >>> search = SolverSequential()
            >>> search.check_states()
            
            >>> search.valid_states = []
            >>> search.check_states()
            Traceback (most recent call last):
                ...
            ContradictionException: No more valid states!
            
            >>> x = sympy.sympify('x')            
            >>> search.valid_states = [{x: 1}, []]
            >>> search.check_states()            
        '''
        if self.valid_states is None:
            return
        elif len(self.valid_states) == 0:
            raise ContradictionException('No more valid states!')
        else:
            return
    
    def sub_var(self, num_var=5):
        ''' Substitute more variables into each system '''
        if num_var is None:
            vars_to_sub, self.variables_to_sub = self.variables_to_sub, []
        else:
            vars_to_sub, self.variables_to_sub = (self.variables_to_sub[:num_var],
                                                  self.variables_to_sub[num_var:])

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
                    new_state.update(to_sub)
                    self.valid_states.append((new_state, new_eqns))
                except ContradictionException:
                    continue
        
        assert not self.variables_subbed.intersection(set(vars_to_sub))
        self.variables_subbed = self.variables_subbed.union(set(vars_to_sub))
        self.check_states()
        
        # Reset the sympy cache when done
        clear_cache()
    
    def add_equation(self, eqn):
        ''' Add an equation to each of the states, substituting in known values
            and checking for contradictions
            
            >>> search = SolverSequential()
            >>> eqns = ['x + y == z', 'x*y + u + v == 3']
            >>> eqns = str_eqns_to_sympy_eqns(eqns)
            >>> for eqn in eqns:
            ...     search.add_equation(eqn)
            ...     print search.variables_to_sub
            [x, y, z]
            [u, v, x, y, z]
        '''
        exclude_var = self.variables_subbed.union(set(self.variables_to_sub))
        for variable in sorted(eqn.atoms(sympy.Symbol), key=str):
            if variable in exclude_var:
                continue
            self.variables_to_sub.append(variable)
        # And re-order
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
                if is_equation(new_eqn):
                    apply_contradictions([new_eqn])
                self.valid_states.append((state_dict, eqns + [new_eqn]))
            except ContradictionException:
                continue
        
        self.check_states()
    
    def get_deductions(self):
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
    
#    from semiprime_tools import num_to_factor_num_qubit
#    from carry_equations_generator import generate_carry_equations
##    prod = 143
#    prod = 56153
##    prod = 4306239659
##    prod = EXPERIMENTS[15].product
# 
#    fact1, fact2 = num_to_factor_num_qubit(prod)
#    equations = generate_carry_equations(fact1, fact2, prod)
#    eqns = filter(is_equation, equations)
#    eqns_copy = eqns[:]
#
#    search = SolverSequential()
#    
#    print 'Num Qubits Start: {}'.format(len(expressions_to_variables(eqns)))    
#    from time import time
#    s = time()    
#    while(len(eqns_copy)):
#        search.add_equation(eqns_copy.pop(0))
##        if len(eqns_copy):
##            search.add_equation(eqns_copy.pop(0))
#
#        search.sub_var(num_var=None)
##        search.sub_var(num_var=4)
#
##        if len(eqns_copy):
##            search.add_equation(eqns_copy.pop())
##            search.sub_var(num_var=None)
##        print search.get_deductions()
#        print '{}\t{}'.format(len(search.valid_states), 
#                              len(search.get_deductions()),)
#    e = time()
#    print 'Solved in {:.3f}s'.format(e-s)