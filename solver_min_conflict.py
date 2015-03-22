# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 09:00:05 2015

@author: Richard
"""

import itertools

import numpy
import sympy

from carry_equations_generator import generate_carry_equations
from semiprime_tools import num_to_factor_num_qubit
from sympy_helper_fns import expressions_to_variables, is_constant
from verification import check_solutions

class MinConSolver(object):
    ''' Class that solves equations by the min-conflict method '''
    
    def __init__(self, equations):
        
        # Filter out any crud
        equations = filter(is_equation, equations)
        self.equations = equations
        self.atoms = numpy.array(sorted(expressions_to_variables(equations), key=str))
        self.num_atoms = len(self.atoms)
        self.state = self.get_random_state()
        
    def get_random_state(self):
        ''' Return a random starting state for the solver '''
        numpy.random.seed(seed=self.num_atoms)
        return numpy.random.random_integers(0, 1, self.num_atoms)
        
    @property
    def state_dict(self):
        ''' Return a state represented by a dictionary, ready for subbing '''
        return dict(zip(self.atoms, self.state))

    @property
    def current_num_conflicts(self):
        return sum(map(is_true, self.equations))

    def run(self, var_step=3, stuck_iter=100, num_iter=10000):
        ''' Run num_iter iterations of the algorithm substituting var_step
            variables at a time
        '''
        print self.state
        num_curr_conflict = 0
        curr_conflict = self.current_num_conflicts
        for i in xrange(num_iter):
            vars_to_sub = self.get_next_var_to_step(num_var=var_step,
                                                    num_const_iter=num_curr_conflict)
            total_num_conflicts = self.step(vars_to_step=vars_to_sub)
            print self.state, total_num_conflicts
            
            # If we've made no progress...            
            if curr_conflict == total_num_conflicts:
                num_curr_conflict += 1
                # We might want to stop...
                if num_curr_conflict == stuck_iter:
                    break
            # Or reset...
            elif curr_conflict > total_num_conflicts:
                num_curr_conflict = 0
                curr_conflict = total_num_conflicts
            # Or run for the hills
            else:
                raise ValueError('wtf?')
            
            # Or we might be done!
            if total_num_conflicts == 0:
                print 'Finished in {} iterations'.format(i)
                break
    
    def _indices_of_variables(self, variables):
        ''' Return the indices of given variables '''
        indices = numpy.searchsorted(map(str, self.atoms), map(str, variables))
        assert all(self.atoms[indices] == variables)
        return indices

    def get_next_var_to_step(self, num_var, num_const_iter=0):
        ''' Pick the next variables to play with '''
#        return numpy.random.choice(self.atoms, size=num_var, replace=False)
        state = self.state_dict    

        # Work out which equations need looking at
        eqn_eval = [eqn.subs(state) for eqn in self.equations]
        false_eqn = numpy.array([not e for e in eqn_eval])
        false_eqn = numpy.array(self.equations)[false_eqn]
        # Now find the variables in these equations
        false_var = list(expressions_to_variables(false_eqn))
        # And give them some extra weight in the sampling        
        false_ind = self._indices_of_variables(false_var)        
        weights = numpy.ones(len(self.atoms))
        weights[false_ind] += 0.1#len(false_ind)
        weights /= sum(weights)
        
        print 'p of picking false var: {}'.format(sum(weights[false_ind]))

#        print weighted_var
        var_to_step = numpy.random.choice(self.atoms, size=num_var, 
                                          replace=False, p=weights)
        return var_to_step
    
    def step(self, vars_to_step):
        ''' Find the best choice of assignments for vars_to_step '''
        # Make sure of type
        vars_to_step = list(vars_to_step)
        
        # Now substitute everything we're not interested in changing
        trimmed_state = self.state_dict
        for variable in vars_to_step:
            trimmed_state.pop(variable)
        
        subbed_eqns = [eqn.subs(trimmed_state) for eqn in self.equations]
        trimmed_eqns = filter(is_equation, subbed_eqns)

        # Now work out how many conflicts we have of stuff that isn't an
        # equation
        pre_determined = filter(lambda x: not is_equation(x), subbed_eqns)
        predet_num_true = sum(map(is_true, pre_determined))
        
        # Work out possible states
        possible_assignments = itertools.product(range(2), repeat=len(vars_to_step))

        # Now find the best one!
        best_satisfactions = 0
        best_assignments = []
        
        for assignment in possible_assignments:
            to_sub = dict(zip(vars_to_step, assignment))
            _eqns = [eqn.subs(to_sub) for eqn in trimmed_eqns]
            assert all(map(lambda x: not is_equation(x), _eqns))
            num_true = sum(map(is_true, _eqns))
            if num_true > best_satisfactions:
                best_assignments = [assignment]
                best_satisfactions = num_true
            elif num_true == best_satisfactions:
                best_assignments.append(assignment)
            else:
                continue
        
        # Now randomly pick the best assignment and do the substitution
        rand_ind = numpy.random.choice(len(best_assignments), 1)
        new_assignment = best_assignments[rand_ind]
        var_indices = self._indices_of_variables(vars_to_step)
        for index, val in zip(var_indices, new_assignment):
            self.state[index] = val
        
        # Number of conflicts is the number of equations - (number true)        
        return len(self.equations) - (predet_num_true + best_satisfactions)

class MinConSolver2(MinConSolver):
    ''' Class that solves equations by the min-conflict method. Uses distance
        from 0 as a conflict, rather than equation True or False    
    '''
    def __init__(self, equations):
        super(MinConSolver2, self).__init__(equations=equations)
        
        self.equations = [e.lhs - e.rhs for e in self.equations]
    
    @property
    def current_num_conflicts(self):
        return sum(map(abs, [eqn.subs(self.state_dict) for eqn in self.equations]))    
    
    def step(self, vars_to_step):
        ''' Find the best choice of assignments for vars_to_step '''
        # Make sure of type
        vars_to_step = list(vars_to_step)
        
        # Now substitute everything we're not interested in changing
        trimmed_state = self.state_dict
        for variable in vars_to_step:
            trimmed_state.pop(variable)
        
        subbed_eqns = [eqn.subs(trimmed_state) for eqn in self.equations]
        trimmed_eqns = filter(lambda x: not is_constant(x), subbed_eqns)

        # Now work out how many conflicts we have of stuff that isn't an
        # equation
        pre_determined = filter(is_constant, subbed_eqns)
        predet_num_conflicts = sum(map(abs, pre_determined))
        
        # Work out possible states
        possible_assignments = itertools.product(range(2), repeat=len(vars_to_step))

        # Now find the best one!
        #TODO properly
        least_conflicts = None
        best_assignments = []
        
        for assignment in possible_assignments:
            to_sub = dict(zip(vars_to_step, assignment))
            _eqns = [eqn.subs(to_sub) for eqn in trimmed_eqns]
            assert all(map(is_constant, _eqns))
            num_conflicts = sum(map(abs, _eqns))
            if (least_conflicts is None) or (num_conflicts < least_conflicts):
                best_assignments = [assignment]
                least_conflicts = num_conflicts
            elif num_conflicts == least_conflicts:
                best_assignments.append(assignment)
            else:
                continue
        
        # Now randomly pick the best assignment and do the substitution
        rand_ind = numpy.random.choice(len(best_assignments), 1)
        new_assignment = best_assignments[rand_ind]
        var_indices = self._indices_of_variables(vars_to_step)
        for index, val in zip(var_indices, new_assignment):
            self.state[index] = val
        
        # Number of conflicts is the number of equations - (number true)        
        return predet_num_conflicts + least_conflicts
    
    def get_next_var_to_step(self, num_var, num_const_iter=0):
        ''' Pick the next variables to play with.
            num_const_iter is the number of iterations we've been stuck on        
        '''
        state = self.state_dict
        
        # Work out which equations need looking at, and add weights to
        # incorrect variables
        weights = (numpy.ones(len(self.atoms)) / len(self.atoms)) * (num_const_iter + 1)

        for eqn in self.equations:
            error = abs(eqn.subs(state))
            variables = list(eqn.atoms(sympy.Symbol))
            var_indices = self._indices_of_variables(variables)
            weights[var_indices] += error / len(eqn.atoms(sympy.Symbol))
            
        # Renormalise
        weights /= sum(weights)

#        print ', '.join(['{:.3f}'.format(w) for w in weights])

#        print weighted_var
        var_to_step = numpy.random.choice(self.atoms, size=num_var, 
                                          replace=False, p=weights)
        return var_to_step

def is_equation(eqn):
    ''' Return True if it is an equation rather than a boolean value.
        If it is False, raise a ContradictionException. We never want anything
        that might be False

        >>> x, y = sympy.symbols('x y')
        >>> eq1 = sympy.Eq(x, y)
        >>> eq2 = sympy.Eq(x, x)
        >>> eq3 = sympy.Eq(x, y).subs(y, x)
        >>> expr = sympy.sympify('x + y')

        >>> is_equation(eq1)
        True
        >>> is_equation(eq2)
        False
        >>> is_equation(eq3)
        False
        >>> is_equation(expr)
        False
    '''
    if sympy.__version__ == '0.7.5':
        return isinstance(eqn, sympy.Equality)
    else:
        return eqn is True
    
def is_true(eqn):
    ''' Return True if it is an equation rather than a boolean value.
        If it is False, raise a ContradictionException. We never want anything
        that might be False

        >>> x, y = sympy.symbols('x y')
        >>> eq1 = sympy.Eq(x, y)
        >>> eq2 = sympy.Eq(x, x)
        >>> eq3 = sympy.Eq(x, y).subs(y, x)
        >>> eq4 = sympy.Eq(1, 0)

        >>> is_true(eq1)
        0
        >>> is_true(eq2)
        1
        >>> is_true(eq3)
        1
        >>> is_true(eq4)
        0
    '''
    if sympy.__version__ == '0.7.5':
        if isinstance(eqn, sympy.Equality):
            return False
        elif isinstance(eqn, sympy.boolalg.BooleanAtom):
            if eqn:
                return 1
            else:
                return 0
        else:
            raise ValueError('Invalid type {}'.format(type(eqn)))
    else:
        return eqn is True

if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    prod = 143
#    prod = 56153
    prod = 4306239659
    fact1, fact2 = num_to_factor_num_qubit(prod)
    equations = generate_carry_equations(fact1, fact2, prod)
    
    system = MinConSolver2(equations)
    system.run(stuck_iter=100, var_step=4)
    
#    check_solutions(prod, system.state_dict, verbose=True)
#    print system.state_dict
#    for e in system.equations: print e, is_equation(e), is_false(e)