# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 00:00:11 2015

@author: Richard
"""

from collections import defaultdict
import itertools
from operator import itemgetter
from random import shuffle

import sympy

from contradiction_exception import ContradictionException
from contradictions import apply_contradictions
from sympy_helper_fns import (max_value, min_value, is_equation,
                              is_constant,
                              num_add_terms, parity, expressions_to_variables,
                              standardise_equation)

class JudgementMixin(object):
    ''' Mixin object that gives access to the judgements '''

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
                    apply_contradictions([_eqn])
                
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


    @staticmethod
    def _are_equations_similar(eqns):
        ''' Given a list of equations, work out their 'alignment' or similarity.
            Equations are similar if they share lots of the same variables or
            are tangled in some way.
            Specifically: score = 
            sum((occurances of each variable - 1)**2) / num_var

            >>> system = EquationSolver()
            >>> equations = ['x == 1 - y', 
            ...              'y == 1 - z',]
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system._are_equations_similar(equations)
            0.333

            >>> equations = ['a == b', 
            ...              'x == 1 - y',]
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system._are_equations_similar(equations)
            0.0

            >>> equations = ['a == 1 + b - c',
            ...              'a + 1 == 3 - 2*b + 5*c',
            ...              'x == 1 - y',]
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system._are_equations_similar(equations)
            0.6

            >>> equations = ['a == 1 + b - c',
            ...              'a + 1 == 3 - 2*b + 5*c',
            ...              'a*b == 1 - b + d',
            ...              'x == 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7',]
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system._are_equations_similar(equations)
            0.75

            Real examples that produce deductions
            >>> equations = ['p3 + 2*q1*q2 + q3 + z23 == 2*z34 + 1', 
            ... 'p3*q5 + p4*q4 + p5*q3 + p6*q2 + 2*q1 + q2*q6 + z68 + z78 == 4*z810 + 8*z811 + 2*z89 + 1']
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system._are_equations_similar(equations)
            0.235
            
            Another real example
            >>> equations = ['p3 + 2*q1*q2 + q3 + z23 == 2*z34 + 4*z35 + 1',
            ...              'q1 + q2 == z23 + 2*z24']
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system._are_equations_similar(equations)
            0.375
        '''
        #TODO improve the scoring system, this is quite crude
        num_atoms = len(expressions_to_variables(eqns))
        atom_freq = defaultdict(int)
        for eqn in eqns:
            for atom in eqn.atoms(sympy.Symbol):
                atom_freq[atom] += 1
        
        score = sum([(v - 1)**2 for v in atom_freq.values()])
        score /= float(num_atoms)
        return round(score, 3)

    def judgement_mini_assumption_multi_eqn(self, eqns, num_var=4, 
                                  coef_transform=lambda x: pow(x, 0.01),
                                  cutoff=0.2):
        ''' Given an equation, assume the most common num_var are 0/1 and see
            if we can get any contradictions.
            coef_transform is a function used to rank variables by taking 
            coef_transform(abs(coef))
            
            Now also check to see if we can get any values that have to be
            equal or unequal - again by checking each case
            
            >>> equations = ['x == 1 - y', 
            ...              'y == 1 - z',]
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system = EquationSolver()
            >>> system.judgement_mini_assumption_multi_eqn(equations, 
            ...                       num_var=3)
            >>> system.deductions
            {x: -y + 1, z: -y + 1, y: -z + 1}

            # Check that it works even when we overload the parameters
            >>> equations = ['1 + x == 2*y + z', 
            ...              'x == 1',]
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system = EquationSolver()
            >>> system.judgement_mini_assumption_multi_eqn(equations, 
            ...                       num_var=100)
            >>> system.deductions
            {x: 1, z: 0, y: 1}


            An example where we couldn't make any deductions by considering the
            equations in isolation
            >>> equations = ['x*y*z == 0',
            ...              '3*x + 3*y == 2 + z + 3*a',]
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system = EquationSolver()

            >>> for eqn in equations:
            ...     system.judgement_mini_assumption(eqn, 
            ...                       num_var=3)
            >>> system.deductions
            {}
            >>> system.judgement_mini_assumption_multi_eqn(equations, 
            ...                       num_var=3)
            >>> system.deductions
            {x: -y + 1}


            A real example
            >>> equations = ['p3 + 2*q1*q2 + q3 + z23 == 2*z34 + 1', 
            ... 'p3*q5 + p4*q4 + p5*q3 + p6*q2 + 2*q1 + q2*q6 + z68 + z78 == 4*z810 + 8*z811 + 2*z89 + 1']
            >>> num_var = 6            
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system = EquationSolver()

            >>> for eqn in equations:
            ...     system.judgement_mini_assumption(eqn, 
            ...                       num_var=num_var)
            >>> system.deductions
            {}
            >>> system.judgement_mini_assumption_multi_eqn(equations, 
            ...                       num_var=num_var)
            >>> system.deductions
            {z811: 0}
            
            Another real example
            >>> equations = ['p3 + 2*q1*q2 + q3 + z23 == 2*z34 + 4*z35 + 1',
            ...              'q1 + q2 == z23 + 2*z24']
            >>> num_var = 4
            >>> equations = str_eqns_to_sympy_eqns(equations)
            >>> system = EquationSolver()

            >>> for eqn in equations:
            ...     system.judgement_mini_assumption(eqn, 
            ...                       num_var=num_var)
            >>> system.deductions
            {}
            >>> system.judgement_mini_assumption_multi_eqn(equations, 
            ...                       num_var=num_var)
            >>> system.deductions
            {z35: 0}
        '''
        if isinstance(eqns, sympy.Equality):
            eqns = [eqns]
        if len(eqns) == 1:
            self.judgement_mini_assumption(eqns[0], num_var=num_var, 
                                           coef_transform=coef_transform)
        elif len(eqns) == 0:
            return
        
        # Now work out if it's worth carrying out the substitution
        if self._are_equations_similar(eqns) < cutoff:
            return        
        
        # Create a dictionary of scores that we're going to use to rank variables
        var_score = defaultdict(int)
        for eqn in eqns:
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
            _eqns = [eqn.subs(to_sub) for eqn in eqns]
            try:
                for _eqn in _eqns:
                    if is_equation(_eqn):
                        _eqn = standardise_equation(_eqn)
                        apply_contradictions([_eqn])
                
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
        for (var1, var2), diff in sorted(difference_grid.iteritems(), key=itemgetter(1)):
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

if __name__ == "__main__":
    import doctest
    from sympy_solver import EquationSolver
    from sympy_helper_fns import str_eqns_to_sympy_eqns
    doctest.testmod()
