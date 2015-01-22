# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 17:33:23 2015

@author: Richard
"""

from collections import defaultdict
from operator import itemgetter
from time import time
import itertools
import math
import sympy

from contradiction_exception import ContradictionException
from sympy_helper_fns import is_constant

ATOMIC_VARIABLES = False

def lexographical_rank_variable(equation_solver):
    ''' Given an EquationSolver, find the next variable that we'd like to
        substitute
        
        >>> a, b, p1, p2, q4 = sympy.symbols('a b p1 p2 q4')
        
        >>> test1 = {str(v): v for v in [a, b]}
        >>> test2 = {str(v): v for v in [p1, p2, q4]}
        >>> lexographical_rank_variable(EquationSolver(variables=test1))
        {b: 0, a: 1}
        >>> lexographical_rank_variable(EquationSolver(variables=test2))
        {p1: 2, p2: 1, q4: 0}
    '''
    vars_ = sorted(equation_solver.unsolved_var, key=str, reverse=True)
    return {var: i for i, var in enumerate(vars_)}


def weighted_frequency_rank_variables(equation_solver, 
                                      atomic_variables=ATOMIC_VARIABLES):
    ''' Given an EquationSolver, find the next variable that we'd like to
        substitute, by summing abs(coef) where coef is the coefficient of any
        term the variable appears in
        
        >>> a, b = sympy.symbols('a b')
        
        >>> eqn1 = sympy.Eq(a + 2*b)
        >>> eqn2 = sympy.Eq(a + b - 1)
        >>> weighted_frequency_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=True)
        defaultdict(<type 'int'>, {b: 3, a: 2})
        
        >>> eqn1 = sympy.Eq(a + 2*a*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> weighted_frequency_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=True)
        defaultdict(<type 'int'>, {b: 3, a: 4})

        >>> eqn1 = sympy.Eq(a - 2*a*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> weighted_frequency_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=True)
        defaultdict(<type 'int'>, {b: 3, a: 4})
        
        >>> eqn1 = sympy.Eq(a - 2*a*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> weighted_frequency_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=False)
        defaultdict(<type 'int'>, {a*b: 2, b: 1, a: 2})
    '''
    var_score = defaultdict(int)
    for eqn in equation_solver.equations:
        for term, coef in (eqn.lhs + eqn.rhs).as_coefficients_dict().iteritems():

            if atomic_variables:
                ## Only use the original variables
                for atom in term.atoms(sympy.Symbol):
                    var_score[atom] += abs(coef)
            else:
                ## Use products of variables for substitutions
                if not is_constant(term):
                    var_score[term] += abs(coef)

    return var_score

def frequency_rank_variables(equation_solver, 
                             atomic_variables=ATOMIC_VARIABLES):
    ''' Given an EquationSolver, find the next variable that we'd like to
        substitute, by summing the number of times a variable appears
        
        >>> a, b = sympy.symbols('a b')
        
        >>> eqn1 = sympy.Eq(a + 2*b)
        >>> eqn2 = sympy.Eq(a + b - 1)
        >>> frequency_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=True)
        defaultdict(<type 'int'>, {b: 2, a: 2})
        
        >>> eqn1 = sympy.Eq(a + 2*a*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> frequency_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=True)
        defaultdict(<type 'int'>, {b: 2, a: 3})

        >>> eqn1 = sympy.Eq(a - 2*a*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> frequency_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=True)
        defaultdict(<type 'int'>, {b: 2, a: 3})
        
        >>> eqn1 = sympy.Eq(a - 2*a*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> frequency_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=False)
        defaultdict(<type 'int'>, {a*b: 1, b: 1, a: 2})
    '''
    var_score = defaultdict(int)
    for eqn in equation_solver.equations:
        for term, coef in (eqn.lhs + eqn.rhs).as_coefficients_dict().iteritems():
            if atomic_variables:            
                for atom in term.atoms(sympy.Symbol):
                    var_score[atom] += 1
            else:
                if not is_constant(term):
                    var_score[term] += 1
    
    return var_score
    
def max_coef_rank_variables(equation_solver, 
                            atomic_variables=ATOMIC_VARIABLES):
    ''' Given an EquationSolver, find the next variable that we'd like to
        substitute, by summing the number of times a variable appears
        
        >>> a, b = sympy.symbols('a b')
        
        >>> eqn1 = sympy.Eq(a + 2*b)
        >>> eqn2 = sympy.Eq(a + b - 1)
        >>> max_coef_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=True)
        defaultdict(<type 'int'>, {b: 2, a: 1})
        
        >>> eqn1 = sympy.Eq(a - 2*a*b + 6*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> max_coef_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=True)
        defaultdict(<type 'int'>, {b: 6, a: 2})

        >>> eqn1 = sympy.Eq(a - 2*a*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> max_coef_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=True)
        defaultdict(<type 'int'>, {b: 2, a: 2})
        
        >>> eqn1 = sympy.Eq(a - 2*a*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> max_coef_rank_variables(EquationSolver(equations=[eqn1, eqn2]), atomic_variables=False)
        defaultdict(<type 'int'>, {a*b: 2, b: 1, a: 2})
    '''
    var_score = defaultdict(int)
    for eqn in equation_solver.equations:
        for term, coef in (eqn.lhs + eqn.rhs).as_coefficients_dict().iteritems():
            if atomic_variables:            
                for atom in term.atoms(sympy.Symbol):
                    var_score[atom] = max(abs(coef), var_score[atom])
            else:
                if not is_constant(term):
                    var_score[term] += max(abs(coef), var_score[term])
    
    return var_score
    
def _rank_dict(term_to_score):
    ''' Given a dict of scores, transform them into a rank 
    
        >>> x1, x2, x3, x4 = sympy.symbols('x1 x2 x3 x4')
        >>> _rank_dict({x1: 2, x2: -1, x3: 5, x4: 2})
        {x3: 1, x4: 2, x1: 3, x2: 4}
    '''
    # Now extract the biggest score
    ranked = sorted(term_to_score.iteritems(), key=itemgetter(1), reverse=True)
    sorted_var = itertools.imap(itemgetter(0), ranked)
    return dict(((term, rank + 1) for rank, term in enumerate(sorted_var)))

def _comb_rank_func((vars_, scores)):
    ''' Transform a (variables, scores) tuple into a single sortable value '''
    return sum([math.pow(score, 0.9) for score in scores])

def get_variables_to_substitute(num_variables, equation_solver,
                                rank_func=weighted_frequency_rank_variables,
                                limit_permutations=1):
    ''' Call the ranking function to rank all of the possible substitutions
        and then return the top results.
        Returns a nested array, where each array represents a list of variables
        to substitute simultaneously.
        
        >>> equations = ['4*x1 + 3*x2 + 2*x3 + 2*x4 + x5 + x6 + x7 - 6']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)

        >>> subs = get_variables_to_substitute(1, system, limit_permutations=None)
        >>> for s in subs: print s
        (x1,)
        (x2,)
        (x3,)
        (x4,)
        (x7,)
        (x6,)
        (x5,)
        
        
        >>> subs = get_variables_to_substitute(2, system, limit_permutations=None)
        >>> for s in subs: print s
        (x2, x1)
        (x3, x1)
        (x4, x1)
        (x3, x2)
        (x7, x1)
        (x2, x4)
        (x6, x1)
        (x7, x2)
        (x3, x4)
        (x1, x5)
        (x2, x6)
        (x3, x7)
        (x2, x5)
        (x3, x6)
        (x7, x4)
        (x3, x5)
        (x6, x4)
        (x4, x5)
        (x7, x6)
        (x7, x5)
        (x6, x5)
    '''
    var_score = rank_func(equation_solver)

    # If we have no variables, return None
    if not len(var_score):
        return None

    # Base it off of variables ranks
    var_score = _rank_dict(var_score)
    
    # Now take the dictionary to a numpy array for fast ranking etc.
    var_score = var_score.items()
    
    # Now make tuples of ((x1, x2, ...), (score1, score2, ...))
    combs_to_sub = list(itertools.combinations(var_score, num_variables))
    combs_to_sub = [zip(*comb) for comb in combs_to_sub]
    
    # Now sort based on the combined score of (score1, score2, score3...)
    combs_to_sub = sorted(combs_to_sub, key=_comb_rank_func)
    combs_to_sub = map(itemgetter(0), combs_to_sub)


    if limit_permutations is None:
        return combs_to_sub
    else:
        return combs_to_sub[:limit_permutations]

def make_simultaneous_assumptions(equation_solver, num_assumptions=3,
                     rank_func=weighted_frequency_rank_variables,
                     limit_permutations=1,
                     verbose=True,
                     return_variables=False):
    ''' Given a system of equations, make a number of assumptions and see if we
        reach a contradiction or not.
        Returns a list of EquationSolver systems representing each possible
        system
        limit_permutations decides whether we just substitute the best
        num_assumptions variables, or every combination of num_assumptions
        variables.

        >>> test_kwargs = {'rank_func': lexographical_rank_variable,
        ...                'verbose': False,
        ...                'limit_permutations': 1}        
        
        1.
        >>> equations = ['x+y']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> assumed_variables = set()
        >>> sols = make_simultaneous_assumptions(system, 
        ...                                      num_assumptions=1, 
        ...                                      **test_kwargs)
        
        >>> for sol in sols: print sol.solutions
        {x: 0, y: 0}

        2.
        >>> equations = ['x + y - 1']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> assumed_variables = set()
        >>> sols = make_simultaneous_assumptions(system, 
        ...                                      num_assumptions=1, 
        ...                                      **test_kwargs)

        >>> for sol in sols: print sol.solutions
        {x: 0, y: 1}
        {x: 1, y: 0}
        
        3.
        >>> equations = ['x + y - 1', 'x + y']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> sols = make_simultaneous_assumptions(system, 
        ...                                      num_assumptions=1, 
        ...                                      **test_kwargs)
        Traceback (most recent call last):
            ...
        ContradictionException: No substitution is consistent


        4. Test making 2 simultaneous substitutions
        Make the equation really long so judgement_n_term doesn't make anything
        too ugly
        >>> equations = ['x1 + x2 + x3 + x4 + x5 + x6 + x7 - 2']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> sols = make_simultaneous_assumptions(system, 
        ...                                      num_assumptions=2,
        ...                                      **test_kwargs)
        >>> for sol in sols: print sol.solutions
        {x1: 0, x2: 0}
        {x1: 1, x2: 0}
        {x1: 0, x2: 1}
        {x3: 0, x7: 0, x2: 1, x6: 0, x4: 0, x1: 1, x5: 0}


        5. Test making all the different combinations of 1 variable
        >>> alt_kwargs = test_kwargs.copy()
        >>> alt_kwargs['limit_permutations'] = None

        Make the equation really long so judgement_n_term doesn't make anything
        too ugly
        >>> equations = ['x1 + x2 + x3 + x4 + x5 + x6 + x7 - 2']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> sols = make_simultaneous_assumptions(system, 
        ...                                      num_assumptions=1,
        ...                                      **alt_kwargs)
        >>> for sol in sols: print sol.solutions
        {x1: 0}
        {x1: 1}
        {x2: 0}
        {x2: 1}
        {x3: 0}
        {x3: 1}
        {x4: 0}
        {x4: 1}
        {x5: 0}
        {x5: 1}
        {x6: 0}
        {x6: 1}
        {x7: 0}
        {x7: 1}

        6. Test making all the different combinations of 2 variables
        NOTE: Uses 5's EquationSolver

        >>> sols = make_simultaneous_assumptions(system, 
        ...                                      num_assumptions=2,
        ...                                      **alt_kwargs)
        >>> for sol in sols[:15]: print sol.solutions
        {x1: 0, x2: 0}
        {x1: 1, x2: 0}
        {x1: 0, x2: 1}
        {x3: 0, x7: 0, x2: 1, x6: 0, x4: 0, x1: 1, x5: 0}
        {x3: 0, x1: 0}
        {x3: 0, x1: 1}
        {x3: 1, x1: 0}
        {x3: 1, x7: 0, x2: 0, x6: 0, x4: 0, x1: 1, x5: 0}
        {x4: 0, x1: 0}
        {x4: 0, x1: 1}
        {x4: 1, x1: 0}
        {x3: 0, x7: 0, x2: 0, x6: 0, x4: 1, x1: 1, x5: 0}
        {x3: 0, x2: 0}
        {x3: 0, x2: 1}
        {x3: 1, x2: 0}
    '''

    vars_to_sub = get_variables_to_substitute(num_variables=num_assumptions,
                                              equation_solver=equation_solver,
                                              rank_func=rank_func,
                                              limit_permutations=limit_permutations)

    equation_solvers = []
    times = []
    
    # Keep a record of what we're doing
    subs_record = []
    
    for var_combination in vars_to_sub:

        if verbose:
            print '\nSubstituting: ' + ', '.join(map(str, var_combination))
        
        for sub_values in itertools.product(xrange(2), repeat=num_assumptions):
            start = time()        
            try:
                solver = equation_solver.copy()
                for var, val in itertools.izip(var_combination, sub_values):
                    solver.update_value(var, val)
                solver.solve_equations()
                equation_solvers.append(solver)
                end = time()
                time_taken = end - start
                times.append(time_taken)
                if verbose:
                    print 'Substitution completed in\t{:.2f}s'.format(time_taken)
                
                # Now add to the log
                subs_record.append(dict(zip(var_combination, sub_values)))
            except ContradictionException:
                end = time()
                time_taken = end - start
                times.append(time_taken)
                if verbose:
                    print 'Contradiction reached in\t{:.2f}s'.format(time_taken)
    
    if not len(equation_solvers):
        raise ContradictionException('No substitution is consistent')
    
    if verbose:
        print '{} substitutions made in {:.2f}s'.format(len(times), sum(times))
    
    if return_variables:
        return equation_solvers, subs_record
    else:
        return equation_solvers


def make_sequential_assumptions(equation_solver, num_assumptions=3, 
                     assumed_variables=None, count_determined=False,
                     rank_func=weighted_frequency_rank_variables):
    ''' Given a system of equations, make a number of assumptions and see if we
        reach a contradiction or not.
        Returns a list of EquationSolver systems representing each possible
        system
        Assumed variables is a list of variables that have been chosen
        If count_determined is False, then if a variable leads to a contradiction
        one way but not the other, then it doesn't count as an assumption.
        If True, a fixed number of iterations occurs
        
        >>> test_kwargs = {'rank_func': lexographical_rank_variable,
        ...                'count_determined': True}        
        
        1.
        >>> equations = ['x+y']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> assumed_variables = set()
        >>> sols = make_sequential_assumptions(system, num_assumptions=1, assumed_variables=assumed_variables, **test_kwargs)
        
        >>> assumed_variables
        set([x])
        >>> for sol in sols: print sol.solutions
        {x: 0, y: 0}
        
        2.
        >>> equations = ['x + y - 1']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> assumed_variables = set()
        >>> sols = make_sequential_assumptions(system, num_assumptions=1, assumed_variables=assumed_variables, **test_kwargs)
        >>> assumed_variables
        set([x])
        >>> for sol in sols: print sol.solutions
        {x: 0, y: 1}
        {x: 1, y: 0}
        
        3.
        >>> equations = ['x + y - 1', 'x + y']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> assumed_variables = set()
        >>> sols = make_sequential_assumptions(system, num_assumptions=1, assumed_variables=assumed_variables, **test_kwargs)
        Traceback (most recent call last):
            ...
        ContradictionException: Contradiction either way when substituting x

        4.
        Note you might need to add more z variables to stop the n_term
        substitution judgement coming in
        >>> equations = ['x + y + z1 + z2 + z3 + z4 + z5 + z6 - 1', 'a + b']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> assumed_variables = set()
        >>> sols = make_sequential_assumptions(system, num_assumptions=2, assumed_variables=assumed_variables, **test_kwargs)
        >>> assumed_variables
        set([x, a])
        
        Note one path has solved a lot more variables than the other
        >>> for sol in sols: print sol.solutions
        {x: 0, b: 0, a: 0}
        {z6: 0, z5: 0, z1: 0, z4: 0, x: 1, z2: 0, y: 0, a: 0, z3: 0, b: 0}
        
        And we can see this in the remaining equations
        >>> for sol in sols: print sol.equations
        [y + z1 + z2 + z3 + z4 + z5 + z6 == 1]
        []

        Now we try the above example with 3 guesses
        >>> assumed_variables = set()
        >>> sols = make_sequential_assumptions(system, num_assumptions=3, assumed_variables=assumed_variables, **test_kwargs)
        >>> assumed_variables
        set([x, y, a])
        
        Note one path has solved a lot more variables than the other
        >>> for sol in sols: print sol.solutions
        {a: 0, x: 0, b: 0, y: 0}
        {z6: 0, z5: 0, z1: 0, z4: 0, x: 0, z2: 0, y: 1, a: 0, z3: 0, b: 0}
        {z6: 0, z5: 0, z1: 0, z4: 0, x: 1, z2: 0, y: 0, a: 0, z3: 0, b: 0}
        
        And we can see this in the remaining equations
        >>> for sol in sols: print sol.equations
        [z1 + z2 + z3 + z4 + z5 + z6 == 1]
        []
        []
        
        
        5.
        Using the same example, we show the fancy count determined feature:
        >>> alt_kwargs = test_kwargs.copy()
        >>> alt_kwargs['count_determined'] = False
        >>> assumed_variables = set()
        >>> sols = make_sequential_assumptions(system, num_assumptions=3, assumed_variables=assumed_variables, **alt_kwargs)
        >>> assumed_variables
        set([x, y, a, z1])
        
        Note one path has solved a lot more variables than the other
        >>> for sol in sols: print sol.solutions
        {a: 0, x: 0, b: 0, y: 0, z1: 0}
        {z6: 0, z5: 0, z1: 1, y: 0, x: 0, z2: 0, z4: 0, a: 0, z3: 0, b: 0}
        {z6: 0, z5: 0, z1: 0, z4: 0, x: 0, z2: 0, y: 1, a: 0, z3: 0, b: 0}
        {z6: 0, z5: 0, z1: 0, z4: 0, x: 1, z2: 0, y: 0, a: 0, z3: 0, b: 0}
        
        And we can see this in the remaining equations
        >>> for sol in sols: print sol.equations
        [z2 + z3 + z4 + z5 + z6 == 1]
        []
        []
        []
        
        
        6.
        Here is an example where we go down one road (a=1), but it's a dead end
        >>> equations = ['a + b - 1', 'a + x + y - 1', 'a + x + 2*y - 2']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> assumed_variables = set()
        >>> sols = make_sequential_assumptions(system, num_assumptions=2, assumed_variables=assumed_variables, **test_kwargs)
        >>> assumed_variables
        set([a])
        >>> for sol in sols: print sol.solutions
        {y: 1, x: 0, b: 1, a: 0}
        
        >>> for sol in sols: print sol.equations
        []
    '''
    if num_assumptions == 0:
        return [equation_solver]

    if assumed_variables is None:
        assumed_variables = set()
    
    next_var = get_variables_to_substitute(1, equation_solver, 
                                           rank_func=rank_func,
                                           limit_permutations=1)
    if next_var is None:
        return [equation_solver]
    else:
        # Now next_var is a nested list of groups of variables to substitute.
        # For sequential substitution we only want 1 variable, so extract it
        next_var = next_var[0][0]

    # Now we know we want to do some work, pass down the arguments
    kwargs = {'rank_func': rank_func,
              'count_determined': count_determined}

    assumed_variables.add(next_var)
    
    # Make systems for each assumption
    eqnsol0 = equation_solver.copy()
    eqnsol0.update_value(next_var, 0)
    eqnsol1 = equation_solver.copy()
    eqnsol1.update_value(next_var, 1)

    # Now let the systems rip!
    try:
        eqnsol0.solve_equations()
        could_be_0 = True
    except ContradictionException as contradiction_0:
        could_be_0 = False
    try:
        eqnsol1.solve_equations()
        could_be_1 = True
    except ContradictionException as contradiction_1:
        could_be_1 = False
    
    # It could be 0
    if could_be_0:
        # Suppose it could be either
        if could_be_1:
            # Then get the results from both cases and concatenate them
            # If we try going down a path and find a contradiction either way,
            # then it can't be that one!
            try:
                path_0 = make_sequential_assumptions(eqnsol0, 
                                          num_assumptions=num_assumptions - 1,
                                          assumed_variables=assumed_variables,
                                          **kwargs)
            except ContradictionException:
                return make_sequential_assumptions(eqnsol1, 
                                        num_assumptions=num_assumptions - int(count_determined),
                                        assumed_variables=assumed_variables,
                                        **kwargs)
            
            try:
                path_1 = make_sequential_assumptions(eqnsol1, 
                                          num_assumptions=num_assumptions - 1,
                                          assumed_variables=assumed_variables,
                                          **kwargs)
            except ContradictionException:
                return make_sequential_assumptions(eqnsol0, 
                                        num_assumptions=num_assumptions - int(count_determined),
                                        assumed_variables=assumed_variables,
                                        **kwargs)
            return path_0 + path_1
    
        # It's not 1, but it could be 0
        else:
            return make_sequential_assumptions(eqnsol0, 
                                    num_assumptions=num_assumptions-int(count_determined),
                                    assumed_variables=assumed_variables,
                                    **kwargs)
    
    # The variable cannot be 0
    else:
        # It must be 1
        if could_be_1:
            return make_sequential_assumptions(eqnsol1, 
                                    num_assumptions=num_assumptions-int(count_determined),
                                    assumed_variables=assumed_variables,
                                    **kwargs)
            
    
        # If it can't be either, then we must have done something wrong
        else:            
            error_str = 'Contradiction either way when substituting {}'.format(next_var)
            raise ContradictionException(error_str)
        

if __name__ == "__main__":
    import doctest
    import sympy
    from sympy_solver import EquationSolver
    doctest.testmod()