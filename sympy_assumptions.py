# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 17:33:23 2015

@author: Richard
"""

from collections import defaultdict

import sympy

from contradiction_exception import ContradictionException

def get_next_lexographical_variable(equation_solver):
    ''' Given an EquationSolver, find the next variable that we'd like to
        substitute
        
        >>> a, b, p1, p2, q4 = sympy.symbols('a b p1 p2 q4')
        
        >>> test1 = {str(v): v for v in [a, b]}
        >>> test2 = {str(v): v for v in [p1, p2, q4]}
        >>> get_next_lexographical_variable(EquationSolver(variables=test1))
        a
        >>> get_next_lexographical_variable(EquationSolver(variables=test2))
        p1
    '''
    vars_ = sorted(equation_solver.unsolved_var, key=str)
    
    if len(vars_):
        return vars_[0]
    else:
        return None

def get_next_important_variable(equation_solver):
    ''' Given an EquationSolver, find the next variable that we'd like to
        substitute, by summing abs(coef) where coef is the coefficient of any
        term the variable appears in
        
        >>> a, b, p1, p2, q4 = sympy.symbols('a b p1 p2 q4')
        
        >>> eqn1 = sympy.Eq(a + 2*b)
        >>> eqn2 = sympy.Eq(a + b - 1)
        >>> get_next_important_variable(EquationSolver(equations=[eqn1, eqn2]))
        b
        
        >>> eqn1 = sympy.Eq(a + 2*a*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> get_next_important_variable(EquationSolver(equations=[eqn1, eqn2]))
        a

        >>> eqn1 = sympy.Eq(a - 2*a*b)
        >>> eqn2 = sympy.Eq(a + b - 20)
        >>> get_next_important_variable(EquationSolver(equations=[eqn1, eqn2]))
        a
    '''
    vars_ = equation_solver.unsolved_var
    
    var_score = defaultdict(int)
    for eqn in equation_solver.equations:
        for term, coef in (eqn.lhs + eqn.rhs).as_coefficients_dict().iteritems():
            for atom in term.atoms(sympy.Symbol):
                var_score[atom] += 1#abs(coef)
    
    # If we have no variables, return None
    if not len(var_score):
        return None

    # Now extract the biggest score
    biggest = None
    biggest_score = -1
    for var, score in var_score.iteritems():
        if score > biggest_score:
            biggest_score = score
            biggest = var
    if biggest is not None:
        assert var in vars_
    return biggest



def make_assumptions(equation_solver, num_assumptions=3, 
                     assumed_variables=None, count_determined=False,
                     var_getter=get_next_important_variable):
    ''' Given a system of equations, make a number of assumptions and see if we
        reach a contradiction or not.
        Returns a list of EquationSolver systems representing each possible
        system
        Assumed variables is a list of variables that have been chosen
        If count_determined is False, then if a variable leads to a contradiction
        one way but not the other, then it doesn't count as an assumption.
        If True, a fixed number of iterations occurs
        
        >>> test_kwargs = {'var_getter': get_next_lexographical_variable,
        ...                'count_determined': True}        
        
        1.
        >>> equations = ['x+y']
        >>> equations = map(sympy.sympify, equations)
        >>> equations = map(sympy.Eq, equations)
        >>> system = EquationSolver(equations)
        >>> assumed_variables = set()
        >>> sols = make_assumptions(system, num_assumptions=1, assumed_variables=assumed_variables, **test_kwargs)
        
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
        >>> sols = make_assumptions(system, num_assumptions=1, assumed_variables=assumed_variables, **test_kwargs)
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
        >>> sols = make_assumptions(system, num_assumptions=1, assumed_variables=assumed_variables, **test_kwargs)
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
        >>> sols = make_assumptions(system, num_assumptions=2, assumed_variables=assumed_variables, **test_kwargs)
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
        >>> sols = make_assumptions(system, num_assumptions=3, assumed_variables=assumed_variables, **test_kwargs)
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
        >>> sols = make_assumptions(system, num_assumptions=3, assumed_variables=assumed_variables, **alt_kwargs)
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
        >>> sols = make_assumptions(system, num_assumptions=2, assumed_variables=assumed_variables, **test_kwargs)
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
    
    next_var = var_getter(equation_solver)
    if next_var is None:
        return [equation_solver]

    # Now we know we want to do some work, pass down the arguments
    kwargs = {'var_getter': var_getter,
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
                path_0 = make_assumptions(eqnsol0, 
                                          num_assumptions=num_assumptions - 1,
                                          assumed_variables=assumed_variables,
                                          **kwargs)
            except ContradictionException:
                return make_assumptions(eqnsol1, 
                                        num_assumptions=num_assumptions - int(count_determined),
                                        assumed_variables=assumed_variables,
                                        **kwargs)
            
            try:
                path_1 = make_assumptions(eqnsol1, 
                                          num_assumptions=num_assumptions - 1,
                                          assumed_variables=assumed_variables,
                                          **kwargs)
            except ContradictionException:
                return make_assumptions(eqnsol0, 
                                        num_assumptions=num_assumptions - int(count_determined),
                                        assumed_variables=assumed_variables,
                                        **kwargs)
            return path_0 + path_1
    
        # It's not 1, but it could be 0
        else:
            return make_assumptions(eqnsol0, 
                                    num_assumptions=num_assumptions-int(count_determined),
                                    assumed_variables=assumed_variables,
                                    **kwargs)
    
    # The variable cannot be 0
    else:
        # It must be 1
        if could_be_1:
            return make_assumptions(eqnsol1, 
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