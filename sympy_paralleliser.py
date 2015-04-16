# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 21:09:13 2015

@author: Richard
"""

import itertools
import multiprocessing
import sys

from sympy.core.cache import clear_cache

from sympy_subs import subs_many
from contradiction_exception import ContradictionException

DEFAULT_NUM_PROCESSES = {'linux': 6,
                         'linux2': 6,
                         'win32': 3,
                         }
# Default to 3
DEFAULT_NUM_PROCESSES = DEFAULT_NUM_PROCESSES.get(sys.platform, 3)

def get_pool(processes=DEFAULT_NUM_PROCESSES):
    ''' Return a pool object '''
    return multiprocessing.Pool(processes=processes)

def _subs_wrap((eqns, sub_tuple)):
    ''' Wrapper function to allow parallelisation '''
    subbed = subs_many(eqns, sub_tuple)
    clear_cache()
    return subbed

def _solve_equations_wrap(solver):
    # Make sure we don't parallelise again
    solver.close_pool()
    solver.parallelise = False
    try:
        solver.solve_equations()
    except ContradictionException:
        return None
    return solver

def paralellised_subs(equations, to_sub, processes=DEFAULT_NUM_PROCESSES, pool=None):
    ''' Take a list of equations and a substitution dict and parallelise the
        substitution.
        If a pool is given, ignore processes

        >>> x, y, z = sympy.symbols('x y z')
        >>> eqns = [x + y - 1,
        ...         x*z - 1,
        ...         x - 1]
        >>> eqns = map(balance_terms, map(sympy.Eq, eqns))
        >>> eqns = [(eqn,) for eqn in eqns]
        >>> subs = {x: 1, z: 2}
        >>> subbed = paralellised_subs(eqns, subs)
        >>> for eqn in subbed: print eqn
        [y + 1 == 1]
        [False]
        [True]
    '''
    if isinstance(to_sub, dict):
        to_sub = to_sub.items()
    
    if pool is None:
        pool = get_pool(processes=processes)
        close_on_exit = True
    else:
        close_on_exit = False

    sub_copies = itertools.repeat(to_sub, len(equations))
    result = pool.map(_subs_wrap, zip(equations, sub_copies))

    # If we made the pool, shut it down again
    if close_on_exit:
        pool.close()
        pool.join()
        pool = None
    return result

def paralellised_solve_equations(solvers, processes=DEFAULT_NUM_PROCESSES, pool=None):
    ''' Take a list of solvers and call each of their solve_equations methods
        If a pool is given, ignore processes

        >>> prod = 143
        >>> dim1, dim2 = num_to_factor_num_qubit(prod)
        >>> eqns = generate_carry_equations(dim1, dim2, prod)
        >>> p1, q1 = sympy.symbols('p1 q1')
        >>> solvers = []
        >>> for pv, qv in itertools.product(range(2), repeat=2):
        ...     solver = SOLVER(eqns)
        ...     solver.add_solution(p1, pv)
        ...     solver.add_solution(q1, qv)
        ...     solvers.append(solver)
        >>> solvers = paralellised_solve_equations(solvers)
        >>> for solver in solvers:
        ...     if solver is None:
        ...         print None
        ...     else:
        ...         print len(solver.unsolved_var)
        None
        0
        0
        None
    '''
    if pool is None:
        pool = get_pool(processes=processes)
        close_on_exit = True
    else:
        close_on_exit = False

    result = pool.map(_solve_equations_wrap, solvers)

    # If we made the pool, shut it down again
    if close_on_exit:
        pool.close()
        pool.join()
        pool = None
    return result
    
if __name__ == "__main__":
    import doctest
    import sympy
    import itertools
    from sympy_helper_fns import balance_terms
    from carry_equations_generator import generate_carry_equations
    from semiprime_tools import num_to_factor_num_qubit
    from solver_hybrid import SolverHybrid
    SOLVER = SolverHybrid
    doctest.testmod()