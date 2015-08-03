"""
Created on Tue Jul 28 17:26:33 2015

@author: richard
"""

from collections import MutableSet, defaultdict, OrderedDict
import re

import sympy

from sympy_helper_fns import expressions_to_variables

def is_carry_variable(variable):
    ''' Return whether a variable is a carry variable or not

        >>> variables = sympy.symbols('x, y, a, p1, p, q57, z, z468')
        >>> for v in variables:
        ...     print v, '\t', is_carry_variable(v)
        x        False
        y        False
        a        False
        p1        False
        p        False
        q57        False
        z        False
        z468        True
    '''
    return re.match('z\d+', str(variable)) is not None

def filter_all(variables):
    ''' Filter out Nones

        >>> for e in eqns:
        ...     print e, '\t', filter_all(e.atoms(sympy.Symbol))
        p1 + q1 == 2*z12 + 1        [z12, p1, q1]
        p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1        [z12, q1, z24, q2, p2, z23, p1]
        p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1        [p2, q1, z35, p1, z34, z23, q2]
    '''
    return filter(None, variables)

def filter_carry_variables(variables):
    ''' Filter out noncarry variables

        >>> for e in eqns:
        ...     print e, '\t', filter_carry_variables(e.atoms(sympy.Symbol))
        p1 + q1 == 2*z12 + 1        [z12]
        p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1        [z12, z24, z23]
        p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1        [z35, z34, z23]
    '''
    return filter(is_carry_variable, variables)

def filter_noncarry_variables(variables):
    ''' Filter out carry variables

        >>> for e in eqns:
        ...     print e, '\t', filter_noncarry_variables(e.atoms(sympy.Symbol))
        p1 + q1 == 2*z12 + 1        [p1, q1]
        p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1        [q1, q2, p2, p1]
        p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1        [p2, q1, p1, q2]
    '''
    return filter(lambda v: not is_carry_variable(v), variables)


class EquationGraph(MutableSet):
    ''' An object for holding equations, that can be used for ordering them
        or providing neighbours etc
    '''

    def __init__(self, equations=None):
        if equations is None:
            equations = set()

        self.equations = set(equations)

    def add_equation(self, eqn):
        self.equations.append(eqn)

    def __iter__(self):
        return iter(self.equations)

    def __getitem__(self, i):
        return self.equations.__getitem__(i)

    def __delitem__(self, i):
        self.equations.__delitem__(i)

    def __contains__(self, i):
        return self.equations.__contains__(i)

    def __len__(self):
        return self.equations.__len__()

    def add(self, obj):
        self.equations.add(obj)

    def discard(self, obj):
        self.equations.discard(obj)

    def variables(self):
        ''' Return a set of the variables that appear in the equations

        >>> print eqn_graph.variables()
        set([z12, q1, q2, z35, z24, p2, z34, z23, p1])
        '''
        return expressions_to_variables(self.equations)


    def variables_to_equations(self):
        ''' Return a {variable: set(eqns that contain variable)} mapping

            >>> for v, eqns in eqn_graph.variables_to_equations().iteritems():
            ...     print v
            ...     for e in eqns: print e
            z35
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            z24
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            z34
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            z12
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            p1 + q1 == 2*z12 + 1
            p2
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            z23
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            q2
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            q1
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            p1 + q1 == 2*z12 + 1
            p1
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            p1 + q1 == 2*z12 + 1
        '''
        variables = self.variables()

        variables_to_eqns = defaultdict(set)
        for v in variables:
            for e in self.equations:
                if v in e.atoms():
                    variables_to_eqns[v].add(e)

        variables_to_eqns = OrderedDict(sorted(variables_to_eqns.iteritems(),
                                               key=lambda x: len(x[1])))

        return variables_to_eqns

    def variables_to_neighbours(self, filter_func=filter_all):
        ''' Return a {variable: neighbours} map, sorted by how many neighbours
            there are

            >>> for v, n in eqn_graph.variables_to_neighbours().iteritems():
            ...     print v, n
            z35 set([q1, p2, z34, z23, z35, p1, q2])
            z24 set([z12, p2, q1, z23, p1, z24, q2])
            z34 set([q1, p2, z34, z23, z35, p1, q2])
            z12 set([z12, p2, q1, z23, p1, z24, q2])
            p2 set([z12, q1, q2, z35, z24, p2, z34, z23, p1])
            z23 set([z12, q1, q2, z35, z24, p2, z34, z23, p1])
            q2 set([z12, q1, q2, z35, z24, p2, z34, z23, p1])
            q1 set([z12, q1, q2, z35, z24, p2, z34, z23, p1])
            p1 set([z12, q1, q2, z35, z24, p2, z34, z23, p1])

            >>> print eqn_graph.variables_to_neighbours().keys()
            [z35, z24, z34, z12, p2, z23, q2, q1, p1]

            >>> filter = filter_carry_variables
            >>> for v, n in eqn_graph.variables_to_neighbours(filter_func=filter).iteritems():
            ...     print v, n, len(filter(n))
            z35 set([q1, p2, z34, z23, z35, p1, q2]) 3
            z24 set([z12, p2, q1, z23, p1, z24, q2]) 3
            z34 set([q1, p2, z34, z23, z35, p1, q2]) 3
            z12 set([z12, p2, q1, z23, p1, z24, q2]) 3
            p2 set([z12, q1, q2, z35, z24, p2, z34, z23, p1]) 5
            z23 set([z12, q1, q2, z35, z24, p2, z34, z23, p1]) 5
            q2 set([z12, q1, q2, z35, z24, p2, z34, z23, p1]) 5
            q1 set([z12, q1, q2, z35, z24, p2, z34, z23, p1]) 5
            p1 set([z12, q1, q2, z35, z24, p2, z34, z23, p1]) 5
        '''
        variables_to_neighbours = []
        for v, s in self.variables_to_equations().iteritems():
            variables_to_neighbours.append((v, expressions_to_variables(s)))

        variables_to_neighbours.sort(key=lambda x: len(filter_func(x[1])))
        return OrderedDict(variables_to_neighbours)

    def smallest_equation_intersecting_root(self, root, filter_func=filter_all):
        ''' Given a set of variables, find the equation with the smallest
            number of variables in addition to the root

            >>> for e in eqn_graph.equations:
            ...     print e
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            p1 + q1 == 2*z12 + 1

            >>> root = set(sympy.symbols('p1 p2'))
            >>> print eqn_graph.smallest_equation_intersecting_root(root)
            p1 + q1 == 2*z12 + 1

            >>> root = set([sympy.symbols('z12')])
            >>> print eqn_graph.smallest_equation_intersecting_root(root)
            p1 + q1 == 2*z12 + 1

            >>> print eqn_graph.smallest_equation_intersecting_root(root, filter_func=filter_carry_variables)
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
        '''
        # Find only equations that intersect our root, or everything if it's empty
        if len(root):
            intersecting_eqns = filter(lambda e: len(e.atoms(sympy.Symbol).intersection(root)), self.equations)
            intersecting_eqns = set(intersecting_eqns)
        else:
            intersecting_eqns = set(self.equations)
        # Now find only equations that have something interesting left to substitute
        candidate_eqns = filter(lambda e: len(filter_func(e.atoms(sympy.Symbol).difference(root))), intersecting_eqns)
        if candidate_eqns:
            return min(candidate_eqns, key=lambda e: len(filter_func(e.atoms(sympy.Symbol).difference(root))))

    def optimal_substitution_path(self, root, filter_func=filter_all):
        ''' Given a root, and a filter for the variables we want to substitute,
            find an optimal path through the equation graph

            >>> root = set(sympy.symbols('p1 p2 z23'))
            >>> for e in eqn_graph_full.optimal_substitution_path(root): print e
            p1 + q1 == 2*z12 + 1
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            p1 + p2*q2 + q1 + z24 + z34 == 2*z45 + 4*z46
            p2 + q2 + z35 + z45 == 2*z56 + 4*z57
            z46 + z56 + 1 == 2*z67
            z57 + z67 == 2*z78 + 1
            z78 == 0

            >>> root = set(sympy.symbols('p1 p2 z35 z34'))
            >>> for e in eqn_graph_full.optimal_substitution_path(root): print e
            p1 + q1 == 2*z12 + 1
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            p1 + p2*q2 + q1 + z24 + z34 == 2*z45 + 4*z46
            p2 + q2 + z35 + z45 == 2*z56 + 4*z57
            z46 + z56 + 1 == 2*z67
            z57 + z67 == 2*z78 + 1
            z78 == 0

            >>> root = set()
            >>> for e in eqn_graph_full.optimal_substitution_path(root,
            ...     filter_func=filter_carry_variables): print e
            z78 == 0
            z57 + z67 == 2*z78 + 1
            z46 + z56 + 1 == 2*z67
            p2 + q2 + z35 + z45 == 2*z56 + 4*z57
            p1*q2 + p2*q1 + z23 + 2 == 2*z34 + 4*z35 + 1
            p1 + p2*q2 + q1 + z24 + z34 == 2*z45 + 4*z46
            p1*q1 + p2 + q2 + z12 == 2*z23 + 4*z24 + 1
            p1 + q1 == 2*z12 + 1
        '''
        path = []

        next_eqn = self.smallest_equation_intersecting_root(root,
                                                            filter_func=filter_func)
        while next_eqn is not None:
            path.append(next_eqn)
            root.update(filter_func(next_eqn.atoms(sympy.Symbol)))
            # Now find everything with something extra to offer            
            next_eqn = self.smallest_equation_intersecting_root(root,
                                                                filter_func=filter_func)
            
            # And add all equations that are entirely contained in the current
            # roots but not added yet
            for eqn in self.equations.difference(set(path)):
                filtered_atoms = set(filter_func(eqn.atoms(sympy.Symbol)))
                if filtered_atoms.issubset(root):
                    path.append(eqn)

        return path

if __name__ == '__main__':
    import doctest
    from carry_equations_generator2 import generate_carry_equations
    prod = 143
    eqns = generate_carry_equations(product=prod)[:3]
    eqn_graph = EquationGraph(eqns)

    eqns_full = generate_carry_equations(product=prod)
    eqn_graph_full = EquationGraph(eqns_full)

    doctest.testmod(globs=globals())