# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 11:14:47 2015

@author: Richard
"""

from collections import Iterable
import itertools
import sympy

from contradiction_exception import ContradictionException
from contradictions import apply_contradictions
from sympy_helper_fns import is_constant, is_one_or_zero, is_simple_binary, num_add_terms

PRUNING_DEFAULT = True

class EquivalenceDict(dict):
    ''' EquivalenceDict uses a graph like algorithm to store and retrieve
        equivalences.
        
        In this world everything maps to itself at the outset
    '''

    def __init__(self, *args, **kwargs):
        ''' Since the checks are bypassed when initialising with equivalences,
            we should do this manually
        '''
        super(EquivalenceDict, self).__init__()


        if 'pruning' in kwargs.keys():
            self.pruning = kwargs.pop('pruning')
        else:
            # This is the default value for pruning
            self.pruning = PRUNING_DEFAULT

        try:
            for arg in args:
                if isinstance(arg, dict):
                    for k, v in arg.iteritems():
                        self[k] = v
                elif isinstance(arg, Iterable):
                    for k, v in arg:
                        self[k] = v
                else:
                    raise NotImplementedError('Unknown type {}'.format(type(args)))
        except Exception as e:
            print args
            raise e
        
        for k, v in kwargs.iteritems():
            self[k] = v
        
    def __setitem__(self, key, value):
        ''' Set the roots of the graphs to be equal 
        
            >>> eq_dict = EquivalenceDict([(1, 3), (3, 6)])
            >>> eq_dict[2] = 4
            
            Nothing has crossed yet, so behave normally
            >>> print eq_dict
            {1: 3, 2: 4, 3: 6}
            
            Now the worlds collide, where 2 roots are equal
            >>> eq_dict[4] = 6
            >>> print eq_dict
            {1: 3, 2: 4, 3: 6, 4: 6}
            
            Now add another to the same tree.
            Note how the root 6 has been made to point to 7
            >>> eq_dict[1] = 7
            >>> print eq_dict
            {1: 3, 2: 4, 3: 6, 4: 6, 6: 7}
            
            Now assign 2 non roots. Note how the root of 1 graph, 11, points
            to the root of the other: 7.
            >>> eq_dict[10] = 11
            >>> eq_dict[10] = 2
            >>> print eq_dict
            {1: 3, 2: 4, 3: 6, 4: 6, 6: 7, 10: 11, 11: 7}
            
            Avoid cyclic things
            >>> eq_dict = EquivalenceDict()
            >>> for k, v in [(1, 2), (2, 1), (1, 2)]:
            ...     eq_dict[k] = v
        '''
        # Key is 'mapped' to value anyway, so don't put it in the underlying
        # data structure
        if key == value:
            return

        # Deal with the roots rather than the actual nodes
        #TODO Work out if we want overloaded getitem or that of the EquivalenceDict
        key = self[key]
        value = self[value]
#        key = EquivalenceDict.__getitem__(self, key)
#        value = EquivalenceDict.__getitem__(self, value)
        
        # If the roots are the same we're also done!        
        if key == value:
            return
        
        super(EquivalenceDict, self).__setitem__(key, value)
    
    def get_old(self, key):
        ''' Access to dict.get

            >>> eq_dict = EquivalenceDict([(1, 3), (3, 6), (6, 9)])
            >>> eq_dict.get_old(1)
            3
            >>> eq_dict.get_old(9)
            
            >>> eq_dict.get_old(5)
            
        '''
        return super(EquivalenceDict, self).get(key)
    
    def get(self, key):
        ''' Override get to include root-finding capabilities. In this world,
            where everything is equivalent to itself, get shouldn't return a
            None, rather, the key itself

            >>> eq_dict = EquivalenceDict([(1, 3), (3, 6), (6, 9)])
            >>> eq_dict.get(1)
            9
            >>> eq_dict.get(9)
            9
            >>> eq_dict.get(5)
            5
        '''
        return EquivalenceDict.__getitem__(self, key)
    
    def __getitem__(self, key):
        ''' Use the graph to find the root of the node
            Fetch the root of the key by fetching subsequent values found 

            NOTE keys that aren't found are still returned, as they are
            equivalent to themselves

            >>> eq_dict = EquivalenceDict([(1, 3), (3, 6), (6, 9)])
            >>> for i in [1, 3, 6, 9, 5]: print eq_dict[i]
            9
            9
            9
            9
            5
            
            Show off the fancy pruning
            >>> eq_dict = EquivalenceDict(pruning=True)
            >>> for k, v in [(1, 2), (4, 5), (1, 4), (1, 9), (1, 11)]:
            ...     eq_dict[k] = v
            ...     print eq_dict
            {1: 2}
            {1: 2, 4: 5}
            {1: 2, 2: 5, 4: 5}
            {1: 5, 2: 5, 4: 5, 5: 9}
            {1: 9, 2: 5, 4: 5, 5: 9, 9: 11}
        '''
        parent = super(EquivalenceDict, self).get(key)
        if parent is None:
            return key
        else:
            if self.pruning:
                super(EquivalenceDict, self).__setitem__(key, self.get(parent))
            # Recursive call to find the root
            # Make sure we call THIS function, not an overwritten version
            return EquivalenceDict.__getitem__(self, parent)


    def update(self, other):
        ''' Check update is using all the extra bells and whistles 
        
            >>> eqn1 = EquivalenceDict([(1, 3), (3, 6)])
            >>> eqn2 = EquivalenceDict([(1, 4)])
            >>> eqn1.update(eqn2)
            >>> print eqn1
            {1: 3, 3: 6, 6: 4}
        '''
        for k, v in other.iteritems():
            self[k] = v

    def copy(self):
        ''' Return a copy of the EquivalenceDict 
            
            >>> eqn = EquivalenceDict([(1,2)], pruning=True)
            >>> copy = eqn.copy()
            >>> print isinstance(copy, EquivalenceDict), copy.pruning
            True True
        '''
        copy = type(self)()
        copy.pruning = self.pruning
        for k, v in self.iteritems():
            copy[k] = v
        return copy
    
class BinaryEquivalenceDict(EquivalenceDict):
    ''' EquivalenceDict uses a graph like algorithm to get equivalences and
        check for consistency.
        
        In this world everything maps to itself at the outset
    '''
    # These are the states variables can be grounded in    
    GROUND_ROOTS = set([sympy.S.Zero, sympy.S.One])

    @staticmethod
    def _check_input_node(node):
        ''' Sympify inputs and check simple_binary (0, 1, x, 1-x) '''
        # Sympify any inputs
        if isinstance(node, int):
            node = sympy.sympify(node)
        if not is_simple_binary(node):
            raise ValueError('{} not allowed in a binary system'.format(node))
        return node

    def __getitem__(self, key):
        ''' Also check for 1-key

            >>> x, y, z = sympy.symbols('x y z')
            >>> eq_dict = BinaryEquivalenceDict([(x, y), (y, 1-z), (1-z, x)])
            >>> print eq_dict
            {x: y, y: -z + 1}
            >>> for i in [x, y, z, 1-x, 1-y, 1-z]: print '{} == {}'.format(i, eq_dict[i])
            x == -z + 1
            y == -z + 1
            z == z
            -x + 1 == z
            -y + 1 == z
            -z + 1 == -z + 1
            
            Show of all of the fancy logic
            >>> eq_dict[1-x] = 1
            >>> print eq_dict
            {x: y, z: 1, y: -z + 1}
            >>> for i in [x, y, z, 1-x, 1-y, 1-z]: print '{} == {}'.format(i, eq_dict[i])
            x == 0
            y == 0
            z == 1
            -x + 1 == 1
            -y + 1 == 1
            -z + 1 == 0

            Now that we check negations, avoid these infinite loops
            >>> eq_dict = BinaryEquivalenceDict({x: 1 - y, y: 1 - x})
            >>> print eq_dict[x]
            -y + 1

    
            Show off the fancy pruning
            >>> x1, x2, x3, x4 = sympy.symbols('x1 x2 x3 x4')
            >>> eq_dict = EquivalenceDict([(x1, x2), (x2, x3), (x3, x4)], pruning=True)
            >>> print eq_dict
            {x3: x4, x1: x2, x2: x3}
            >>> assert eq_dict[x2] == x4
            >>> print eq_dict
            {x3: x4, x1: x2, x2: x4}
            >>> assert eq_dict[x1] == x4
            >>> print eq_dict
            {x3: x4, x1: x4, x2: x4}
        '''
        key = self._check_input_node(key)
        
        # We only fetch things for single variables
        if num_add_terms(key) == 2:
            return 1 - BinaryEquivalenceDict.__getitem__(self, 1 - key)

        value = super(BinaryEquivalenceDict, self).__getitem__(key)

        if key == value:
            return value

        # If value isn't a ground state, maybe (not value) is.
        if value not in BinaryEquivalenceDict.GROUND_ROOTS:
            alt_value = super(BinaryEquivalenceDict, self).__getitem__(1 - value)
            return 1 - alt_value

        # If we can't do anything funky, never mind
        return value

    def __setitem__(self, key, value):
        ''' Set the roots of the graphs to be equal.
            Also:
            check for consistency
            shuffle the keys so the distinct roots are always roots
        
            Since we are performing the old method too, we need to repeat the
            tests with variable numbers instead
            >>> x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11 = sympy.symbols(
            ... ' '.join(['x{}'.format(i) for i in xrange(1, 12)]))
            >>> eq_dict = BinaryEquivalenceDict([(x1, x3), (x3, x6)])
            
            Check binary conditions are enforced            
            >>> eq_dict[x2] = 4
            Traceback (most recent call last):
                ...
            ValueError: 4 not allowed in a binary system
           
            
            Nothing has crossed yet, so behave normally
            >>> eq_dict[x2] = x4 
            >>> print eq_dict
            {x3: x6, x1: x3, x2: x4}
            
            Now the worlds collide, where 2 roots are equal
            >>> eq_dict[x4] = x6
            >>> print eq_dict
            {x3: x6, x4: x6, x1: x3, x2: x4}

            Now add another to the same tree.
            Note how the root 6 has been made to point to 7
            >>> eq_dict[x1] = x7
            >>> print eq_dict
            {x6: x7, x3: x6, x4: x6, x1: x3, x2: x4}
            
            Now assign 2 non roots. Note how the root of 1 graph, 11, points
            to the root of the other: 7.
            >>> eq_dict[x10] = x11
            >>> eq_dict[x10] = x2
            >>> print eq_dict
            {x3: x6, x10: x11, x2: x4, x6: x7, x11: x7, x4: x6, x1: x3}
            
            
            Now for the interesting binary/sympy related tests

            Check ground states are always roots            
            >>> eq_dict = BinaryEquivalenceDict([(0, x1)])
            >>> print eq_dict
            {x1: 0}
            >>> eq_dict[0] = x2
            >>> print eq_dict
            {x1: 0, x2: 0}
            
            Check non-monic terms are always roots.
            We want this so that we can substitute single variables out
            >>> eq_dict = BinaryEquivalenceDict([(1 - x1, x2)])
            >>> print eq_dict
            {x1: -x2 + 1}
            >>> eq_dict[0] = x2
            >>> print eq_dict
            {x1: -x2 + 1, x2: 0}
            >>> eq_dict[x3] = 1 - x1
            >>> print eq_dict
            {x3: 0, x1: -x2 + 1, x2: 0}
            >>> print eq_dict.items()
            [(x3, 0), (x1, 1), (x2, 0)]

            Check the roots of the above system
            >>> for expr in [0, 1, x1, x2, x3, 1-x1]: print expr, eq_dict[expr]
            0 0
            1 1
            x1 1
            x2 0
            x3 0
            -x1 + 1 0

            Make some obvious improvements
            >>> eq_dict = BinaryEquivalenceDict([(1 - x1, 1 - x2)])
            >>> print eq_dict
            {x1: x2}
            
            Check obvious contradiction            
            >>> eq_dict = BinaryEquivalenceDict()
            >>> eq_dict[0] = 1
            Traceback (most recent call last):
                ...
            ContradictionException: 0 != 1
            
            Check conflicting ground states
            >>> eq_dict = BinaryEquivalenceDict([(x1, 0), (x2, 1), (x1, x2)])
            Traceback (most recent call last):
                ...
            ContradictionException: 0 != 1

            Check conflicting variable choice
            >>> eq_dict = BinaryEquivalenceDict([(x1, x3), (x2, 1 - x3), (x1, x2)])
            Traceback (most recent call last):
                ...
            ContradictionException: x3 != -x3 + 1
            
            

        '''
        key = self._check_input_node(key)
        value = self._check_input_node(value)
        
        # Key is 'mapped' to value anyway, so don't put it in the underlying
        # data structure
        if key == value:
            return

        # Deal with the roots rather than the actual nodes
        key = BinaryEquivalenceDict.__getitem__(self, key)
        value = BinaryEquivalenceDict.__getitem__(self, value)

        # If roots are equal, they are already connected.
        if key == value:
            return

        # Catch inequality in variable space
        if key == 1 - value:
            raise ContradictionException('{} != {}'.format(key, value))

        # NOTE we already know key and value aren't equal
        if key in BinaryEquivalenceDict.GROUND_ROOTS:
            if value in BinaryEquivalenceDict.GROUND_ROOTS:
                # Now check that we haven't got a contradiction if we've grounded both
                # variables
                raise ContradictionException('{} != {}'.format(key, value))
            else:
                # If key is grounded, but value isn't then swap them over
                key, value = value, key
            
        # Double negate to we only ever get variable: equivalence
        if num_add_terms(key) == 2:
            key, value = 1 - key, 1 - value
        
        super(BinaryEquivalenceDict, self).__setitem__(key, value)

    def update(self, other):
        ''' Check update is using all the extra bells and whistles 
        
            >>> a, b, c, d = sympy.symbols('a b c d')            
            >>> eqn1 = BinaryEquivalenceDict([(a, b), (b, c)])
            >>> eqn2 = BinaryEquivalenceDict([(a, d)])
            >>> eqn1.update(eqn2)
            >>> print eqn1
            {c: d, b: c, a: b}
            
            Check they're all connected explicitly
            >>> for u, v in itertools.combinations([a, b, c, d], 2):
            ...     assert eqn1[u] == eqn1[v]
        '''
        super(BinaryEquivalenceDict, self).update(other)

    def copy(self):
        ''' Check copy method
    
            >>> eqn = BinaryEquivalenceDict([(1, sympy.Symbol('x'))], pruning=True)
            >>> copy = eqn.copy()
            >>> print isinstance(copy, BinaryEquivalenceDict), copy.pruning
            True True
        '''
        return super(BinaryEquivalenceDict, self).copy()

    ## Make sure iterables, keys, values and items all do the right things
    def __iter__(self):
        ''' Override __iter__ so that we can iterate over the implicit mappings
            in the usual way. Also test that the other methods also match this
            
            >>> a, b, c, d, e = sympy.symbols('a b c d e')
            >>> eq_dict = BinaryEquivalenceDict([(a, b), (c, 1-d), (d, 1-e)])

            >>> print eq_dict.items()
            [(c, e), (a, b), (d, -e + 1)]

            >>> for i in eq_dict: print i
            c
            a
            d
            
            We want the solutions to be as nice as possible
            >>> from solver_hybrid import SolverHybrid
            >>> from sympy_helper_fns import is_equation
            >>> from carry_equations_generator import generate_carry_equations
            >>> from semiprime_tools import num_to_factor_num_qubit

            >>> prod = 143
            >>> fact1, fact2 = num_to_factor_num_qubit(prod)
            >>> equations = generate_carry_equations(fact1, fact2, prod)
            >>> equations = filter(is_equation, equations)
            >>> system = SolverHybrid(equations)
            >>> system.solve_equations()
            >>> vars = sorted(system.variables.values(), key=str)
            >>> sorted_sol = sorted(system.solutions.items(), key=lambda x: str(x[0]))
            >>> for s in sorted_sol: print s
            (p1, q2)
            (p2, -q2 + 1)
            (q1, -q2 + 1)
            (z12, 0)
            (z23, 0)
            (z24, 0)
            (z34, 1)
            (z35, 0)
            (z45, 1)
            (z46, 0)
            (z56, 1)
            (z57, 0)
            (z67, 1)
            >>> for v in vars: print v, system.solutions[v]
            p1 q2
            p2 -q2 + 1
            q1 -q2 + 1
            q2 q2
            z12 0
            z23 0
            z24 0
            z34 1
            z35 0
            z45 1
            z46 0
            z56 1
            z57 0
            z67 1
            >>> for k in system.solutions: print k, system.solutions[k]
            z56 1
            z12 0
            q1 -q2 + 1
            z34 1
            z35 0
            z45 1
            z24 0
            z67 1
            z46 0
            p2 -q2 + 1
            z23 0
            z57 0
            p1 q2
        '''
        seen = set()
        for node in super(BinaryEquivalenceDict, self).__iter__():
            if node in self.GROUND_ROOTS:
                raise ValueError()
            node = node.atoms(sympy.Symbol)
            assert len(node) == 1
            node = node.pop()
            if node in seen:
                continue
            else:
                seen.add(node)
                yield node

    def iterkeys(self):
        return self.__iter__()

    def itervalues(self):
        ''' Override __iter__ so that we can iterate over the implicit mappings
            in the usual way. Also test that the other methods also match this
            
            >>> a, b, c, d, e = sympy.symbols('a b c d e')
            >>> eq_dict = BinaryEquivalenceDict([(a, b), (c, d), (d, 1-e), (e, 1)])

            >>> for i in eq_dict.itervalues(): print i
            0
            1
            b

            >>> for i in eq_dict.values(): print i
            0
            1
            b
        '''
        seen = set()
        seen_add = seen.add
        return (self[x] for x in self.iterkeys() if not (self[x] in seen or seen_add(self[x])))
    
    def iteritems(self):
        return ((k, self[k]) for k in self.iterkeys())

    def keys(self):
        return list(self.iterkeys())

    def values(self):
        return list(self.itervalues())

    def items(self):
        return list(self.iteritems())


if __name__ == "__main__":
    import doctest
    
    # Turn off pruning so the tests are clearer and stable
    PRUNING_DEFAULT = False
    doctest.testmod()