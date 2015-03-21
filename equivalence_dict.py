# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 11:14:47 2015

@author: Richard
"""

from collections import Iterable
import sympy

from contradiction_exception import ContradictionException
from contradictions import apply_contradictions
from sympy_helper_fns import is_constant, is_one_or_zero, is_simple_binary

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

        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.iteritems():
                    self[k] = v
            elif isinstance(arg, Iterable):
                for k, v in arg:
                    self[k] = v
            else:
                raise NotImplementedError('Unknown type {}'.format(type(args)))
        
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
        '''
        # Key is 'mapped' to value anyway, so don't put it in the underlying
        # data structure
        if key == value:
            return
        # Deal with the roots rather than the actual nodes        
        key = self.get(key)
        value = self.get(value)
        
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
        ''' Override get to include root-finding capabilities
        
            >>> eq_dict = EquivalenceDict([(1, 3), (3, 6), (6, 9)])
            >>> eq_dict.get(1)
            9
            >>> eq_dict.get(9)
            9
            >>> eq_dict.get(5)
            5
        '''
        parent = super(EquivalenceDict, self).get(key)
        if parent is None:
            # It's not mapped to anything, so return 
            return key
        else:
            return self[parent]
    
    def __getitem__(self, key):
        ''' Use the graph to find the root of the node
            Fetch the root of the key by fetching subsequent values found 

            NOTE keys that aren't found are still returned, as they are
            equivalent to themselves

            >>> eq_dict = EquivalenceDict([(1, 3), (3, 6), (6, 9)])
            >>> eq_dict[1]
            9
            >>> eq_dict[3]
            9
            >>> eq_dict[6]
            9
            >>> eq_dict[9]
            9
            >>> eq_dict[5]
            5
        '''
        parent = super(EquivalenceDict, self).get(key)
        if parent is None:
            return key
        else:
            #TODO update to grandparents to flatten tree
#            super(EquivalenceDict, self).__setitem__(key)
            # Recursive call to find the root
            return self[parent]

class BinaryEquivalenceDict(EquivalenceDict):
    ''' EquivalenceDict uses a graph like algorithm to get equivalences and
        check for consistency.
        
        In this world everything maps to itself at the outset
    '''
    # These are the states variables can be grounded in    
    GROUND_ROOTS = set([sympy.S.Zero, sympy.S.One])

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
            >>> eq_dict = BinaryEquivalenceDict()
            >>> eq_dict[0] = 1
            Traceback (most recent call last):
                ...
            ContradictionException: 0 != 1
            
            # Check conflicting ground states
            >>> eq_dict = BinaryEquivalenceDict([(x1, 0), (x2, 1), (x1, x2)])
            Traceback (most recent call last):
                ...
            ContradictionException: 0 != 1

            # Check conflicting variable choice
            >>> eq_dict = BinaryEquivalenceDict([(x1, x3), (x2, 1 - x3), (x1, x2)])
            Traceback (most recent call last):
                ...
            ContradictionException: x3 != -x3 + 1
        '''
        if isinstance(key, int):
            key = sympy.sympify(key)
        if isinstance(value, int):
            value = sympy.sympify(value)
        for v in [key, value]:
            if is_constant(v) and (not is_one_or_zero(v)):
                raise ValueError('{} not allowed in a binary system'.format(v))
            assert is_simple_binary(v)

        # Key is 'mapped' to value anyway, so don't put it in the underlying
        # data structure
        if key == value:
            return

        # Deal with the roots rather than the actual nodes        
        key = self.get(key)
        value = self.get(value)
        
        # If roots are equal, they are already connected.
        if key == value:
            return

        # NOTE we already know key and value aren't equal
        if key in BinaryEquivalenceDict.GROUND_ROOTS:
            if value in BinaryEquivalenceDict.GROUND_ROOTS:
                # Now check that we haven't got a contradiction if we've grounded both
                # variables
                raise ContradictionException('{} != {}'.format(key, value))
            
            else:
                # We have that key is grounded and value isn't, so we should
                # swap them over so that ground states are always roots
                key, value = value, key
        
        apply_contradictions([sympy.Eq(key, value)])
        if key == 1 - value:
            raise ContradictionException('{} != {}'.format(key, value))
        
        super(BinaryEquivalenceDict, self).__setitem__(key, value)

if __name__ == "__main__":
    import doctest
    doctest.testmod()