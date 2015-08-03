# -*- coding: utf-8 -*-
"""
Created on Mon Aug 03 12:08:08 2015

    A home for misc functions and tools that have no other home
    
@author: Richard
"""

def unique_array_stable(array):
    ''' Given a list of things, return a new list with unique elements with
        original order preserved (by first occurence)
        
        >>> print unique_array_stable([1, 3, 5, 4, 7, 4, 2, 1, 9])
        [1, 3, 5, 4, 7, 2, 9]
    '''
    seen = set()
    seen_add = seen.add
    return [x for x in array if not (x in seen or seen_add(x))]