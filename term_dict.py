# -*- coding: utf-8 -*-
"""
Created on Mon Jun 01 20:03:59 2015

@author: Richard
"""


from collections import defaultdict, MutableMapping
import operator

PRUNING = True

class TermDict(MutableMapping):
    ''' TermDict is a dictionary designed to represent a polynomial, allowing
        for constant lookup time of terms
    '''
    
    def __init__(self, default_type=int):
        self.default_type = default_type
        self._term_dict = defaultdict(default_type)
    
    def __delitem__(self, key):
        del self._term_dict[key]
    
    def __getitem__(self, key):
        value = self._term_dict[key]
        if PRUNING and (value == self.default_type()):
            self._term_dict.pop(key, None)
        return value

    def __setitem__(self, key, value):
        self._term_dict[key] = value
 
        if (PRUNING and (value == self.default_type())):
            self._term_dict.pop(key, None)
    
    def __iter__(self):
        return self._term_dict.__iter__()
    
    def __len__(self):
        return len(self._term_dict)
    
    def __repr__(self):
        return self._term_dict.__repr__()
        sorted_ = sorted(self._term_dict.items(), key=lambda x: len(str(x[0])))
        return '\n'.join(['{}: {}'.format(k, v) for k, v in sorted_])

    def copy(self):
        copy = TermDict(self.default_type)
        copy._term_dict = self._term_dict.copy()
        return copy


    ## Arithmetic functions
    def __additive_op__(self, other, operation):
        ''' Support addition 
        
            >>> td = TermDict()
            >>> td['x'] = 1
            >>> print td
            x: 1
            >>> td['x'] -= 1
            >>> print td
            <BLANKLINE>
            
            >>> td['x'] = td['y'] = td['xy'] = 1
            >>> print td
            y: 1
            x: 1
            xy: 1

            >>> td = td + 1
            >>> print td
            y: 2
            x: 2
            xy: 2

            >>> print td - 2
            <BLANKLINE>
            
            >>> print td + td
            y: 4
            x: 4
            xy: 4
        '''
        result = self.copy()

        if isinstance(other, (int, float)):
            for _key, _value in self.iteritems():
                result[_key] = operation(_value, other)
            return result

        if isinstance(other, TermDict):
            for _key, _value in other.iteritems():
                result[_key] = operation(self[_key], _value)
            return result

        raise ValueError('Unknown type of other {}: {}'.format(other, type(other)))
    
    def __add__(self, other):
        return self.__additive_op__(other, operator.add)
    
    def __sub__(self, other):
        return self.__additive_op__(other, operator.sub)


if __name__ == '__main__':
    import doctest
    doctest.testmod()