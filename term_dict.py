# -*- coding: utf-8 -*-
"""
Created on Mon Jun 01 20:03:59 2015

@author: Richard
"""


from collections import defaultdict, MutableMapping


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
        if value == self.default_type():
            self._term_dict.pop(key)
        return value

    def __setitem__(self, key, value):
        if value == self.default_type():
            if key in self.keys():
                self._term_dict.pop(key)
        else:
            self._term_dict[key] = value
    
    def __iter__(self):
        return self._term_dict.__iter__()
    
    def __len__(self):
        return len(self._term_dict)
    
    def __repr__(self):
        return self._term_dict.__repr__()
        
    def copy(self):
        copy = TermDict(self.default_type)
        copy._term_dict = self._term_dict.copy()
        return copy
    
if __name__ == '__main__':
    td = TermDict()
    td['x'] = 1
    print td
    td['x'] -= 1
    print td