# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 16:00:24 2015

@author: Richard
"""

from contradiction_exception import ContradictionException
from sympy_helper_fns import (expressions_to_variables, str_eqns_to_sympy_eqns,
                              is_equation)

eqns = '''q2 + q7 + z89 + 1 == 2*q2*q7 + 2*z910 + 4*z911
    2*z910 + 4*z911 == z89 + 2
    z810 == - z910 - 2*z911 + 1'''.split('\n')

eqns = str_eqns_to_sympy_eqns(eqns)

variables = list(expressions_to_variables(eqns))

valid_solns = []

for sol in itertools.product(range(2), repeat=len(variables)):
    to_sub = dict(zip(variables, sol))
    try:
        _eqns = [e.subs(to_sub) for e in eqns]
        _eqns = filter(is_equation, _eqns)
        if len(_eqns) > 0:
            raise ContradictionException(str(_eqns))
        
        valid_solns.append(sol)
    except ContradictionException as e:
        continue
    
print len(valid_solns)
print variables
for s in valid_solns: print s