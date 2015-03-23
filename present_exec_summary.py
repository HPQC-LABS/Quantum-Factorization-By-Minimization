# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:01:47 2015

@author: Richard
"""
import itertools
import sympy
from sympy_solver import EquationSolver
from sympy_helper_fns import *
from verification import get_target_solutions
from cfg_sympy_solver import (EXPERIMENTS, QUBIT_REDUCTION_ID,
                              EXPERIMENTS_21)
from semiprime_tools import num_to_factor_num_qubit
from sympy_solver import EquationSolver
from verification import check_solutions, check_substitutions, get_target_factors, get_target_solutions
from carry_equations_generator import generate_carry_equations
from contradiction_exception import ContradictionException

prod = 365375409332725729551031220428655686168881402197
p, q = get_target_factors(prod)

assert p * q == prod
print '{} =\n\t{}x\n\t{}\n'.format(prod, p, q)
print '{} =\n\t{}x\n\t{}\n'.format(bin(prod)[2:], bin(p)[2:], bin(q)[2:])


eqns = '''q2 + q7 + z89 + 1 == 2*q2*q7 + 2*z910 + 4*z911
2*z910 + 4*z911 == z89 + 2
z810 + z910 + 2*z911 == 1'''.split('\n')
eqns = str_eqns_to_sympy_eqns(eqns)

## Work out the solutions by brute force
print 'Brute force Python solutions'
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


print len(valid_solns), 'solutions'
print variables
for s in valid_solns: print s


## Print the prolog solutions
print '\n'
print '''Prolog solutions:
Q2 = 0
Q7 = 1
Z810 = 0
Z89 = 0
Z910 = 1
Z911 = 0

Q2 = 1
Q7 = 0
Z810 = 0
Z89 = 0
Z910 = 1
Z911 = 0
'''



## Now get the solutions from the beginning starting with python
python_sol = get_target_solutions(prod)
for v in sorted(python_sol.keys(), key=str):
    print '{} = {}'.format(v, python_sol[v])

num_fact1, num_fact2 = num_to_factor_num_qubit(prod)
equations = generate_carry_equations(num_fact1, num_fact2, prod)
#from ut_generate_prolog_script import get_script, write_script
#write_script(get_script(equations), 'full_eqns_80x80_23.pl')
print '''Prolog solutions from all equations:
P1 = 0
P10 = 0
P11 = 0
P12 = 0
P13 = 0
P14 = 0
P15 = 0
P16 = 0
P17 = 0
P18 = 0
P19 = 0
P2 = 0
P20 = 0
P21 = 0
P22 = 0
P23 = 0
P24 = 0
P25 = 0
P26 = 0
P27 = 0
P28 = 0
P29 = 0
P3 = 1
P30 = 0
P31 = 0
P32 = 0
P33 = 0
P34 = 0
P35 = 0
P36 = 0
P37 = 0
P38 = 0
P39 = 0
P4 = 1
P40 = 0
P41 = 0
P42 = 0
P43 = 0
P44 = 0
P45 = 0
P46 = 0
P47 = 0
P48 = 0
P49 = 0
P5 = 0
P50 = 0
P51 = 0
P52 = 0
P53 = 0
P54 = 0
P55 = 0
P56 = 0
P57 = 0
P58 = 0
P59 = 0
P6 = 0
P60 = 0
P61 = 0
P62 = 0
P63 = 0
P64 = 0
P65 = 0
P66 = 0
P67 = 0
P68 = 0
P69 = 0
P7 = 1
P70 = 0
P71 = 0
P72 = 0
P73 = 0
P74 = 0
P75 = 0
P76 = 0
P77 = 0
P78 = 0
P8 = 0
P9 = 0
Q1 = 0
Q10 = 0
Q11 = 0
Q12 = 0
Q13 = 0
Q14 = 0
Q15 = 0
Q16 = 0
Q17 = 0
Q18 = 0
Q19 = 0
Q2 = 1
Q20 = 0
Q21 = 0
Q22 = 0
Q23 = 0
Q24 = 0
Q25 = 0
Q26 = 0
Q27 = 0
Q28 = 0
Q29 = 0
Q3 = 1
Q30 = 0
Q31 = 0
Q32 = 0
Q33 = 0
Q34 = 0
Q35 = 0
Q36 = 0
Q37 = 0
Q38 = 0
Q39 = 0
Q4 = 1
Q40 = 0
Q41 = 0
Q42 = 0
Q43 = 0
Q44 = 0
Q45 = 0
Q46 = 0
Q47 = 0
Q48 = 0
Q49 = 0
Q5 = 0
Q50 = 0
Q51 = 0
Q52 = 0
Q53 = 0
Q54 = 0
Q55 = 0
Q56 = 0
Q57 = 0
Q58 = 0
Q59 = 0
Q6 = 0
Q60 = 0
Q61 = 0
Q62 = 0
Q63 = 0
Q64 = 0
Q65 = 0
Q66 = 0
Q67 = 0
Q68 = 0
Q69 = 0
Q7 = 0
Q70 = 0
Q71 = 0
Q72 = 0
Q73 = 0
Q74 = 0
Q75 = 0
Q76 = 0
Q77 = 0
Q78 = 0
Q8 = 0
Q9 = 0
'''