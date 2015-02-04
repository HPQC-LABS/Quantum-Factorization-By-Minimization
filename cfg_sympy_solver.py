# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 13:31:21 2015

@author: Richard
"""

from schaller_equations_generator import generate_schaller_equations
from carry_equations_generator import generate_carry_equations

from objective_function_helper import (equations_to_vanilla_coef_str, 
                                       equations_to_sum_coef_str,
                                       equations_to_recursive_schaller_coef_str)

# Tuples of equation generator, objective function generator
QUBIT_REDUCTION_ID =   {0: (generate_carry_equations, equations_to_vanilla_coef_str),
                        1: (generate_schaller_equations, equations_to_sum_coef_str),
                        2: (generate_carry_equations, equations_to_recursive_schaller_coef_str)
                        }

EXPERIMENTS = {
    1: (4, 4, 143),
    2: (8, 8, 56153),
    3: (17, 17, 4299161663),
    4: (21, 21, 1099532599387),
    5: (24, 24, 70368895172689),
    6: (17, 17, 4296409109),
    7: (17, 17, 4306239659),
    8: (17, 17, 4345168637),
    9: (17, 17, 4314890543),
    10: (17, 17, 4297981997),
    11: (21, 21, 1099526307889),
    12: (17, 17, 4301127773),
    13: (45, 45, 309485009822787627980424653),
    14: (45, 45, 309485009821943203050291389),
    15: (51, 51, 1267650600228508624673600186743),
    100: (165, 165, 1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139),
}    
