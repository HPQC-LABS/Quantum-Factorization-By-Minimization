# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 13:31:21 2015

@author: Richard
"""

from collections import namedtuple

from schaller_equations_generator import generate_schaller_equations
from carry_equations_generator import generate_carry_equations

from objective_function_helper import (equations_to_vanilla_coef_str, 
                                       equations_to_sum_coef_str,
                                       equations_to_recursive_schaller_coef_str)

# File containing the list of semiprimes
SEMIPRIME_FILENAME = 'large_semiprimes.txt'
# File containing the answers
FACTORS_FILENAME = 'large_semiprimes_factors.txt'
# Name of the dict that verification.py will use to find the solutions
FACTOR_DICT_FILENAME = 'factor_dict'

# Tuples of equation generator, objective function generator
QUBIT_REDUCTION_ID =   {0: (generate_carry_equations, equations_to_vanilla_coef_str),
                        1: (generate_schaller_equations, equations_to_sum_coef_str),
                        2: (generate_carry_equations, equations_to_recursive_schaller_coef_str)
                        }

Experiment = namedtuple('Experiment', ('digits_multiplicand_1',
                                       'digits_multiplicand_2',
                                       'product',
                                       'num_qubits_expected'))

EXPERIMENTS = {
    1: Experiment(4, 4, 143, 1),
    2: Experiment(8, 8, 56153, 35),
    3: Experiment(17, 17, 4299161663, 1),
    4: Experiment(21, 21, 1099532599387, 1),
    5: Experiment(24, 24, 70368895172689, 0),
    6: Experiment(17, 17, 4296409109, 1),
    7: Experiment(17, 17, 4306239659, 1),
    8: Experiment(17, 17, 4345168637, 1),
    9: Experiment(17, 17, 4314890543, 1),
    10: Experiment(17, 17, 4297981997, 1),
    11: Experiment(17, 17, 4301127773, 1),
    12: Experiment(21, 21, 1099526307889, 0),
    13: Experiment(21, 21, 1099551473989, 1),
    14: Experiment(21, 21, 1099585029317, 1),
    15: Experiment(45, 45, 309485009822787627980424653, 1),
    16: Experiment(45, 45, 309485009821943203050291389, 1),
    17: Experiment(51, 51, 1267650600228508624673600186743, 1),
    18: Experiment(51, 51, 1267650600228402790082356974917, 1),
    100: Experiment(165, 165, 1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139, -1),
}    

PRESENTS = {
    #30x30_009
    1: Experiment(30, 30, 288230493189572837, 6),
    # 40x40_001
    2: Experiment(40, 40, 302231454937742154138497, 6),
    # 40x40_016
    3: Experiment(40, 40, 302231455142251316913389, 7),
    # 80x80_23
    4: Experiment(100, 100, 365375409332725729551031220428655686168881402197, 7),
}

EXPERIMENTS_21 = {
    1: Experiment(21, 21, 1099526307889, 0),
    2: Experiment(21, 21, 1099551473989, 1),
    3: Experiment(21, 21, 1099585029317, 1),
    4: Experiment(21, 21, 1099589223769, 0),
    5: Experiment(21, 21, 1099597612433, 1),
    6: Experiment(21, 21, 1099618585129, 0),
    7: Experiment(21, 21, 1099635361109, 1),
    8: Experiment(21, 21, 1099639557193, 1),
    9: Experiment(21, 21, 1099664721601, 1),
    10: Experiment(21, 21, 1099677306109, 1),
    11: Experiment(21, 21, 1099698279521, 1),
    12: Experiment(21, 21, 1099723448393, 1),
    13: Experiment(21, 21, 1099731839761, 0),
    14: Experiment(21, 21, 1099740228649, 5),
    15: Experiment(21, 21, 1099748617937, 1),
    16: Experiment(21, 21, 1099752803077, 1),
    17: Experiment(21, 21, 1099752812581, 1),
    18: Experiment(21, 21, 1099769592277, 1),
    19: Experiment(21, 21, 1099777982209, 21),
    20: Experiment(21, 21, 1099811527037, 1),
    21: Experiment(21, 21, 1099819917733, 1),
    22: Experiment(21, 21, 1099874457257, 1),
}

EXTENDED_EXPERIMENTS = EXPERIMENTS.copy()
EXTENDED_EXPERIMENTS.pop(100)
for exp_set in [EXPERIMENTS_21, ]:
    for exp_num, params in exp_set.iteritems():
        pre_digit = ''.join(map(str, params[:2]))
        extended_num = int(pre_digit + str(exp_num))
        assert EXTENDED_EXPERIMENTS.get(extended_num) is None
        EXTENDED_EXPERIMENTS[extended_num] = params