# -*- coding: utf-8 -*-
"""
Created on Fri Mar 06 10:42:11 2015

@author: Richard
"""

import os

from PrimeFactorizationEquationGenerator import num_to_factor_num_qubit

# Check the batch_scipts folder exists
BATCH_FOLDER = 'batch_scripts'
if not os.path.exists(BATCH_FOLDER):
    os.makedirs(BATCH_FOLDER)

input_filename = 'large_semiprimes.txt'
nums = open(input_filename, 'r').read()

# Turn the text into numbers
nums = map(int, nums.split('\n'))

# This holds the lines we want to print
lines = []

# Used to count the position of the prime in its category
count = 0
prev_num_digit = None

# f is the file object
f = None
for num in nums:
    num_d1, num_d2 = num_to_factor_num_qubit(num)
    if (num_d1, num_d2) != prev_num_digit:
        if lines:
            filename = 'batch_experiments_{}x{}.sh'.format(*prev_num_digit)
            filename = os.path.join(BATCH_FOLDER, filename)
            f = open(filename, 'w')
            f.write('#!/bin/bash\n\n')
            f.write('cd /home/tanburn/Quantum-Factorization-By-Minimization/\n\n')
            f.writelines(lines)
            f.close()
            lines = []
        count = 1
        prev_num_digit = (num_d1, num_d2)
        lines.append('\n\n## {}x{} semiprimes\n\n'.format(num_d1, num_d2))
    else:
        count += 1
    count_str = str(count).rjust(3, '0')
    entry_str = 'python ../PrimeFactorizationEquationGenerator.py {prod} > ../{num_d1}x{num_d2}/{num_d1}x{num_d2}_{count}_0.txt\n'.format(prod=num, num_d1=num_d1, num_d2=num_d2, count=count_str)
    lines.append(entry_str)
