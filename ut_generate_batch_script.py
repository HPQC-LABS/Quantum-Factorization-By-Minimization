# -*- coding: utf-8 -*-
"""
Created on Fri Mar 06 10:42:11 2015

@author: Richard
"""

import os

from semiprime_tools import num_to_factor_num_qubit
from cfg_sympy_solver import SEMIPRIME_FILENAME, SEMIPRIMES_HAMMING_TEMPLATE

BATCH_FOLDER = 'batch_scripts'
BATCH_FOLDER_HAMMING = 'batch_scripts_hamming'
RESULTS_FOLDER = 'results'
input_filename = SEMIPRIME_FILENAME


DIMENSIONS = range(20, 260, 10)

BATCH_TEMPLATE = 'batch_experiments_{}x{}.sh'
BATCH_TEMPLATE_HAMMING = 'batch_experiments_{}x{}.sh'

RUN_BATCH_FILENAME = 'run_batch.sh'
RUN_BATCH_HAMMING_FILENAME = 'run_batch_hamming.sh'
PRINT_RESULTS_FILENAME = 'print_results.sh'
COPY_RESULTS_FILENAME = 'copy_results.sh'

BASH_STR = '#!/bin/sh'

HOME_DIR = '/home/tanburn/Quantum-Factorization-By-Minimization/'
BATCH_DIR = os.path.join(HOME_DIR, BATCH_FOLDER)
BATCH_HAMMING_DIR = os.path.join(HOME_DIR, BATCH_FOLDER_HAMMING)
RESULTS_DIR = os.path.join(HOME_DIR, RESULTS_FOLDER)

def generate_batch_scripts_hamming(dim='20'):
    dim = str(dim)
    if not os.path.exists(BATCH_FOLDER_HAMMING):
        os.makedirs(BATCH_FOLDER_HAMMING)
    
    filename = SEMIPRIMES_HAMMING_TEMPLATE.format(dim, dim)
    exps = open(filename, 'r').read()

    # Turn the text into tuples
    exps = exps.split('\n')
    exps = filter(None, exps)

    # Used to count the position of the prime in its category
    count = 0

    batch_filename = BATCH_TEMPLATE_HAMMING.format(dim, dim)
    batch_filename = os.path.join(BATCH_FOLDER_HAMMING, batch_filename)
    with open(batch_filename, 'w') as f:    
        f.write(BASH_STR + '\n\n')
        f.write('cd {HOME_DIR}\n\n'.format(HOME_DIR=HOME_DIR))

        for exp in exps:
            prod, dist, p, q = exp.strip('(').strip(')').split(', ')
            count += 1
            count_str = str(count).rjust(3, '0')
            entry_str = 'python PrimeFactorizationEquationGenerator.py {prod} > {RES_FOL}/{num_d1}x{num_d2}/{num_d1}x{num_d2}_{count}_0.txt\n'.format(prod=prod, num_d1=dim, num_d2=dim, count=count_str, RES_FOL=RESULTS_DIR)
            f.write(entry_str)
    f.close()

    return True

def generate_batch_scripts():
# Check the batch_scipts folder exists
    if not os.path.exists(BATCH_FOLDER):
        os.makedirs(BATCH_FOLDER)
    
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
                filename = BATCH_TEMPLATE.format(*prev_num_digit)
                filename = os.path.join(BATCH_FOLDER, filename)
                f = open(filename, 'w')
                f.write(BASH_STR + '\n\n')
                f.write('cd {HOME_DIR}\n\n'.format(HOME_DIR=HOME_DIR))
                f.writelines(lines)
                f.close()
                lines = []
            count = 1
            prev_num_digit = (num_d1, num_d2)
            lines.append('\n\n## {}x{} semiprimes\n\n'.format(num_d1, num_d2))
        else:
            count += 1
        count_str = str(count).rjust(3, '0')
        entry_str = 'python PrimeFactorizationEquationGenerator.py {prod} > {RES_FOL}/{num_d1}x{num_d2}/{num_d1}x{num_d2}_{count}_0.txt\n'.format(prod=num, num_d1=num_d1, num_d2=num_d2, count=count_str, RES_FOL=RESULTS_DIR)
        lines.append(entry_str)

    return True

def generate_run_batch_scripts(output_filename=RUN_BATCH_FILENAME,
                               batch_folder=BATCH_DIR, 
                               batch_template=BATCH_TEMPLATE):
    ''' Generate a file that will open a new screen and run a given batch file 
    
    '''
    f = open(output_filename, 'w')
    f.write(BASH_STR + '\n\n')

    f.write('cd {BATCH_DIR}\n\n'.format(BATCH_DIR=batch_folder))

    template = batch_template.format('$1', '$1')

    script = '''
if [ "$1" == "" ]; then
    exit 1
fi

EXP_NAME="exp$1"

screen -S $EXP_NAME -d -m

BATCH_FILE="batch_experiments_$1x$1.sh"

screen -S $EXP_NAME -X stuff ". ./$BATCH_FILE\\r"'''.format(temp=template)

    f.write(script)
    
    f.close()

def generate_run_batch_hamming_scripts(output_filename=RUN_BATCH_HAMMING_FILENAME,
                               batch_folder=BATCH_HAMMING_DIR, 
                               batch_template=BATCH_TEMPLATE_HAMMING):
    ''' Generate a file that will open a new screen and run a given batch file 
    
    '''
    f = open(output_filename, 'w')
    f.write(BASH_STR + '\n\n')

    f.write('cd {BATCH_DIR}\n\n'.format(BATCH_DIR=batch_folder))

    template = batch_template.format('$1', '$1')

    script = '''
if [ "$1" == "" ]; then
    exit 1
fi

EXP_NAME="hexp$1"

screen -S $EXP_NAME -d -m

BATCH_FILE="batch_experiments_$1x$1.sh"

screen -S $EXP_NAME -X stuff ". ./$BATCH_FILE\\r"'''.format(temp=template)

    f.write(script)
    
    f.close()

def generate_print_results(dimensions):
    ''' Generate the script that puts the grep output into a file in a given folder '''
    f = open(PRINT_RESULTS_FILENAME, 'w')
    f.write(BASH_STR + '\n\n')

    f.write('cd {HOME_DIR}\n\n'.format(HOME_DIR=HOME_DIR))

    f.write('SEARCH_STRING="End"\n\n')

    template = '''echo "*** {dim} ***"\necho\ngrep $SEARCH_STRING {RES_FOL}/{dim}x{dim}/*\necho\n'''
    for dim in dimensions:
        f.write(template.format(dim=dim, RES_FOL=RESULTS_FOLDER))
    f.close()

def generate_copy_results(dimensions):
    ''' Generate the script that copies the grep output into 
        the public_html folder    
    '''
    f = open(COPY_RESULTS_FILENAME, 'w')
    f.write(BASH_STR + '\n\n')
    f.write('cd {HOME_DIR}\n\n'.format(HOME_DIR=HOME_DIR))

    template = '''cp -r {RES_FOL}/{dim}x{dim}/ ../public_html/\n'''
    for dim in dimensions:
        f.write(template.format(dim=dim, RES_FOL=RESULTS_FOLDER))

    f.write('./{print_res_file} > ../public_html/results_summary.txt\n'.format(print_res_file=PRINT_RESULTS_FILENAME))
    
    f.close()


if __name__ == '__main__':
    generate_batch_scripts()
    generate_run_batch_scripts()
    generate_run_batch_hamming_scripts()
    generate_print_results(DIMENSIONS)
    generate_copy_results(DIMENSIONS)
    for dim in xrange(20, 310, 10):
        generate_batch_scripts_hamming(dim)
    