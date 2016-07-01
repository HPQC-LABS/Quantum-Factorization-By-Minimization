# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 17:58:02 2015

    Tools for parsing results printouts, both of hybrid solvers and sage 
    outputs, and turning them into a summary spreadsheet, useful for plotting.

@author: Richard
"""

from collections import defaultdict
import os
import re


from semiprime_tools import factor_binary_differences, num_to_factor_num_qubit
from verification import VERIFICATION_FAILURE_MESSAGE, VERIFICATION_SUCCESS_MESSAGE
from sympy_helper_fns import str_eqns_to_sympy_eqns

grep_str = 'grep -r "End" * -A 7'

experiment_path_re = '(?P<dim1>\d*)x(?P<dim2>\d*)\/(?P=dim1)x(?P=dim2)_(?P<id_>\d*)_0.txt'
line_re = experiment_path_re + '(?P<delim>[-:])(?P<data>.*)'
line_re = re.compile(line_re)

LINE_START_TO_ATTRIBUTE = {
    'Num Qubits End:'
}

SAGE_TIMINGS_DIR = 'sage_timings'

CROP = slice(0, -3)

#t = experiment_path.match('100x100/100x100_047_0.txt:Num Qubits End: 1080')
#t = experiment_path.match('12x12')
#print t
#if t is not None:
#    print t.groups()

def _calc_prop0(prod):
    if isinstance(prod, (int, long)):
        prod = bin(prod)[2:]
    return '{:.3f}'.format(prod.count('0') * 1.0 / len(prod))

class ResultSummary(object):
    
    def __init__(self, dim1, dim2, id_):
        self.dim1 = dim1
        self.dim2 = dim2
        self.id_ = id_
        
        self.num_qubits_end = ''
        self.time_elapsed = ''
        self.verified = False
        self.hamming_distance = ''
        self.product = ''
        self.p = ''
        self.q = ''

        self.sage_time = '-1.0'
    
    def check_consistency(self):
        
        assert self.verified        
        
        if (self.p and self.q):
            
            p = int(self.p, 2)
            q = int(self.q, 2)            
            
            if self.hamming_distance:
                assert self.hamming_distance == str(factor_binary_differences(p, q))
        
            if self.product:
                assert p * q == int(self.product, 2)        
    
    def parse_data(self, data_blob):
        ''' Given a bit of text, work out what we want to do with it '''
        if data_blob.startswith('Num Qubits End:'):
            self.num_qubits_end = data_blob.split(' ')[-1]
        elif data_blob.startswith('Solved in '):
            self.time_elapsed = re.match('Solved in ([0-9\.]*)s', data_blob).group(1)
        elif data_blob == VERIFICATION_FAILURE_MESSAGE:
            raise ValueError(VERIFICATION_FAILURE_MESSAGE)
        elif data_blob == VERIFICATION_SUCCESS_MESSAGE:
            self.verified = True
        elif data_blob.startswith('Hamming distance: '):
            self.hamming_distance = re.match('Hamming distance: (\d*)', data_blob).group(1)
        elif re.match('[01]* =', data_blob):
            self.product = data_blob[:-2]
        elif re.match('[01]* x', data_blob):
            self.p = data_blob[:-2]
        elif re.match('[01]*', data_blob):
            self.q = data_blob
        elif data_blob == '':
            return
        else:
            raise ValueError('Unable to parse {}'.format(data_blob))
    
    @classmethod
    def from_text(cls, text):

        # Strip away any new line crud        
        text = text.strip()
        lines = text.split('\n')
        
        # Process the first line to assert the dimension and ID
        line = lines.pop(0)
        matched = line_re.match(line)
        
        # For the first line, do all the special stuff
        dim1, dim2, id_, delim, data = matched.groups()
        assert delim == ':'
        
        summary = cls(dim1, dim2, id_)
        summary.parse_data(data)
        
        for line in lines:
            _dim1, _dim2, _id_, _delim, data = line_re.match(line).groups()
            assert (dim1, dim2, id_, '-') == (_dim1, _dim2, _id_, _delim)
            summary.parse_data(data)
        
        summary.check_consistency()
        return summary
    
    @property
    def summary_tuple(self):
        return (self.dim1, self.dim2, self.id_, self.num_qubits_end, 
                self.time_elapsed, self.hamming_distance, self.p0_prop,
                self.q0_prop, _calc_prop0(int(self.product)), 
                SAGE_TIMINGS_DICT.get(str(int(self.product, 2)), '-1.0'), 
                self.product, self.p, self.q)[CROP]
    
    def to_text(self):
        pass

    def to_csv(self):
        return ','.join(self.summary_tuple)

    @property
    def p0_prop(self):
        return '{:.3f}'.format(1.0 * self.p.count('0') / len(self.p))
    
    @property
    def q0_prop(self):
        return '{:.3f}'.format(1.0 * self.q.count('0') / len(self.q))

def grep_to_result_summaries(filename):
    f = open(filename, 'r')
    text = f.read()
    exps = text.split('--')
    summaries = []
    for e in exps:
        res_sum = ResultSummary.from_text(e)
        summaries.append(res_sum)
    return summaries

def summaries_to_csv(summaries, filename):
    with open(filename, 'w') as f:
        header = ['dim1', 'dim2', 'id', 'num_qubits_end', 'time_elapsed',
                  'hamming_distance', 'p0_prop', 'q0_prop', 'prod0_prop', 'sage_time', 'product', 'p', 'q'][CROP]
        f.write(','.join(header + ['\n']))
        for summary in summaries:
            f.write(summary.to_csv())
            f.write('\n')
    return True

### Sage time extraction

def get_sage_timing_dict():
    files = os.listdir(SAGE_TIMINGS_DIR)
#    for f in files:
#        print f, type(f)
    sage_times = {}
    for filepath in files:
        if re.match('timings_\d*x\d*.txt', filepath):
            with open(os.path.join(SAGE_TIMINGS_DIR, filepath), 'r') as file_:
                for line in file_.readlines():
                    match = re.match('\((?P<time>\d*\.\d*),(?P<prod>\d*),.*\)', line)
                    if match:
                        data = match.groupdict()
                        sage_times[data['prod']] = data['time']
    return sage_times

SAGE_TIMINGS_DICT = get_sage_timing_dict()


### Equation extraction

def extract_equations(filename):
    with open(filename, 'r') as file_:
        eqns = []
        _add = False
        for line in file_.readlines():
            if line == 'Final Equations\n':
                _add = True
                continue
            if _add:
                if '==' in line:
                    eqns.append(line)
                else:
                    _add = False
    return str_eqns_to_sympy_eqns(eqns)

def extract_qubit_profile(filename):
    ''' 
        Extract the number of n-qubit interactions from a coefficient matrix file
    '''
    with open(filename, 'r') as file_:
        qubit_profile = defaultdict(int)
        for line in file_.readlines():
            num_qubits = len(line.strip().split(' ')) - 1            
            qubit_profile[num_qubits] += 1
    return qubit_profile
    

### Utilities
def recursive_result_search():
    ''' Return all files in a folder with .txt extension '''
    files = []
    for _ in os.walk('results'):
        for file_ in _[2]:
            if file_.endswith('.txt'):
                files.append(os.path.join(_[0], file_))
    return files

if __name__ == '__main__':
    input_ = 'grep_out.txt'
    output = 'summary_csv.txt'
    sums = grep_to_result_summaries(input_)
    summaries_to_csv(sums, output)

