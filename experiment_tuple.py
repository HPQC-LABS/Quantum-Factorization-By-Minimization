# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 12:03:22 2015

@author: Richard
"""

from collections import namedtuple

Experiment = namedtuple('Experiment', ('digits_multiplicand_1',
                                       'digits_multiplicand_2',
                                       'product',
                                       'num_qubits_expected'))

Experiment2 = namedtuple('Experiment2', ('p',
                                       'q',
                                       'product',
                                       'num_qubits_expected'))