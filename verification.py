# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 18:35:14 2014

@author: Richard
"""

import math

import GenerateTableOutput
import GenerateCarry
import EquationHandler
import sys
from sympy_solver import EquationSolver



def factorise(n):
    ''' Return list of factors '''
    factors = []
    i = 2
    while True:
        if n == 1:
            break
        while not n % i:
            factors.append(i)
            n /= i
        i += 1
    factors.reverse()
    return factors

def binary_factorisation(n):
    ''' Return the binary factorisation of a number '''
    factors = factorise(n)
    bin_fact = map(bin, factors)
    return bin_fact


def print_factorisation(n):
    ''' Print the binary factorisation of a number '''
    bin_fact = binary_factorisation(n)
    max_len = len(max(bin_fact, key=len))
    fact_str = '\n'.join([b_f[2:].rjust(max_len) for b_f in bin_fact])
    print '{}\n={}\n{}'.format(n, bin(n), fact_str)


def print_experiment_to_file(n, num_digits_factor):
    output_filename = 'results\{}.txt'.format(n)
    
    old_std = sys.stdout
    sys.stdout = open(output_filename, 'w')

    digitsInMultiplicand1 = num_digits_factor
    digitsInMultiplicand2 = num_digits_factor
    product = n
    
    if product < 10**16:
        print_factorisation(product)

    multiplier = []
    multiplication = []
    carry = []
    
    #digitsInMultiplicand1 will always have the greater number of digits
    if digitsInMultiplicand1 < digitsInMultiplicand2:
    	temp = digitsInMultiplicand1
    	digitsInMultiplicand1 = digitsInMultiplicand2
    	digitsInMultiplicand2 = temp
    
    binPrime = bin(product)[2:]
    if (digitsInMultiplicand1 + digitsInMultiplicand2) > len(binPrime):
    	for i in range (0, ((digitsInMultiplicand1 + digitsInMultiplicand2)-len(binPrime))):
    		binPrime = "0" + binPrime
    
    #Generate multipliers based on digits
    #	They take form 1,p2,p1,1 and 1,q2,q1,1
    #	This code will have to be rewritten to support >2 multiplicands
    strP = []
    strQ = []
    for i in range(1,digitsInMultiplicand1-1):
    	strP.append("p" + str(i))
    for i in range(1, digitsInMultiplicand2-1):
    	strQ.append("q" + str(i))
    strP.append("1")
    strP.insert(0, "1")
    strQ.append("1")
    strQ.insert(0, "1")
    multiplier.append(strP)
    multiplier.append(strQ)
    
    #Generate intermediary Multiplication row values
    #	This is the result of multiplying p by every bit of q
    for i in strQ:
    	temp = []
    	for j in strP:
    		if i == "1":
    			temp.append(j)
    		else:
    			if j == "1":
    				temp.append(i)
    			else:
    				temp.append(j + i)
    	multiplication.append(temp)
    
    #Find Carry row values
    myParams = [digitsInMultiplicand1, digitsInMultiplicand2]
    carry = GenerateCarry.CarryGenerator(myParams)
    
    #Generate Output
    myParams = [multiplier, multiplication, carry, binPrime]
    formattedCols = GenerateTableOutput.FormatOutput(myParams)
    #print ""
    
    #Generate Equations
    myParams = [formattedCols, carry]
    eqns = EquationHandler.GenerateEquations(myParams)
    
    system = EquationSolver.from_params(eqns)
    system.solve_equations(verbose=True)
    
    sys.stdout.close()
    sys.stdout = old_std

def do_batch(filename, num_digits_factor):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    for p in lines:
        print p
        p = int(p)
        try:
            print_experiment_to_file(p, num_digits_factor=num_digits_factor)
        except:
            print 'Error thrown for {}!!!'.format(p)

if __name__ == '__main__':
#    print_factorisation(143)
#    print_factorisation(70368895172689)
#    print_experiment_to_file(143, 4)
#    print_experiment_to_file(4299161663, 17)
    do_batch('17x17semiprimesFirst100.txt', 17)