#!/usr/bin/env python

import GenerateTableOutput
import GenerateCarry
import EquationHandler
import sys

from verification import print_binary_factorisation, check_solutions

__author__ = "Nathaniel Bryans"
__credits__ = ["Nathaniel Bryans", "Nikesh Dattani"]
__version__ = "0.0.6"
__status__ = "Prototype"

#Declare variables
#digitsInMultiplicand1 = 5
#digitsInMultiplicand2 = 5
#product = 493
#product = 899
#product = 841
#product = 551

#digitsInMultiplicand1 = 20
#digitsInMultiplicand2 = 20
#product = 1099511627775
#a = [37, 41, 43, 47, 53, 59, 61]
#b = [67, 71,73,79,83,89,97,101,107,109,113,127]
#m = []
#for i in a:
	#for j in b:
	#	m.append(i*j)

#for r in m:
digitsInMultiplicand1 = 7
digitsInMultiplicand2 = 6
#product = r
#product = 99
product = 3737

OutputFileName = "output.txt"

#digitsInMultiplicand1 = 5
#digitsInMultiplicand2 = 4
#product = 403

exp = 100
if exp == 1:
    digitsInMultiplicand1 = 4
    digitsInMultiplicand2 = 4
    product = 143
if exp == 2:
    digitsInMultiplicand1 = 8
    digitsInMultiplicand2 = 8
    product = 56153
if exp == 3:
    digitsInMultiplicand1 = 17
    digitsInMultiplicand2 = 17
    product = 4299161663

if exp == 6:
    digitsInMultiplicand1 = 17
    digitsInMultiplicand2 = 17
    product = 4296409109

if exp == 7:
    digitsInMultiplicand1 = 17
    digitsInMultiplicand2 = 17
    product = 4306239659

if exp == 8:
    digitsInMultiplicand1 = 17
    digitsInMultiplicand2 = 17
    product = 4345168637

if exp == 9:
    digitsInMultiplicand1 = 17
    digitsInMultiplicand2 = 17
    product = 4314890543

if exp == 10:
    # 3 term test
    digitsInMultiplicand1 = 17
    digitsInMultiplicand2 = 17
    product = 4297981997

if exp == 4:
    digitsInMultiplicand1 = 21
    digitsInMultiplicand2 = 21
    product = 1099532599387

if exp == 11:
    # 3 term test 2
    digitsInMultiplicand1 = 21
    digitsInMultiplicand2 = 21
    product = 1099526307889

if exp == 5:
    digitsInMultiplicand1 = 24
    digitsInMultiplicand2 = 24
    product = 70368895172689

if exp == 12:
    digitsInMultiplicand1 = 17
    digitsInMultiplicand2 = 17
    product = 4301127773

if exp == 13:
    digitsInMultiplicand1 = 45
    digitsInMultiplicand2 = 45
    product = 309485009822787627980424653

if exp == 14:
    digitsInMultiplicand1 = 45
    digitsInMultiplicand2 = 45
    product = 309485009821943203050291389

if exp == 15:
    digitsInMultiplicand1 = 51
    digitsInMultiplicand2 = 51
    product = 1267650600228508624673600186743


if exp == 100:
    digitsInMultiplicand1 = 165
    digitsInMultiplicand2 = 165
    product = 1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139


#We can override the digit and product values above using arguments
if len(sys.argv) > 2:
	digitsInMultiplicand1 = int(sys.argv[1])
	digitsInMultiplicand2 = int(sys.argv[2])
	product = int(sys.argv[3])
	OutputFileName = str(sys.argv[4])

numMultipliers = 2	#eventually this will be dynamically generated by output

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

print binPrime

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

if 0:
    import EquationSolver
    myParams = [eqns, carry, strP, strQ, product, OutputFileName]
    EquationSolver.SolveEquation(myParams)
else:
    from sympy_solver import EquationSolver
    from time import time
    s = time()
    # None means the result will be printed to screen
    output = None#OutputFileName
    system = EquationSolver.from_params(eqns, output_filename=output, 
                                        log_deductions=False)
    system.solve_equations(verbose=True)
    system.print_summary()
    print 'Solved in {}s'.format(time() - s)
#    try:
#        coef_filename = OutputFileName.replace('.txt', '_coef.txt')
#        system.objective_function_to_file(coef_filename)
#    except:
#        print 'Failed to write the coefficient'

    check_solutions(product, system.solutions.copy(), verbose=True)

#    sim1 = system.simplified_system()
#    sim1.update_value(sim1.get_var('p1'), 1)
#    sim1.solve_equations(verbose=True)
