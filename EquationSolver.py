#!/usr/bin/env python

from EquationGraphNode import Node
from EquationSolverHelpers import EquationSolverState
from EquationSolverHelpers import VarValueAssignmentPicker
import ReHandler
import Tests
import sys
import pdb
import itertools
from copy import deepcopy
from copy import copy

__author__ = "Nathaniel Bryans"
__credits__ = ["Nathaniel Bryans", "Nikesh Dattani"]
__version__ = "0.0.91"
__status__ = "Prototype"

contradiction = 0	#changed to 1 when we encounter a contradiction

def GenerateMultiplicationTuple(params):
	nodeleft = Node(params[0][:params[1]])
	noderight = Node(params[0][params[1]:])
	nodeparent = Node("x")
	nodeparent.AddChild(nodeleft)
	nodeparent.AddChild(noderight)
	return nodeparent

def SolveEquation(params):
	#params[0][0][i] contains all leftsides
	#params[0][1][i] contains all rightsides
	#params[1][i] contains carrys
	#params[2] contains strP
	#params[3] contains strQ
	#params[4] contains product
	#params[5] contain OutputFileName
	
	columnEqns = []
	zValues = {}
	varValues = {}
	eqnValues = {}
	J1sApplied = []
	
	strP = params[2]
	strQ = params[3]
	OutputFileName = params[5]
	
	global contradiction
	
	#Build Equation Tree
	for i in range(0, len(params[0][0])):
		#for each eqn create leftside
		leftnode = Node("+")
		for j in params[0][0][i]:
			if j:			
				tempnode = Node(j)
				leftnode.AddChild(tempnode)
		
		#create rightside
		rightnode = Node("+")
		for j in params[0][1][i]:
			if j:
				tuple = ReHandler.MatchProductWithNatural(j)
				if tuple:
					myParams = [j, tuple]
					tempnode = GenerateMultiplicationTuple(myParams)
					rightnode.AddChild(tempnode)
				else:
					tempnode = Node(j)
					rightnode.AddChild(tempnode)
				
		#append together
		parentNode = Node("=")
		parentNode.AddLeft(leftnode)
		parentNode.AddRight(rightnode)
		
		#add to eqn collection
		columnEqns.append(parentNode)
	
	
	for i in params[1]:
		for j in i:
			keyStr = "z" + str(j[0]) + str(j[1])
			zValues[keyStr] = -1	

	
	#Tests.TestModifiedJ1Test3(varValues, eqnValues, J1sApplied)
	#Tests.TestJ2Test2(zValues)
	#Tests.TestJ4Test3(varValues, zValues)
	#Tests.TestJ5(varValues, strP, strQ) 
	#Tests.TestJ6(varValues)
	#Tests.TestJ7Test12(varValues, zValues)
	
	#exit()
	#print "press enter to continue:",
	#filename = raw_input()
	
	solved = 0
	depth = params[5]
	
	for i in columnEqns:
		hz = HandleZero(i, zValues, varValues)
		if hz == 1:
			columnEqns.remove(i)
	
	columnEqnsBeginningStr = []
	columnEqnsBeginning = deepcopy(columnEqns)
	for i in range(0, len(columnEqns)):
		columnEqnsBeginningStr.append(columnEqns[i].PrintGraph())

	X = RunAndCheck(zValues, varValues, eqnValues, J1sApplied, columnEqns, columnEqnsBeginningStr, strP, strQ)

	valuesChecked = ['', '']
	
	if X == 1:
		solved = 1 #Managed to solve without substituting any variables
		print "Solved"
	
	
	remainingQubits = CountQubitsRemaining(zValues, varValues, strP, strQ)
	possibleQubits = MaximumQubits(zValues, strP, strQ)
	#if (remainingQubits < possibleQubits/4) | (remainingQubits < 30):
	#if (remainingQubits < possibleQubits/4):
	f = open(OutputFileName, 'a')
	#f.write('-----------------------------------------\n')
	f.write(str(params[4]) + ' ')
	f.write('possibleQubit: ' + str(possibleQubits) + ' remainingQubits: ' + str(remainingQubits) + '\n')
	for i in columnEqns:
		f.write(i.PrintGraph())
		f.write('\n')
	f.write("Z values:\n")
	for i in zValues.keys():
		f.write(str(i) + " " + str(zValues[i])+"\n")
	#f.write("varValues:\n")
	#f.write(varValues+"\n")
	f.close()
	
	

		
def RunAndCheck(zValues, varValues, eqnValues, J1sApplied, columnEqns, columnEqnsStart, strP, strQ):
	global contradiction
	count = 0
	consecutiveInaction = 0
	solved = 0
	done = 0
	count = 0
	while solved == 0 and done == 0:
		count += 1
		StepsToSolveEquation(zValues, varValues, eqnValues, J1sApplied, 
			columnEqns, strP, strQ)
		
		if CompareOldAndNewEqns(columnEqnsStart, columnEqns) > 0:
			consecutiveInaction = 0
			columnEqnsStart = []
			for i in range(0, len(columnEqns)):
				columnEqnsStart.append(columnEqns[i].PrintGraph())
		else:
			consecutiveInaction += 1
			
			#Check if done iteration
		solved = 1	
		for i in zValues.keys():
			if zValues[i] == -1:	#There is still a zValue to cancel out
				solved = 0	
		if solved == 0:
			solved = 1
			for i in range(1, len(strP)-1):
				if strP[i] not in varValues:	#There is a value of p remaining to find
					solved = 0
			for i in range(1, len(strQ)-1):
				if strQ[i] not in varValues:	#There is a value of q remaining to find
					solved = 0
		if contradiction == 1:
			solved = 0
			done = 1
		
		if consecutiveInaction > 1:
			done = 1
		count += 1
		
	if solved == 1:
		#Run the steps to solve equation ONE more time to look for contradiction
		StepsToSolveEquation(zValues, varValues, eqnValues, J1sApplied, 
			columnEqns, strP, strQ)
		if contradiction == 1:
			#There was a contradiction after all
			contradiction = 0
			return 0
		#print "Number of times running StepsToSolveEquation: " + str(count)
		return 1
	if done == 1:
		if contradiction == 1:
			contradiction = 0
		return 0
			

def SetNextLowestUnassignedZValue(zValues, assignedValue, curLowestToSet):
	#print "here"
	for i in sorted(zValues):
		if zValues[i] == -1:
			curLowestToSet = curLowestToSet - 1
			if curLowestToSet == 0:
				zValues[i] = assignedValue
				#print i
				return 1 
			elif curLowestToSet == -1:
				return 1
	return 0
	
def StepsToSolveEquation2(zValues, varValues, eqnValues, J1sApplied, 
	columnEqns, strP, strQ):
	global contradiction
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	for i in columnEqns:
		ModifiedJudgement1(i, varValues, eqnValues, J1sApplied)
		Judgement2(i, zValues)
		Judgement3(i, zValues)
		Judgement4(i, varValues, zValues)
		Judgement7(i, varValues, zValues)
	Judgement5(varValues, strP, strQ)
	Judgement6(varValues)
	Judgement8(varValues, strP, strQ)
	CheckForVarAssignmentContradictions(varValues, strP, strQ)
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	# for i in columnEqns:
		# Judgement7(i, varValues, zValues)
	# SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	# CheckForVarAssignmentContradictions(varValues, strP, strQ)
	# for i in columnEqns:
		# Judgement7(i, varValues, zValues)
	# SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	# CheckForVarAssignmentContradictions(varValues, strP, strQ)
	
def StepsToSolveEquation(zValues, varValues, eqnValues, J1sApplied, 
	columnEqns, strP, strQ):	
	global contradiction
	
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	
	#Judgement 1
	for i in columnEqns:
		ModifiedJudgement1(i, varValues, eqnValues, J1sApplied)
		
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	
	#Judgement 2
	for i in columnEqns:
		Judgement2(i, zValues)
		
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	
	#Judgement 3
	for i in columnEqns:
		Judgement3(i, zValues)
		
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	
	#Judgement 4
	for i in columnEqns:
		Judgement4(i, varValues, zValues)

	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	
	#Judgement 5
	Judgement5(varValues, strP, strQ)
	
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)	
	
	#Judgement 6
	Judgement6(varValues)
	
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	
	#Judgement 8
	Judgement8(varValues, strP, strQ)
	
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	
	#Judgement 7
	for i in columnEqns:
		Judgement7(i, varValues, zValues)
		
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	
	for i in columnEqns:
		Judgement7(i, varValues, zValues)
		
	SubstitueVariables(columnEqns, zValues, varValues, eqnValues)
	
	for i in columnEqns:
		Judgement7(i, varValues, zValues)
		
	CheckForVarAssignmentContradictions(varValues, strP, strQ)
	
def PrintEquations(columnEqns):
	print ""
	for i in columnEqns:
		print i.PrintGraph()
	print ""
	
def CompareOldAndNewEqns(list1, list2):	
	strings1 = []
	strings2 = []
	
	for i in list1:
		strings1.append(i)
	for j in list2:
		strings2.append(j.PrintGraph())
		
	if len(strings1) != len(strings2):
		return 0
		
	for i in strings1:
		match = 0
		for j in strings2:
			if i == j:
				match = 1
				break
		if match == 0:
			return 1
	return 0

def SubstitueVariables(columnEqns, zValues, varValues, eqnValues):
	eqnsToDelete = []
	for i in columnEqns:
		SubstituteVarsForValues(i, zValues, varValues, eqnValues)
		
		IntermediatePassRules(i)
		
		x = HandleZero(i, zValues, varValues)
		if x == 1:
			columnEqns.remove(i)
			

def HandleZero(eqn, zValues, varValues):
	global contradiction
	if eqn.ContainsLoneZeroOnRight() == 1:		
		leftterms = eqn.left.GetTerms()
		for i in leftterms:
			if ReHandler.IsInt(i) > 0 and int(i) != 0:
				#print "FOUND A CONTRADICTION 1"
				contradiction = 1
			else:
				#if its a z value set z
				if ReHandler.MatchZVariable(i):
					SetZValueAndCheckForContradiction(zValues, i, 0)
				else:
					SetVarValueAndCheckForContradiction(varValues, str(i), 0)
		return 0
	eqn.RemoveZeroes()
	#SHOULD DELETE EQN, not removeZeroes
	return 0
		
def ModifiedJudgement1(eqn, varValues, eqnValues, J1sApplied):
	global contradiction
	numTermsOnRight = eqn.right.CountTerms()
	numTermsOnLeft = eqn.left.CountTerms()
	sumOfNonVariablesOnLeft = eqn.left.CountNonVariables()
	sumOfNonVariablesOnRight = eqn.right.CountNonVariables()
	
	if sumOfNonVariablesOnLeft == 0 and numTermsOnRight == 1 and sumOfNonVariablesOnRight == 1 and numTermsOnLeft >= 2:
		#find leftterms
		leftterms = eqn.left.GetTerms()
		eqnTuple = tuple(leftterms)
		eqnValues[eqnTuple] = 1
		for i in leftterms:
			for j in leftterms:
				if i != j:
					string1 = str(i) + str(j)
					string2 = str(j) + str(i)
					if string1 not in J1sApplied and string2 not in J1sApplied:	
						#print "Applying J1 to equation " + eqn.PrintGraph()
						J1sApplied.append(string1)
						J1sApplied.append(string2)
						for i in varValues.keys():
							if str(i) == str(string1):
								if varValues[string1] == 1:
									#print "WE FOUND A CONTRADICTION 9"
									contradiction = 1
							if str(i) == str(string2):
								if varValues[string2] == 1:
									#print "WE FOUND A CONTRADICTION 9"
									contradiction = 1
						if contradiction == 0:
							SetVarValueAndCheckForContradiction(varValues, string1, 0)
							SetVarValueAndCheckForContradiction(varValues, string2, 0)

	elif sumOfNonVariablesOnLeft == 0 and numTermsOnRight == 1 and sumOfNonVariablesOnRight == 1 and numTermsOnLeft == 1:
		#eqn simply says one variable equals 1, set that variable to 1
		leftterms = eqn.left.GetTerms()
		string1 = str(leftterms[0])
		#if string1 not in J1sApplied:
		J1sApplied.append(string1)
		SetVarValueAndCheckForContradiction(varValues, string1, 1)
			#Remove equation?

def Judgement2(eqn, zValueDict):
	if eqn.value != "=":
		#error
		print "ERROR in judgement 2"
		exit()
	
	numTermOnLeft = eqn.left.CountTermsContainingVariables()
	sumOfNonVariablesOnRight = eqn.right.CountNonVariables()
	sumOfNonVariablesOnLeft = eqn.left.CountNonVariables()
	largestMultiplier = eqn.right.FindLargestMultiplier(-1, 0)
	
	while largestMultiplier > 0:
		#perform check
		if largestMultiplier > numTermOnLeft - (sumOfNonVariablesOnRight-sumOfNonVariablesOnLeft):
			#set associated z variable to 0
			#remove z variable and multiplier from graph
			#print "Applying J2 to equation " + eqn.PrintGraph()
			x = eqn.right.RetrieveAndDeleteTermWithLargestMultiplier(largestMultiplier)
			for i in x:
				tuple = ReHandler.MatchZVariable(i)
				if tuple:
					SetZValueAndCheckForContradiction(zValueDict, tuple, 0)
					
		#calculate next largest multiplier
		largestMultiplier = eqn.right.FindLargestMultiplier(largestMultiplier, 0)
		

def Judgement3(eqn, zValueDict):
	#Depending on our implementation of J3 it may be possible to eliminate '-'
	#signs in the equation trees by simply not moving terms that would cause a -
	#from left to right UNTIL we can apply J3
	if eqn.value != "=":
		#error
		print "ERROR in judgement 3"
		exit()
	
	#if a single z term on right with a natural number multiplier
	#if sum of NonVariables on right minus NonVariables from left is negative
	sumOfNonVariablesOnRight = eqn.right.CountNonVariables()
	sumOfNonVariablesOnLeft = eqn.left.CountNonVariables()
	numberOfZTermsOnRight = eqn.right.CountTermsContainingVariables()
	
	if numberOfZTermsOnRight == 1:
		if sumOfNonVariablesOnRight-sumOfNonVariablesOnLeft < 0:			
			largestMultiplier = eqn.right.FindLargestMultiplier(-1, 0)
			if largestMultiplier > 1:
				#We have met the conditions for J3
				#print "Applying J3 to equation " + eqn.PrintGraph()
				#Have to move terms from left to right
				#Have to set variable to be 1
				#Have to have terms on right cancel out
				#MAY have to add a minus sign here

				x = eqn.right.RetrieveAndDeleteTermWithLargestMultiplier(largestMultiplier)
				for i in x:
					tuple = ReHandler.MatchZVariable(i)
					if tuple:
						SetZValueAndCheckForContradiction(zValueDict, tuple, 1)
					else:
						tempnode = Node(i)
						eqn.right.AddChild(tempnode)

				eqn.BalanceNaturalNumbers()

# if x + y + .. = a where n (number of terms on left) = a, then x=1, y=1, ...
def Judgement4(eqn, varValues, zValueDict):
	global contradiction
	sumOfNonVariablesOnLeft = eqn.left.CountNonVariables()
	sumOfNonVariablesOnRight = eqn.right.CountNonVariables()
	numberOfTermsOnRightContainingVariables = eqn.right.CountTermsContainingVariables()
	numberOfTermsOnLeftContainingVariables = eqn.left.CountTermsContainingVariables()
	numberOfTermsOnLeft = eqn.left.CountTerms()
	numberOfTermsOnRight = eqn.right.CountTerms()
	
	if numberOfTermsOnRightContainingVariables == 0:
		if sumOfNonVariablesOnLeft + numberOfTermsOnLeftContainingVariables == sumOfNonVariablesOnRight:
			#print "Applying J4 to equations " + eqn.PrintGraph()		
			SetLeftTermsToValue(eqn, 1, varValues, zValueDict)
		elif sumOfNonVariablesOnLeft + numberOfTermsOnLeftContainingVariables < sumOfNonVariablesOnRight:
			#print "WE FOUND A CONTRADICTION 2"
			contradiction = 1
		

def Judgement5(varValues, strP, strQ):
	for i in varValues.keys():
		if varValues[i] == 0:
			#Check that there is only a single letter in the key, otherwise it already has form p1q1
			numLetter = ReHandler.CountLettersInString(i)
			if numLetter == 1:
				#Apply J5
				term1 = i
				for i in range(1, len(strP)-1):
					term2 = strP[i]
					string1 = term1 + term2
					string2 = term2 + term1
					SetVarValueAndCheckForContradiction(varValues, string1, 0)
					SetVarValueAndCheckForContradiction(varValues, string2, 0)
				for i in range(1, len(strQ)-1):
					term2 = strQ[i]
					string1 = term1 + term2
					string2 = term2 + term1
					SetVarValueAndCheckForContradiction(varValues, string1, 0)
					SetVarValueAndCheckForContradiction(varValues, string2, 0)

#if x = 1, then xy = y					
def Judgement8(varValues, strP, strQ):
	for j in varValues.keys():
		if varValues[j] == 1:
			numLetter = ReHandler.CountLettersInString(j)
			if numLetter == 1:
				term1 = j
				for i in range(1, len(strP)-1):
					term2 = strP[i]
					string1 = term1 + term2
					string2 = term2 + term1
					if string1 not in varValues:
						varValues[string1] = term2
					if string2 not in varValues:
						varValues[string2] = term2
					
				for i in range(1, len(strQ)-1):
					term2 = strQ[i]
					string1 = term1 + term2
					string2 = term2 + term1
					if string1 not in varValues:
						varValues[string1] = term2
					if string2 not in varValues:
						varValues[string2] = term2
						
def CheckForVarAssignmentContradictions(varValues, strP, strQ):
	global contradiction
	for j in varValues.keys():
		if varValues[j] == 0:
			numLetter = ReHandler.CountLettersInString(j)
			if numLetter == 1:	#Ensures only a single p or q value
				term1 = j
				for i in range(1, len(strP)-1):
					term2 = strP[i]
					string1 = term1 + term2
					string2 = term2 + term1	
					if string1 in varValues:
						if varValues[string1] == 1:
							contradiction = 1
					if string2 in varValues:
						if varValues[string2] == 1:
							contradiction = 1								
				for i in range(1, len(strQ)-1):
					term2 = strQ[i]
					string1 = term1 + term2
					string2 = term2 + term1				
					if string1 in varValues:
						if varValues[string1] == 1:
							contradiction = 1
					if string2 in varValues:
						if varValues[string2] == 1:
							contradiction = 1
						
def Judgement6(varValues):
	listOfVariables = []
	for i in varValues.keys():
		if varValues[i] == 1:
			numLetter = ReHandler.CountLettersInString(i)
			if numLetter == 1:
				listOfVariables.append(i)
	
	for i in listOfVariables:
		for j in listOfVariables:
			string1 = i + j
			string2 = j + i
			SetVarValueAndCheckForContradiction(varValues, string1, 1)
			SetVarValueAndCheckForContradiction(varValues, string2, 1)
				
def Judgement7(eqn, varValues, zValueDict):
	global contradiction
	'A decision tree of cases Ive encountered and figured out, may add things later'
	sumOfNonVariablesOnLeft = eqn.left.CountNonVariables()
	sumOfNonVariablesOnRight = eqn.right.CountNonVariables()

	if sumOfNonVariablesOnLeft > sumOfNonVariablesOnRight and eqn.right.CountTerms() == 1 and eqn.right.CountTermsContainingVariables() == 0:
		#Before recombining, there is only a single term on the right that is overpowered
		#by the term on the left.
		#print "WE FOUND A CONTRADICTION 3"
		contradiction = 1
	else:	
		constant = sumOfNonVariablesOnRight - sumOfNonVariablesOnLeft
		
		#In Case we have NonVariablesOnBothSides, Recombine Equations
		eqn.BalanceNaturalNumbers()
		
		if constant > 0:			
			#Constant on Right
			countTermsOnLeft = eqn.left.CountTerms()
			countTermsOnRight = eqn.right.CountTerms()
			if countTermsOnLeft == constant:
				SetLeftTermsToValue(eqn, 1, varValues, zValueDict)
				SetAllRightTermsToValue(eqn, 0, zValueDict)				
			elif countTermsOnLeft > constant:				
				if countTermsOnRight == 2:					
					variableOnRight = eqn.right.FindLargestMultiplier(-1, 0)
					if variableOnRight+constant > countTermsOnLeft:
						SetLargestRightTermToValue(eqn,variableOnRight, 0, zValueDict)	
						
			elif countTermsOnLeft < constant:
				if countTermsOnRight == 1:
					#print "WE FOUND A CONTRADICTION 4"
					contradiction = 1
		elif constant < 0:
			#Constant on Left
			constant = abs(constant)
			countVariablesOnRight = eqn.right.CountTermsContainingVariables()
			countTermsOnRight = eqn.right.CountTerms()
			if countTermsOnRight == 1 and countVariablesOnRight == 1:
				#single variable term on right
				rightVariable = eqn.right.FindLargestMultiplier(-1,0)
				if rightVariable < constant:
					#print "WE FOUND A CONTRADICTION 5"
					contradiction = 1
				elif rightVariable == constant:
					#set right variable to 1
					SetLargestRightTermToValue(eqn, rightVariable, 1, zValueDict)				
					#set left variable to 0
					SetLeftTermsToValue(eqn, 0, varValues, zValueDict)
				elif rightVariable > constant:
					maxPossibleOnLeft = constant + eqn.left.CountTermsContainingVariables()
					if maxPossibleOnLeft < rightVariable:
						#print "WE FOUND A CONTRADICTION 6"
						contradiction = 1
					elif maxPossibleOnLeft == rightVariable:
						SetLargestRightTermToValue(eqn, rightVariable, 1, zValueDict)
						SetLeftTermsToValue(eqn, 1, varValues, zValueDict)
			elif countTermsOnRight > 1 and countVariablesOnRight > 1 and countTermsOnRight == countVariablesOnRight:
				#multiple variable terms on right
				countTermsOnLeft = eqn.left.CountTerms()
				sumOfRightVars = GetSumOfAllRightVars(eqn)
				if countTermsOnLeft == 1:
					#single term (constant) on left
					if sumOfRightVars == constant:
						SetAllRightTermsToValue(eqn, 1, zValueDict)
					elif sumOfRightVars > constant:
						IfSingleCombinationSummingToParameterSetValuesAppropriately(eqn, constant, zValueDict)
					elif sumOfRightVars < constant:
						#print "WE FOUND A CONTRADICTION 7"
						contradiction = 1
				elif countTermsOnLeft > 1:
					#multiple terms (constant + variable) on left
					if AllRightCombinationsGuaranteedToBeEven(eqn) == 1:
						if constant%2 == 0:
							countVariablesOnLeft = eqn.left.CountTermsContainingVariables()
							if countVariablesOnLeft == 1:
								SetLeftTermsToValue(eqn, 0, varValues, zValueDict)
					if eqn.right.CountTerms() == 2:
						rightTerms = GetAllRightVars(eqn)
						higherValue = rightTerms[1]
						lowerValue = rightTerms[0]
						if rightTerms[0] != rightTerms[1]:
							if rightTerms[0] > rightTerms[1]:
								higherValue = rightTerms[0]
								lowerValue = rightTerms[1]
								constantPlusLeftVar = constant + eqn.left.CountTermsContainingVariables()
								if constant > lowerValue and constantPlusLeftVar < lowerValue+higherValue :
									SetLargestRightTermToValue(eqn, higherValue, 1, zValueDict)
									SetLargestRightTermToValue(eqn, lowerValue, 0, zValueDict)
								elif constant > lowerValue and constantPlusLeftVar == lowerValue+higherValue:
									SetLargestRightTermToValue(eqn, higherValue, 1, zValueDict)
									#Triple check this assumption
								
								
def SubstituteVarsForValues(eqn, zValueDict, varValues, eqnValues):
	for i in zValueDict.keys():
		if zValueDict[i] != -1:
			eqn.ChangeValue(i, zValueDict[i])
	
	for i in varValues.keys():
		eqn.ChangeValue(i, varValues[i])
		
	for i in eqnValues.keys():
		if eqn.ContainsArgumentVariables(i) == 1:
			if len(i) <= eqn.left.CountTerms() and eqn.right.CountTermsContainingVariables() > 0:
				eqn.left.ReplaceArgumentVariablesWithValue(i, eqnValues[i])
	eqn.CleanUpMultiples()

def SubstitutezVarsForValues(eqn, zValueDict):
	for i in zValueDict.keys():
		if zValueDict[i] != -1:
			eqn.ChangeValue(i, zValueDict[i])
	eqn.CleanUpMultiples()

	
def IntermediatePassRules(eqn):
	# print "IntermediatePassRules"
	# #Cancel 1's if on both sides
	# #combine terms on right side (ie a 1 and a 1 will become a 2)
	# #move terms from left to right IF they allow us to take advantage of J3
	if eqn.left.ContainsOne() == 1 and eqn.right.ContainsOne() == 1:
		eqn.RemoveSingleOneEachSide()
		
def SetLeftTermsToValue(eqn, value, varValues, zValueDict):
	leftTerms = eqn.left.GetTerms()
	for i in leftTerms:
		if ReHandler.IsInt(i) < 0:
			#its a variable
			if ReHandler.MatchZVariable(i):
				SetZValueAndCheckForContradiction(zValueDict, i, value)
			else:
				SetVarValueAndCheckForContradiction(varValues, i, value)

def SetZValueAndCheckForContradiction(zValueDict, key, value):
	global contradiction
	if zValueDict[key] != -1 and zValueDict[key] != value:
		#print "WE FOUND A CONTRADICTION 8"
		contradiction = 1
	else:
		zValueDict[key] = value

def SetVarValueAndCheckForContradiction(varValues, key, value):
	global contradiction
	if key in varValues:
		if varValues[key] != -1 and varValues[key] != value:
			contradiction = 1
	else:
		varValues[key] = value
		
def SetLargestRightTermToValue(eqn, largestMultiplier, value, zValueDict):
	x = eqn.right.RetrieveAndDeleteTermWithLargestMultiplier(largestMultiplier)
	for i in x:
		tuple = ReHandler.MatchZVariable(i)
		if tuple:
			SetZValueAndCheckForContradiction(zValueDict, tuple, value)
		else:
			if value != 0:				
				tempnode = Node(i) #checked
				eqn.right.AddChild(tempnode)

def SetAllRightTermsToValue(eqn, value, zValueDict):
	largestMultiplier = eqn.right.FindLargestMultiplier(-1, 0)
	
	while largestMultiplier > 0:
		x = eqn.right.RetrieveAndDeleteTermWithLargestMultiplier(largestMultiplier)
		for i in x:
			tuple = ReHandler.MatchZVariable(i)
			if tuple:
				SetZValueAndCheckForContradiction(zValueDict, tuple, value)
			else:
				if value != 0:
					tempnode = Node(i)
					eqn.right.AddChild(tempnode)
		largestMultiplier = eqn.right.FindLargestMultiplier(largestMultiplier, 0)
		
def GetSumOfAllRightVars(eqn):
	allRightVars = GetAllRightVars(eqn)
	total = 0
	for i in allRightVars:
		total += int(i)
	return total
	
def GetAllRightVars(eqn):
	allRightVars = []
	largestMultiplier = eqn.right.FindLargestMultiplier(-1, 0)
	
	while largestMultiplier > 0:
		allRightVars.append(largestMultiplier)
		largestMultiplier = eqn.right.FindLargestMultiplier(largestMultiplier, 0)
		
	return allRightVars

def AllRightCombinationsGuaranteedToBeEven(eqn):
	allRightVars = GetAllRightVars(eqn)
	for i in allRightVars:
		if (int(i))%2 == 1:
			return 0
	return 1
	
def IfSingleCombinationSummingToParameterSetValuesAppropriately(eqn, parValue, zValueDict):
	correctCombinations = FindCombinationsSummingToParameter(eqn, parValue)
	
	if len(correctCombinations) == 1:
		largestMultiplier = eqn.right.FindLargestMultiplier(-1, 0)
		
		while largestMultiplier > 0:
			x = eqn.right.RetrieveAndDeleteTermWithLargestMultiplier(largestMultiplier)
			thisTermIsInCombination = 0
			for i in x:
				tuple = ReHandler.MatchZVariable(i)
				if not tuple:
					for j in correctCombinations[0]:
						if str(i) == str(j):
							thisTermIsInCombination = 1
							
			if thisTermIsInCombination == 1:
				value = 1
			else:
				value = 0
				
			for i in x:
				tuple = ReHandler.MatchZVariable(i)
				if tuple:
					SetZValueAndCheckForContradiction(zValueDict, tuple, value)
				else:
					if value != 0:					
						tempnode = Node(i) #checked
						eqn.right.AddChild(tempnode)
			largestMultiplier = eqn.right.FindLargestMultiplier(largestMultiplier, 0)

def FindCombinationsSummingToParameter(eqn, parValue):
	correctCombinations = []
	allRightVars = GetAllRightVars(eqn)
	varCombinations = GetAllCombinations(allRightVars)
	
	for i in varCombinations:
		total = 0
		for j in i:
			total += j
		if total == parValue:
			correctCombinations.append(i)
	return correctCombinations
	
def GetAllCombinations(vars):
	varCombinations = []
	for i in range(0, len(vars)):
		x = itertools.combinations(vars, i+1)
		for j in x:
			varCombinations.append(j)
	return varCombinations
	
def CountQubitsRemaining(zValues, varValues, strP, strQ):
	count = 0
	for i in zValues.keys():
		if zValues[i] == -1:
			count += 1
	
	for i in range(1, len(strP)-1):
		if strP[i] not in varValues:
			count += 1
			
	for i in range(1, len(strQ)-1):
		if strQ[i] not in varValues:
			count += 1
	
	return count

def MaximumQubits(zValues, strP, strQ):
	count = 0
	count += len(zValues)
	count += len(strP)-2
	count += len(strQ)-2
	return count
	
def FinalValueCheck(originalEqns, zValues, varValues):
	copyOfOriginalEqns = deepcopy(originalEqns)
	global contradiction

	blankEqnValues = {}
	newVarValues = {}
	
	
	for i in copyOfOriginalEqns:
		#SubstitutezVarsForValues(i, zValues )
		SubstituteVarsForValues(i, zValues, varValues, blankEqnValues )
	for i in copyOfOriginalEqns:
		i.CleanUpMultiples()
		
	for i in copyOfOriginalEqns:
		x = i.left.CountNonVariables() + i.left.CountTermsContainingVariables()
		y = i.right.CountNonVariables()
		if x < y:
			contradiction = 1
		
	if contradiction == 1:
		contradiction = 0
		return 1
	else:
		for i in copyOfOriginalEqns:
			print i.PrintGraph()
	return 0
	