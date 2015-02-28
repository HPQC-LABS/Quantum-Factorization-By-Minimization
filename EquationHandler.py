#!/usr/bin/env python

__author__ = "Nathaniel Bryans"

def GenerateEquations(params):
	numCol = params[0][0]
	colWidth = params[0][1]
	intermediaryOutput = params[0][3]
	answerOutput = params[0][5]
	leftsides = []
	rightsides = []
	
	carry = []
	for i in params[1]:
		for j in i:
			carry.append(j)
			
	for i in reversed(range(0, numCol)):
		#To the left side add entries in the intermediary and carry FROM
		#To the right side add entries in the product and carry TO (time 2^(to-from))
		leftside = []
		rightside = []
		#Read Intermediaries
		for j in intermediaryOutput:
			leftside.append(j[i*colWidth:i*colWidth+colWidth].strip())
			
		#Read Product
		rightside.append(answerOutput[0][i*colWidth:i*colWidth+colWidth].strip())
		
		#Read Carries FROM
		for j in carry:
			if j[1] == numCol-i-1:
				leftside.append("z" + str(j[0]) + str(j[1]))
			if j[0] == numCol-i-1:
				diff = int(j[1])-int(j[0])
				#print diff
				multiplier = 2**diff
				rightside.append(str(multiplier) + "z" + str(j[0]) + str(j[1]))
		
		leftsides.append(leftside)
		rightsides.append(rightside)
		
		#Format Equations
		leftStr = ""
		rightStr = ""
		for i in leftside:
			if i:
				leftStr = leftStr + i + " + " 
		leftStr = leftStr[:-3]
		for i in rightside:
			if i:
				rightStr = rightStr + i + " + "
		rightStr = rightStr[:-3]
#		print leftStr, " = ", rightStr
		
	return [leftsides, rightsides]