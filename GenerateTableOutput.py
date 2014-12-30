#!/usr/bin/env python

__author__ = "Nathaniel Bryans"

def FormatOutput(params):
	#params[0] contains multiplicands
	#params[1] contains multiplier intermediaries
	#params[2] contains carrys
	
	multiplicandsOutput = []
	intermediaryOutput = []
	carryOutput = []
	answerOutput = []
	numCol = len(params[0][0]) + len(params[0][1])
	colWidth = FindLargestMultiplierLength(params[1]) + 1
	
	#Generate Multiplicand strings
	for i in range(0, len(params[0])):
		multiplicandsOutput.append(" "*colWidth*numCol)
		count = 0
		for j in params[0][i]:
			strToInsert = " "*(colWidth-len(j)) + j
			colToPlace = numCol-count-1
			startPlace = colToPlace*colWidth
			#multiplicandsOutput[i][startPlace:startPlace+5] = strToInsert[0:5]
			multiplicandsOutput[i] = multiplicandsOutput[i][:startPlace] + strToInsert + multiplicandsOutput[i][startPlace+colWidth:]
			count += 1

	#Generate Intermediary Multiplier strings
	outsideCount = 0
	for i in range(0, len(params[1])):
		intermediaryOutput.append(" "*colWidth*numCol)
		count = 0
		for j in params[1][i]:
			strToInsert = " "*(colWidth-len(j)) + j
			colToPlace = numCol-count-1-outsideCount
			startPlace = colToPlace*colWidth
			intermediaryOutput[i] = intermediaryOutput[i][:startPlace] + strToInsert + intermediaryOutput[i][startPlace+colWidth:]
			count+= 1
		outsideCount += 1
		
	#Generate Carry strings
	#	Carry's are represented as (x,y) pairs. where x is the col the carry bit is from and y is the col it's going to
	#	Therefore, each string goes in the y col
	for i in range(0, len(params[2])):
		carryOutput.append(" "*colWidth*numCol)
		for j in params[2][i]:
			zStr = "z" + str(j[0]) + str(j[1])
			strToInsert = " "*(colWidth-len(zStr)) + zStr
			colToPlace = numCol-j[1]-1
			startPlace=colToPlace*colWidth
			carryOutput[i] = carryOutput[i][:startPlace] + strToInsert + carryOutput[i][startPlace+colWidth:]	
		
	#Generate/Format prime "answer"
	answerOutput.append(" "*colWidth*numCol)
	for i in range (0, len(params[3])):
		strToInsert = " "*(colWidth-1) + params[3][i]
		colToPlace = i
		startPlace = colToPlace*colWidth
		answerOutput[0] = answerOutput[0][:startPlace] + strToInsert + answerOutput[0][startPlace+colWidth:]
	
	for i in multiplicandsOutput:
		print i
	print "-"*colWidth*numCol
	for i in intermediaryOutput:
		print i
	print "-"*colWidth*numCol
	for i in carryOutput:
		print i
	print "-"*colWidth*numCol
	print answerOutput[0]
	
	output = []
	output.append(numCol)
	output.append(colWidth)
	output.append(multiplicandsOutput)
	output.append(intermediaryOutput)
	output.append(carryOutput)
	output.append(answerOutput)
	
	return output	#useful later for parsing the equation

def FindLargestMultiplierLength(params):
	largest = -1
	for i in params:
		for j in i:
			if len(str(j)) > largest:
				largest = len(str(j))
	return largest	