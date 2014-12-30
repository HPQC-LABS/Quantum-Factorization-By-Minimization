#!/usr/bin/env python

__author__ = "Nathaniel Bryans"

def CarryGenerator(params):
	topDigits = params[0]
	bottomDigits = params[1]

	print topDigits

	li = []
	carry = []
	product = []
	maxBinMulStrLen = int(bottomDigits) + int(topDigits) - 1

	for i in range(0,int(bottomDigits)):
		myParams = [int(topDigits), i, maxBinMulStrLen]
		li.append(BuildBinaryMultiplicationString(myParams))
		
	for i in range (maxBinMulStrLen, 0, -1):
		countOnes = 0
		for j in li:
			if j[i-1:i] == "1":
				countOnes += 1
		#add previous carry to countOnes
		if carry:
			countOnes += int(carry[-1])
			
		binRep = bin(countOnes)[2:] #strip the "0b" from the start
		product.append(binRep[-1])
		
		x = bin(countOnes)[2:-1]
		if x:
			z = int(x,2)
			carry.append(int(x,2))#
		else:
			carry.append(0);
			
	#make sure final carry is included
	product.append(carry[-1])

	pairs = GetCarryOrderedPairs([carry])
	
	carryGroupings = []
	maxRowDifference = 0
	
	for i in pairs:
		difference = i[1]-i[0]
		if  difference > maxRowDifference:
			for j in range(0, difference-maxRowDifference):
				carryGroupings.append([])
			maxRowDifference = difference
		carryGroupings[difference-1].append(i)	
	
	return carryGroupings	
	
def GetCarryOrderedPairs(params):
	#input is the carry list
	#output is a list of all carry pairs (x,y) where x
	#	is the source and y is the destination column
	Carrys = []
	curCol = 0
	for i in params[0]:
		binRep = bin(i)[2:]
		if binRep != "0":
			for j in range(0, len(binRep)):
				#Carrys[curCol] = j+1+curCol
				item = curCol,j+1+curCol
				Carrys.append(item)
		curCol += 1
	return Carrys		

def BuildBinaryMultiplicationString(params):
	str1 = "0" * (params[2] - params[0] - params[1])
	str2 = "1" * params[0]
	str3 = "0" * params[1]
	str = str1+str2+str3
	return str