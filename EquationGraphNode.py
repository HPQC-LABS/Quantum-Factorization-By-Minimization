#!/usr/bin/env python

import sys
import pdb
import ReHandler

__author__ = "Nathaniel Bryans"
__credits__ = ["Nathaniel Bryans", "Nikesh Dattani"]
__version__ = "0.0.4"
__status__ = "Prototype"

class Node:
	'Node for an equation solving graph'
	
	def __init__(self, value):
		self.value = value
		self.children = []
		self.left = ""		#only graph head has defined left
		self.right = "" 	#only graph head has defined right
		
	def AddChild(self, child):
		self.children.append(child)
		
	def AddLeft(self, child):
		self.left = child
	
	def AddRight(self, child):
		self.right = child
			
	def GetTerms(self):
		values = []
		if self.value == "=":
			return self.left.GetTerms().append(self.right.GetTerms())
		if not self.children:
			return self.value
		elif self.value == "x":
			#Concatenate variable and constant
			for i in range(0, len(self.children)):
				if ReHandler.IsInt(self.children[i].value) > 0:
					number = str(self.children[i].value)
				else:
					variable = str(self.children[i].value)
			totalString = number+variable
			values.append(totalString)
			return values
		elif self.value == "+":
			for i in range(0, len(self.children)):
				values.append(self.children[i].GetTerms())
			return values
	
	def CountTerms(self):
		numTerms = 0
		if self.value == "=": 	#head of graph
			return self.left.CountTerms() + self.right.CountTerms()
		if not self.children:
			return 1
		if self.value == "x":
			return 1
		else:
			for i in range(0, len(self.children)):
				numTerms += self.children[i].CountTerms()
			return numTerms
			
	def CountTermsContainingVariables(self):
		numTerms = 0
		if self.value == "=":
			return self.left.CountTermsContainingVariables() + self.right.CountTermsContainingVariables()
		if not self.children:
			if ReHandler.IsInt(self.value) > 0:
				return 0
			else:
				return 1
		if self.value == "x":
			return 1
		if self.value == "+":
			for i in range(0, len(self.children)):
				numTerms += self.children[i].CountTermsContainingVariables()
			return numTerms

	def CountNonVariables(self):
		total = 0
		if self.value == "=":	#head of list
			return self.left.CountNonVariables() + self.right.CountNonVariables()
		if not self.children:		#leaf
			x = ReHandler.IsInt(self.value)
			if x > 0:
				return x
			else:
				return 0
		else:						#not leaf
			if self.value == "x":
				return 0
			else:
				for i in range(0, len(self.children)):
					total += self.children[i].CountNonVariables()
				return total
	
	def TrimNonVariables(self):
		if self.value == "=":
			self.left.TrimNonVariables()
			self.right.TrimNonVariables()
		if not self.children: #leaf
			return ReHandler.IsInt(self.value)
		if self.value == "x":
			return 0
		if self.value == "+":
			numToTrim = []
			for i in range(0, len(self.children)):
				x = self.children[i].TrimNonVariables()
				if x > 0:
					numToTrim.append(i)
			count = 0
			for i in numToTrim:
				self.children.remove(self.children[i-count])
				count += 1		
		
	def BalanceNaturalNumbers(self):
		if self.value != "=":
			#error
			print "Error, we must start at with the root node"
			exit()
		
		numOnLeft = self.left.CountNonVariables()
		numOnRight = self.right.CountNonVariables()
		self.left.TrimNonVariables()
		self.right.TrimNonVariables()
				
		#Add new node to right
		newNum = numOnRight-numOnLeft
		
		if newNum > 0:
			tempnode = Node(newNum)
			self.right.AddChild(tempnode)
		elif newNum < 0:
			#WHY ADD NEGATIVES WHEN WE CAN JUST ADD A POSITIVE TO THE LEFT SIDE
			tempnode = Node(abs(newNum))
			self.left.AddChild(tempnode)
				
	def PrintGraph(self):
		string = ""
		if self.value == "=":	#head of list
			string = string + self.left.PrintGraph() + " " + str(self.value) + " " + self.right.PrintGraph()
			return string
		else:
			if not self.children:
				return str(self.value)
			else:
				for i in range(0, len(self.children)-1):
					string = string + self.children[i].PrintGraph() + " " + self.value + " "
				string = string + self.children[len(self.children)-1].PrintGraph()
				return string
				
	def FindLargestMultiplier(self, maxNum, maxFound):
		if maxNum == -1:	#we want to find largest
			maxNum = sys.maxint
		if self.value == "=":	#We only want this called on left/right branch, not both at once			
			return 0
		else:
			if not self.children:	#If it's a leaf
				x = ReHandler.IsInt(self.value)
				if x > 0 and x < maxNum and x > maxFound:
					return x
				else:
					return 0
			else:					#If it's NOT a leaf
				if self.value == "x":
					for i in range(0, len(self.children)):
						y = self.children[i].FindLargestMultiplier(maxNum, maxFound)
						if y > maxFound:
							maxFound = y
					return maxFound
				elif self.value == "+":
					for i in range(0, len(self.children)):
						if self.children[i].value == "x":
							y = self.children[i].FindLargestMultiplier(maxNum, maxFound)
							if y > maxFound:
								maxFound = y
					return maxFound					
	
	def RetrieveAndDeleteTermWithLargestMultiplier(self, multiplier):
		if self.value == "=":			
			return ""
		else:
			foundMultiplier = 0
			leafList = []
			if not self.children:	#If it's a leaf	
				if self.value == str(multiplier):
					return 1
				else:
					return 0
			else:					#If it's NOT a leaf
				if self.value == "x":
					for i in range(0, len(self.children)):
						if self.children[i].RetrieveAndDeleteTermWithLargestMultiplier(multiplier) != 0:
							#This is an "x" and we found the multiplier attached to it
							#return z value and delete this branch
							foundMultiplier = 1
					if foundMultiplier > 0:
						for i in range(0, len(self.children)):
							leafList.append(self.children[i].value)
					return leafList
				elif self.value == "+":
					for i in range (0, len(self.children)):
						if self.children[i].value == "x":
							x = self.children[i].RetrieveAndDeleteTermWithLargestMultiplier(multiplier)
							if x:
								self.children.remove(self.children[i])
								return x
					return leafList
		return leafList
	
	def ChangeValue(self, valueOld, valueNew):
		'searches for the oldValue (should be a variable) and replaces if with the new one (should be constant)'
		if self.value == "=":
			self.left.ChangeValue(valueOld, valueNew)
			self.right.ChangeValue(valueOld, valueNew)
		if not self.children:
			return
		if self.value == "+":
			numToRemove = []
			for i in range(0, len(self.children)):				
				if self.children[i].value == "x":		
					numToRemoveFromChild = []
					for j in range(0, len(self.children[i].children)):				
						if self.children[i].children[j].value == valueOld:
							if str(valueNew) == "0":
								j = len(self.children[i].children)
								numToRemove.append(i) #take out the whole x branch because of '0'
							elif str(valueNew) == "1":
								numToRemoveFromChild.append(j)								
							else:								
								self.children[i].children[j].value = valueNew #checked
					for x in numToRemoveFromChild:
						self.children[i].children.remove(self.children[i].children[x])
				elif self.children[i].value == valueOld:
					if str(valueNew) != "0":
						self.children[i].value = valueNew #checked
					else:
						numToRemove.append(i)
			for x in numToRemove:
				self.children.remove(self.children[x])
				if len(self.children) == 0:
					newZeroNode = Node("0")
					self.AddChild(newZeroNode)
	
	def CleanUpMultiples(self):
		'removes Xs from tree where necessary'
		if self.value == "=":
			self.left.CleanUpMultiples()
			self.right.CleanUpMultiples()
		if not self.children:
			return
		if self.value == "+":
			for i in range(0, len(self.children)):
				if self.children[i].value == "x":
					if len(self.children[i].children) == 1:
						tempNode = Node(self.children[i].children[0].value)
						self.AddChild(tempNode)
						self.children.remove(self.children[i])
		
	def ContainsArgumentVariables(self, varTuple):
		'given a set of variables, returns true if the left side contains all those variables'
		if self.value == "=":
			return self.left.ContainsArgumentVariables(varTuple)
		if self.value == "+":
			for i in varTuple:
				found = 0
				for j in range(0, len(self.children)):
					if str(self.children[j].value) == i:
						found = 1
				if found != 1:
					return 0
			return 1
		return 0
	
	def ReplaceArgumentVariablesWithValue(self, varTuple, value):
		'given, for example, p1 + q1 = 1, replaces p1 + q1 with 1'
		#print "******************** Replacing Left Side variables********************"
		#print varTuple
		if self.value == "=":
			return self.left.ReplaceLeftSideVariablesWithValue(varTuple, value)
		if self.value == "+":
			for i in varTuple:
				numToRemove = -1
				for j in range(0, len(self.children)):
					if str(self.children[j].value) == str(i):
						numToRemove = j
				self.children.remove(self.children[numToRemove])
			newNode = Node(value)
			self.AddChild(newNode)
			
	def RemoveZeroes(self):
		if self.value == "=":
			self.left.RemoveZeroes()
			self.right.RemoveZeroes()
		if not self.children:
			return
		if self.value == "+":
			numToRemove = []
			for i in range(0, len(self.children)):
				if self.children[i].value == "0":
					numToRemove.append(i)
			for i in numToRemove:
				self.children.remove(self.children[i])
				
	def ContainsZero(self):
		if self.value == "=":
			if self.left.ContainsZero() == 1:
				return 1
			return self.right.ContainsZero()
		if self.value == "+":
			for i in range(0, len(self.children)):
				if self.children[i].value == "0":
					return 1
			return 0
	
	def ContainsLoneZeroOnRight(self):
		if self.value == "=":
			return self.right.ContainsLoneZeroOnRight()
		if self.value == "+":
			if len(self.children) == 0:
				return 1
			elif len(self.children) == 1:
				if str(self.children[0].value) == "0":
					return 1
				else:
					return 0
	
	def ContainsOne(self):
		if self.value == "=":
			return self.left.ContainsOne() + self.right.ContainsOne()
		if self.value == "+":
			for i in range(0, len(self.children)):
				if str(self.children[i].value) == "1":
					return 1
			return 0
	
	def RemoveSingleOneEachSide(self):
		indexToRemove = -1
		if self.value == "=":
			self.left.RemoveSingleOneEachSide()
			self.right.RemoveSingleOneEachSide()
		if self.value == "+":
			for i in range(0, len(self.children)):
				if str(self.children[i].value) == "1":
					indexToRemove = i
			self.children.remove(self.children[indexToRemove])
			if len(self.children) == 0:
				tempZeroNode = Node("0")
				self.AddChild(tempZeroNode)