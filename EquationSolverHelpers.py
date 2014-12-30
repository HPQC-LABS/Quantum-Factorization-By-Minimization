#!/usr/bin/env python

import sys
import pdb
from copy import deepcopy
from copy import copy

__author__ = "Nathaniel Bryans"
__credits__ = ["Nathaniel Bryans", "Nikesh Dattani"]
__version__ = "0.0.91"
__status__ = "Prototype"

class EquationSolverState:
	'Node for an equation solving graph'
	
	def __init__(self, columnEqns, zValues, varValues, eqnValues, 
		J1sApplied):
		self.columnEqns = deepcopy(columnEqns)
		self.zValues =  copy(zValues)
		self.varValues = copy(varValues)
		self.eqnValues = copy(eqnValues)
		self.J1sApplied = copy(J1sApplied)
		
	def GetState(self, columnEqns, zValuest, varValues, eqnValues, 
		J1sApplied):		
		columnEqns = self.columnEqns
		zValuest = self.zValues
		varValues = self.varValues
		eqnValues = self.eqnValues
		J1sApplied = self.J1sApplied
		
class VarValueAssignment:
	def __init__(self, varString, value):
		self.varString = varString
		self.value = value
		
	def PrintEntry(self):
		print str(self.varString) + ": " + str(self.value)
		
class VarValueAssignmentPicker:
	def __init__ (self, zValues, strP, strQ, varValues, numVarToSet):
		self.VarValueAssignments = []
		self.done = 0
		for i in range(1, len(strP)-1):
			if strP[i] not in varValues:
				x = VarValueAssignment(strP[i], -1)
				self.VarValueAssignments.append(x)
		for i in range(1, len(strQ)-1):
			if strQ[i] not in varValues:
				x = VarValueAssignment(strQ[i], -1)
				self.VarValueAssignments.append(x)
		for i in zValues:
			if zValues[i] == -1:
				x = VarValueAssignment(i, -1)
				self.VarValueAssignments.append(x)
		
		print "NUM VAR= " + str(len(self.VarValueAssignments))
		numVarToSet = min(numVarToSet, len(self.VarValueAssignments))
		
		for i in range(0, numVarToSet):
			self.VarValueAssignments[i].value = 0
	
	def MoreVarValuesToTry(self):
		if self.done == 1:
			return 1
		return 0
		

	def recurAlg(self, dict, depth, maxDepth):
		y = 0
		if depth > 1:
			y = self.recurAlg(dict, depth-1, maxDepth)
		
		if y == 2:
			return 2	#complete
		
		x = self.findXSetIndex(dict, maxDepth-depth)
		
		if y != 1:	
			if dict[x].value == 0:
				dict[x].value = 1
				if depth != 1:	#if not the furthest depth
					self.resetZeros(dict, depth, maxDepth)
				return 1
			else:
				if x < len(dict)-1:
					
					if dict[x+1].value == -1:
						dict[x].value = -1
						dict[x+1].value = 0
						if depth == 1:
							return 1
						else:
							self.resetZeros(dict, depth, maxDepth)
							return 1
					else:
						if depth == maxDepth:
							return 2
				else:
					return 0
		if y == 1:
			return 1
		
	def GetNextVarValue(self, zValues, varValues, depth):
		x = self.recurAlg(self.VarValueAssignments, depth, depth)
		
		if x == 2:
			self.done = 1
		#print "----------------------------------------"
		for i in range(0, len(self.VarValueAssignments)):
			if self.VarValueAssignments[i].value != -1:
				#if the first character is 'z'
				if self.VarValueAssignments[i].varString[0] == 'z':
					zValues[self.VarValueAssignments[i].varString] = self.VarValueAssignments[i].value
				else:
					varValues[self.VarValueAssignments[i].varString] = self.VarValueAssignments[i].value
				#print self.VarValueAssignments[i].varString + ": " + str(self.VarValueAssignments[i].value)
		
				
	def resetZeros(self, dict, depth, maxDepth):
		x = self.findXSetIndex(dict, maxDepth-depth)
		if x != -1:
			for i in range(x+1, len(dict)):
				dict[i].value = -1
			for i in range(0, depth-1):
				dict[x+1+i].value = 0
	
	def findXSetIndex(self, dict, x):
		count = 0
		for i in range(0, len(dict)):
			if dict[i].value != -1:
				if count == x:
					return i
				count += 1
		return -1
