#!/usr/bin/env python

from EquationGraphNode import Node
import EquationSolver

def TestModifiedJ1(varValue, eqnValue, J1sApplied):
	print "TEST ModJ1"
	node1 = Node("p10")
	node2 = Node("q1")
	node7 = Node("p2")
	node3 = Node("+")
	node3.AddChild(node1)
	node3.AddChild(node2)
	node3.AddChild(node7)
	node4 = Node("+")
	node5 = Node("1")
	node4.AddChild(node5)
	node6 = Node("=")
	node6.AddLeft(node3)
	node6.AddRight(node4)
	print node6.PrintGraph()
	print varValue
	print eqnValue
	EquationSolver.ModifiedJudgement1(node6, varValue, eqnValue, J1sApplied)
	print varValue
	print eqnValue

def TestModifiedJ1Test2(varValue,eqnValue, J1sApplied):
	print "Test ModJ1 2"
	node1 = Node("p1")
	node3 = Node("+")
	node3.AddChild(node1)
	node4 = Node("+")
	node5 = Node("1")
	node4.AddChild(node5)
	node6 = Node("=")
	node6.AddLeft(node3)
	node6.AddRight(node4)
	print node6.PrintGraph()
	print varValue
	EquationSolver.ModifiedJudgement1(node6, varValue, eqnValue, J1sApplied)
	print varValue

def TestModifiedJ1Test3(varValue,eqnValue, J1sApplied):
	print "Test ModJ1 3"
	node1 = Node("p12")
	node7 = Node("1")
	node3 = Node("+")
	node3.AddChild(node1)
	node3.AddChild(node7)
	node4 = Node("+")
	node5 = Node("1")
	node4.AddChild(node5)
	node6 = Node("=")
	node6.AddLeft(node3)
	node6.AddRight(node4)
	print node6.PrintGraph()
	print varValue
	EquationSolver.ModifiedJudgement1(node6, varValue, eqnValue, J1sApplied)
	print varValue	
	
def TestModifiedJ1Test4(varValue, eqnValue, J1sApplied):
	print "Test ModJ1 4"
	node1 = Node("p12")
	node2 = Node("q2")
	node3 = Node("+")
	node4 = Node("1")
	node5 = Node("+")
	node6 = Node("=")
	
	node3.AddChild(node1)
	node3.AddChild(node2)
	node5.AddChild(node4)
	node6.AddLeft(node3)
	node6.AddRight(node5)
	
	print node6.PrintGraph()
	print varValue
	print eqnValue
	EquationSolver.ModifiedJudgement1(node6, varValue, eqnValue, J1sApplied)
	print node6.PrintGraph()
	print varValue
	print eqnValue
	
def TestJ2(zValueDict):
	print "Test J2"
	node1 = Node("p12")
	node2 = Node("p2")
	node3 = Node("p3")
	node4 = Node("p14")
	node5 = Node("+")
	node5.AddChild(node1)
	node5.AddChild(node2)
	node5.AddChild(node3)
	node5.AddChild(node4)
	node6 = Node("4")
	node7 = Node("z12")
	node8 = Node("x")
	node8.AddChild(node6)
	node8.AddChild(node7)
	node9 = Node("1")
	node10 = Node("+")
	node10.AddChild(node8)
	node10.AddChild(node9)
	node11 = Node("=")
	node11.AddLeft(node5)
	node11.AddRight(node10)
	print node11.PrintGraph()
	print zValueDict
	EquationSolver.Judgement2(node11, zValueDict)
	print node11.PrintGraph()
	print zValueDict

def TestJ2Test2(zValueDict):
	print "Test J2"
	node1 = Node("p12")
	node2 = Node("p2")
	node3 = Node("p3")
	node4 = Node("p4")
	node12 = Node("p5")
	node5 = Node("+")
	node5.AddChild(node1)
	node5.AddChild(node2)
	node5.AddChild(node3)
	node5.AddChild(node4)
	node5.AddChild(node12)
	node6 = Node("4")
	node7 = Node("z12")
	node8 = Node("x")
	node8.AddChild(node6)
	node8.AddChild(node7)
	node9 = Node("1")
	node10 = Node("+")
	node10.AddChild(node8)
	node10.AddChild(node9)
	node11 = Node("=")
	node11.AddLeft(node5)
	node11.AddRight(node10)
	print node11.PrintGraph()
	print zValueDict
	EquationSolver.Judgement2(node11, zValueDict)
	print node11.PrintGraph()
	print zValueDict

def TestJ3(zValueDict):
	print "TEST J3"
	node1 = Node("p12")
	node2 = Node("q1")
	node4 = Node("1")
	node3 = Node("+")
	node3.AddChild(node1)
	node3.AddChild(node2)
	node3.AddChild(node4)
	node5 = Node("z12")
	node7 = Node("2")
	node6 = Node("x")
	node8 = Node("+")
	node6.AddChild(node7)
	node6.AddChild(node5)
	node8.AddChild(node6)
	nodeparent = Node("=")
	nodeparent.AddLeft(node3)
	nodeparent.AddRight(node8)
	print nodeparent.PrintGraph()
	print zValueDict
	EquationSolver.Judgement3(nodeparent, zValueDict)
	print nodeparent.PrintGraph()
	print zValueDict

def TestJ4(varValues, zValueDict):
	node1 = Node("p11")
	node2 = Node("z12")
	node3 = Node("+")
	node3.AddChild(node1)
	node3.AddChild(node2)
	node4 = Node("2")
	node5 = Node("+")
	node5.AddChild(node4)
	node6 = Node("=")
	node6.AddLeft(node3)
	node6.AddRight(node5)
	print node6.PrintGraph()
	print varValues
	print zValueDict
	EquationSolver.Judgement4(node6, varValues, zValueDict)
	print node6.PrintGraph()
	print varValues
	print zValueDict

def TestJ3(varValues, zValueDict):
	node1 = Node("p1")
	node2 = Node("z12")
	node3 = Node("+")
	node3.AddChild(node1)
	node3.AddChild(node2)
	node4 = Node("5")
	node5 = Node("+")
	node5.AddChild(node4)
	node6 = Node("=")
	node6.AddLeft(node3)
	node6.AddRight(node5)
	print node6.PrintGraph()
	print varValues
	print zValueDict
	EquationSolver.Judgement4(node6, varValues, zValueDict)
	print node6.PrintGraph()
	print varValues
	print zValueDict
	
def TestJ4Test3(varValues, zValueDict):
	node1 = Node("p12")
	node2 = Node("z12")
	node7 = Node("3")
	node3 = Node("+")
	node3.AddChild(node1)
	node3.AddChild(node2)
	node3.AddChild(node7)
	node4 = Node("5")
	node5 = Node("+")
	node5.AddChild(node4)
	node6 = Node("=")
	node6.AddLeft(node3)
	node6.AddRight(node5)
	print node6.PrintGraph()
	print varValues
	print zValueDict
	EquationSolver.Judgement4(node6, varValues, zValueDict)
	print node6.PrintGraph()
	print varValues
	print zValueDict
	
def TestJ5(varValues,strP, strQ):
	print "TestJ5"
	string = "p12"
	varValues[string] = 0
	print varValues
	EquationSolver.Judgement5(varValues, strP, strQ)
	print varValues
	
def TestJ6(varValues):
	string1 = "p1"
	string2 = "q12"
	string3 = "p2"
	string4 = "q2"
	varValues[string1] = 1
	varValues[string2] = 1
	varValues[string3] = 0
	varValues[string4] = 1
	print varValues
	EquationSolver.Judgement6(varValues)
	print varValues
	EquationSolver.Judgement6(varValues)
	print varValues
	
def TestJ7(varValues, zValues):
	'4 = 2z12 + 4z23'
	node1 = Node("4")
	node2 = Node("+")
	node3 = Node("2")
	node4 = Node("z12")
	node5 = Node("x")
	node6 = Node("4")
	node7 = Node("z23")
	node8 = Node("x")
	node9 = Node("+")
	node10 = Node("=")
	
	node2.AddChild(node1)
	node5.AddChild(node3)
	node5.AddChild(node4)
	node8.AddChild(node6)
	node8.AddChild(node7)
	node9.AddChild(node5)
	node9.AddChild(node8)
	node10.AddLeft(node2)
	node10.AddRight(node9)
	
	print node10.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node10, varValues, zValues)
	print node10.PrintGraph()
	print varValues
	print zValues
	
def TestJ7Test2(varValues, zValues):
	'6 = 2z12 + 4z23'
	node1 = Node("6")
	node2 = Node("+")
	node3 = Node("2")
	node4 = Node("z12")
	node5 = Node("x")
	node6 = Node("4")
	node7 = Node("z23")
	node8 = Node("x")
	node9 = Node("+")
	node10 = Node("=")
	
	node2.AddChild(node1)
	node5.AddChild(node3)
	node5.AddChild(node4)
	node8.AddChild(node6)
	node8.AddChild(node7)
	node9.AddChild(node5)
	node9.AddChild(node8)
	node10.AddLeft(node2)
	node10.AddRight(node9)
	
	print node10.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node10, varValues, zValues)
	print node10.PrintGraph()
	print varValues
	print zValues
	
def TestJ7Test3(varValues, zValues):
	'7 = 2z12 + 4z23'
	node1 = Node("7")
	node2 = Node("+")
	node3 = Node("2")
	node4 = Node("z12")
	node5 = Node("x")
	node6 = Node("4")
	node7 = Node("z23")
	node8 = Node("x")
	node9 = Node("+")
	node10 = Node("=")
	
	node2.AddChild(node1)
	node5.AddChild(node3)
	node5.AddChild(node4)
	node8.AddChild(node6)
	node8.AddChild(node7)
	node9.AddChild(node5)
	node9.AddChild(node8)
	node10.AddLeft(node2)
	node10.AddRight(node9)
	
	print node10.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node10, varValues, zValues)
	print node10.PrintGraph()
	print varValues
	print zValues
	
def TestJ7Test4(varValues, zValues):
	'4 + p1 = 2z12 + 4z23'
	node1 = Node("4")
	node11 = Node("p12")
	node2 = Node("+")
	node3 = Node("2")
	node4 = Node("z12")
	node5 = Node("x")
	node6 = Node("4")
	node7 = Node("z23")
	node8 = Node("x")
	node9 = Node("+")
	node10 = Node("=")
	
	node2.AddChild(node1)
	node2.AddChild(node11)
	node5.AddChild(node3)
	node5.AddChild(node4)
	node8.AddChild(node6)
	node8.AddChild(node7)
	node9.AddChild(node5)
	node9.AddChild(node8)
	node10.AddLeft(node2)
	node10.AddRight(node9)
	
	print node10.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node10, varValues, zValues)
	print node10.PrintGraph()
	print varValues
	print zValues
	
def TestJ7Test5(varValues, zValues):
	'3 + p1 + p2 = 2z12 + 4z23'
	node1 = Node("3")
	node11 = Node("p11")
	node12 = Node("p2")
	node2 = Node("+")
	node3 = Node("2")
	node4 = Node("z12")
	node5 = Node("x")
	node6 = Node("4")
	node7 = Node("z23")
	node8 = Node("x")
	node9 = Node("+")
	node10 = Node("=")
	
	node2.AddChild(node1)
	node2.AddChild(node11)
	node2.AddChild(node12)
	node5.AddChild(node3)
	node5.AddChild(node4)
	node8.AddChild(node6)
	node8.AddChild(node7)
	node9.AddChild(node5)
	node9.AddChild(node8)
	node10.AddLeft(node2)
	node10.AddRight(node9)
	
	print node10.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node10, varValues, zValues)
	print node10.PrintGraph()
	print varValues
	print zValues
	
def TestJ7Test6(varValues, zValues):
	'3 + p1 + p2 + p3 = 2z12 + 4z23'
	node1 = Node("3")
	node11 = Node("p1")
	node12 = Node("p2")
	node13 = Node("p13")
	node2 = Node("+")
	node3 = Node("2")
	node4 = Node("z12")
	node5 = Node("x")
	node6 = Node("4")
	node7 = Node("z23")
	node8 = Node("x")
	node9 = Node("+")
	node10 = Node("=")
	
	node2.AddChild(node1)
	node2.AddChild(node11)
	node2.AddChild(node12)
	node2.AddChild(node13)
	node5.AddChild(node3)
	node5.AddChild(node4)
	node8.AddChild(node6)
	node8.AddChild(node7)
	node9.AddChild(node5)
	node9.AddChild(node8)
	node10.AddLeft(node2)
	node10.AddRight(node9)
	
	print node10.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node10, varValues, zValues)
	print node10.PrintGraph()
	print varValues
	print zValues
	
def TestJ7Test7(varValues, zValues):
	'p1 + p2 + z12 = 3 + 4z23'
	node1 = Node("p11")
	node2 = Node("p2")
	node3 = Node("z12")
	node4 = Node("+")
	node5 = Node("3")
	node6 = Node("4")
	node7 = Node("z23")
	node8 = Node("x")
	node9 = Node("+")
	node10 = Node("=")
	
	node4.AddChild(node1)
	node4.AddChild(node2)
	node4.AddChild(node3)
	node8.AddChild(node6)
	node8.AddChild(node7)
	node9.AddChild(node8)
	node9.AddChild(node5)
	node10.AddLeft(node4)
	node10.AddRight(node9)
	
	print node10.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node10, varValues, zValues)
	print node10.PrintGraph()
	print varValues
	print zValues
	
def TestJ7Test8(varValues, zValues):
	'p1 + p2 = 1 + 2z23'
	node1 = Node("p1")
	node2 = Node("p12")
	node4 = Node("+")
	node5 = Node("1")
	node6 = Node("2")
	node7 = Node("z23")
	node8 = Node("x")
	node9 = Node("+")
	node10 = Node("=")
	
	node4.AddChild(node1)
	node4.AddChild(node2)
	node8.AddChild(node6)
	node8.AddChild(node7)
	node9.AddChild(node8)
	node9.AddChild(node5)
	node10.AddLeft(node4)
	node10.AddRight(node9)
	
	print node10.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node10, varValues, zValues)
	print node10.PrintGraph()
	print varValues
	print zValues

def TestJ7Test9(varValues, zValues):
	'3 + p10 = 2z12'
	node1 = Node("3")
	node2 = Node("p1")
	node3 = Node("+")
	node4 = Node("2")
	node5 = Node("z12")
	node6 = Node("x")
	node7 = Node("+")
	node8 = Node("=")
	
	node3.AddChild(node1)
	node3.AddChild(node2)
	node6.AddChild(node4)
	node6.AddChild(node5)
	node7.AddChild(node6)
	node8.AddLeft(node3)
	node8.AddRight(node7)
	
	print node8.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node8, varValues, zValues)
	print node8.PrintGraph()
	print varValues
	print zValues
	
def TestJ7Test10(varValues, zValues):
	'2 + p10 = 2z12'
	node1 = Node("2")
	node2 = Node("p10")
	node3 = Node("+")
	node4 = Node("2")
	node5 = Node("z12")
	node6 = Node("x")
	node7 = Node("+")
	node8 = Node("=")
	
	node3.AddChild(node1)
	node3.AddChild(node2)
	node6.AddChild(node4)
	node6.AddChild(node5)
	node7.AddChild(node6)
	node8.AddLeft(node3)
	node8.AddRight(node7)
	
	print node8.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node8, varValues, zValues)
	print node8.PrintGraph()
	print varValues
	print zValues
	
def TestJ7Test11(varValues, zValues):
	'3 + p10 = 4z12'
	node1 = Node("3")
	node2 = Node("p10")
	node3 = Node("+")
	node4 = Node("4")
	node5 = Node("z12")
	node6 = Node("x")
	node7 = Node("+")
	node8 = Node("=")
	
	node3.AddChild(node1)
	node3.AddChild(node2)
	node6.AddChild(node4)
	node6.AddChild(node5)
	node7.AddChild(node6)
	node8.AddLeft(node3)
	node8.AddRight(node7)
	
	print node8.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node8, varValues, zValues)
	print node8.PrintGraph()
	print varValues
	print zValues
	
def TestJ7Test12(varValues, zValues):
	'2 + p1 = 4z12' 
	print "Test J7-12. We expect a contradiction"
	node1 = Node("2")
	node2 = Node("p1")
	node3 = Node("+")
	node4 = Node("4")
	node5 = Node("z12")
	node6 = Node("x")
	node7 = Node("+")
	node8 = Node("=")
	
	node3.AddChild(node1)
	node3.AddChild(node2)
	node6.AddChild(node4)
	node6.AddChild(node5)
	node7.AddChild(node6)
	node8.AddLeft(node3)
	node8.AddRight(node7)
	
	print node8.PrintGraph()
	print varValues
	print zValues
	EquationSolver.Judgement7(node8, varValues, zValues)
	print node8.PrintGraph()
	print varValues
	print zValues
	

	