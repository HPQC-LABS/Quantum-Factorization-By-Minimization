#!/usr/bin/env python

import re

def MatchZVariable(param):
		pattern = re.compile('^z+')
		m = pattern.match(str(param))
		if m:
			return param
		return ""
		
def MatchProductWithNatural(params):
	pattern = re.compile('^\d+\w+')
	m = pattern.match(params)
	if m:
		pattern2 = re.compile('^\d+')
		m2 = pattern2.match(params)
		if m2:
			return m2.end()
	return ""
	
def CountLettersInString(param):
	y = len(re.findall('[a-z]', param))
	return y
	
def IsInt(string1):
		try:
			x = int(string1)
			return x
		except ValueError:
			return -1