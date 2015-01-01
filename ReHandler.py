#!/usr/bin/env python

import re

def get_coefficient(expr):
    ''' Return the leading coefficient of a string
        >>> get_coefficient('23')
        '23'
        >>> get_coefficient('34z1')
        '34'
        >>> get_coefficient('2q2p91')
        '2'
        >>> get_coefficient('z3')
        '1'
        >>> get_coefficient('1')
        '1'
    '''
    if isinstance(expr, int):
        return str(expr)
    coef = re.match('([0-9]*)(.*)', expr).groups()[0]
    if coef == '':
        return '1'
    return coef

def get_variables(expr):
    ''' Return the variables in a string 
    
    >>> get_variables('z2')
    ['z2']
    >>> get_variables('32')
    []
    >>> get_variables('q23p7')
    ['q23', 'p7']
    >>> get_variables('4qp9')
    ['q', 'p9']
    '''
    coef, variables = re.match('([0-9]*)(.*)', expr).groups(0)
    variables = re.findall('([a-z][0-9]*)', variables)
    return variables

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
    try:
        y = len(re.findall('[a-z]', param))
    except:
        print param
    return y
	
def IsInt(string1):
		try:
			x = int(string1)
			return x
		except ValueError:
			return -1


if __name__ == "__main__":
    import doctest
    doctest.testmod()