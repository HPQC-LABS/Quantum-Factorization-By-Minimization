"""
Created on Fri Dec 26 12:35:16 2014

Solve a system of equations with binary variables

@author: Richard Tanburn
"""
from copy import deepcopy
from collections import defaultdict
import inspect
import itertools
import sympy

import ReHandler

__author__ = "Richard Tanburn"
__credits__ = ["Richard Tanburn", "Nathaniel Bryans", "Nikesh Dattani"]
__version__ = "0.0.1"
__status__ = "Prototype"


class ContradictionException(Exception):
    ''' Raised when conflicting values are found '''
    pass

# Parent equation

class EquationSolver(object):
    ''' Solver of equations '''

    def __init__(self, equations=None, variables=None):
        if equations is None:
            equations = []
        if variables is None:
            variables = {}
        # Dict of string tag to variable instance
        self.variables = variables
        # List of Equation instances to be solved
        self.equations = equations
        # Set of deductions we have made
        self.deductions = {}

        # And keep a nested dictionary of who made them
        self.deduction_record = defaultdict(lambda : defaultdict(list))


        # Final solutions
        self.final_equations = None
        # Final variables
        self.final_variables = None
        # Final solutions
        self.solutions = None

    def copy(self):
        ''' Return a new instance of itself '''
        copy = EquationSolver(deepcopy(self.equations), deepcopy(self.variables))
        copy.deductions = deepcopy(self.deductions)
        return copy

    @classmethod
    def from_params(cls, params):
        ''' Create an instance from the outpu of whatever came before
            params[0][i] contains all leftsides
            params[1][i] contains all rightsides
        '''
        equations, variables = parse_equations(params)
        return cls(equations, variables)

    @property
    def deductions_as_equations(self):
        ''' Return deductions as a list of equations '''
        new_equations = []
        for lhs, rhs in self.deductions.iteritems():
            new_equations.append(sympy.Eq(lhs, rhs))
        return sorted(new_equations, key=lambda x: str(x))

    def print_deduction_log(self):
        ''' Print the judgements and the deductions they have made '''
        for judgement, ded_info in self.deduction_record.iteritems():
            print '\n', judgement
            for eqn, deds in ded_info.iteritems():
                eqn_str = str(eqn).ljust(25)
                ded_str = map(lambda (x, y) : '{}={}'.format(x, y), deds)
                ded_str = ', '.join(ded_str)
                print '{}\t=>\t{}'.format(eqn_str, ded_str)

    def solve_equations(self, max_iter=50, verbose=False):
        ''' Solve a system of equations
        '''
        #prev_ded = len(self.deductions) + 1
        #while len(self.deductions) > prev_ded:
        #    prev_ded = len(self.deductions)

        start_qubits = len(self.variables)

        len_ded = len(self.deductions)
        # The number of iterations in which we've made no new deductions
        num_constant_iter = 0

        for i in xrange(max_iter):
            self.equations = self.clean_equations(self.equations)
            self.apply_judgements(self.equations)
            # If we've run out of equations, we've simplified everything we can
            if len(self.equations) == 0:
                break

            if len(self.deductions) == len_ded:
                num_constant_iter += 1
                if num_constant_iter > 3:
                    break
            else:
                len_ded = len(self.deductions)

        # Final clean and go, for good luck
        self.equations = self.clean_equations(self.equations)
        self.apply_judgements(self.equations)
        self.equations = self.clean_equations(self.equations)



        self.reformulate_equations()

        if verbose:
            self.print_deduction_log()

            print

            print 'Unsimplified equations'
            for e in self.equations:
                print e
            print 'Deductions'
            for e in self.deductions_as_equations:
                print e

            print 'Solns'
            print self.solutions
            print 'Final variables'
            print self.final_variables

            print

            print 'Final Equations'
            for e in sorted(self.final_equations, key=lambda x: str(x)):
                print e

            print

            print 'Num Qubits: {} -> {}'.format(start_qubits,
                                                len(self.final_variables))

#            print 'Final equation'
#            print self.objective_function



    def reformulate_equations(self):
        ''' Reformulate the final equations from the deductions.
            Return the final equations and a dictionary of calculated values
        '''
        solved_var = {}

        # First trawl through the deductions for definite solutions
        for expr, val in self.deductions.iteritems():
#            if not is_one_or_zero(val):
#                continue
            latoms = expr.atoms()
            if len(latoms) == 1:
                solved_var[expr] = val

        # Now trawl through for equations which involve the unsolved variables
        final_equations = []
        unsolved_var = set(self.variables.values()).difference(solved_var.keys())

        for eqn in itertools.chain(self.equations,
                                   self.deductions_as_equations):
            expr, val = eqn.lhs, eqn.rhs
            latoms = expr.atoms()
            ratoms = set([val]) if isinstance(val, int) else val.atoms()
            if len(latoms.intersection(unsolved_var)) or len(ratoms.intersection(unsolved_var)):
                eqn = sympy.Eq(expr, val)
                # Substitute all of the solved variables
                eqn = eqn.subs(solved_var)

                if eqn == True:
                    continue

                eqn = remove_binary_squares(eqn)

                # If the final equation is x=constant, then we're also done!!
#                if len(eqn.atoms(sympy.Symbol)) == 1:
#                    assert eqn.rhs.is_constant()
#                    self.update_value(eqn.lhs, eqn.rhs)
#                    solved_var[eqn.lhs] = eqn.rhs
#                    unsolved_var.difference_update([eqn.lhs])
#                else:
                final_equations.append(eqn)

#        self.final_variables = set()
#        for eq in final_equations:
#            for atom in eq.atoms():
#                if not is_one_or_zero(atom):
#                    self.final_variables.add(atom)

        final_equations = sorted(set(final_equations), key=lambda x: str(x))

        self.final_variables = unsolved_var

        self.final_equations = final_equations
        self.solutions = solved_var

    @property
    def objective_function(self):
        ''' Return the final objective function, using self.final_equations '''
        if self.final_equations is None:
            return None

        obj_func = sympy.sympify(0)
        for eqn in self.final_equations:
            if eqn.rhs == 1:
                obj_func += (eqn.lhs - eqn.rhs) ** 2

        obj_func = obj_func.expand().simplify()

#        # Remove squares
        obj_func = remove_binary_squares(obj_func)

        return obj_func.simplify()

    def objective_function_to_file(self, filename=None):
        ''' Write the objective function to a file, or printing if None.
            Also include the dictionary of variable number to original variable
        '''
        out = expression_to_coef_string(self.objective_function)

        if filename is None:
            print out
        else:
            f = open(filename, 'a')
            f.write(out)
            f.close()

    def clean_equations(self, eqns):
        ''' Remove True equations and simplify '''
        cleaned = []
        for eqn in eqns:

            clean = balance_constant(eqn) #sympy.simplify(eqn)
            clean = divide_2(clean)

            # every square is itself
            clean = remove_binary_squares(clean)

            if clean == True:
                continue

            # Substitute in all values
            clean = clean.subs(self.deductions)

            if clean == True:
                continue

            # every square is itself
            clean = remove_binary_squares(clean)

            cleaned.append(clean)

        # Now add any equations where LHS = RHS1, LHS = RHS2 then RHS1 = RHS2
        for eqn1, eqn2 in itertools.combinations(cleaned, 2):
            if (eqn1.lhs == eqn2.lhs) and (eqn1.rhs != eqn2.rhs):
                cleaned.append(sympy.Eq(eqn1.rhs, eqn2.rhs))
                print 'Equation added! {}, {}'.format(eqn1, eqn2)

        # Match wild variables in two term expressions
        cleaned = self.clean_two_term_sub(cleaned)

        return cleaned


    def update_value(self, expr, value):
        ''' Update the global dictionary and check for contradictions

            >>> system = EquationSolver()
            >>> x = sympy.symbols('x')
            >>> system.update_value(x, 0)
            >>> system.deductions
            {x: 0}
        '''
        # Update the deductions process dictionary
        judgement, eqn = get_judgement()

        #print '{} -> {}'.format(expr, value)

        # If expr = 2*x and value ==0, then we can get rid of the 2
        if value == 0:
            expr = expr.as_coeff_Mul()[1]

        current_val = self.deductions.get(expr)

        # If value already maps to expr, avoid the cycle!
        if self.deductions.get(value) == expr:
            return

        # No possible conflict
        if current_val is None:
            self.deductions[expr] = value
            self.deduction_record[judgement][eqn].append((expr, value))

        # If we already know a value of this family
        elif is_one_or_zero(current_val):
            # If we've found a numeric value
            if is_one_or_zero(value):
                if current_val != value:
                    raise ContradictionException('{} is already set to {} != {}'.format(expr,
                                                 current_val, value))
                else:
                    return
            # We know what we're trying to update to is not a numeric, so update it too
            else:
                self.update_value(value, current_val)
                self.deductions[expr] = current_val
                self.deduction_record[judgement][eqn].append((expr, value))

        # Current_val is symbolic
        else:
            if is_one_or_zero(value):
                # Perform another error check
                cc_val = self.deductions.get(current_val)
                if is_one_or_zero(cc_val) and (cc_val != value):
                    raise ContradictionException(
                            '{} is already set to {} != {}'.format(
                            current_val,
                            cc_val,
                            value))
                self.deductions[current_val] = value
                self.deductions[expr] = value
                self.deduction_record[judgement][eqn].append((current_val, value))
                self.deduction_record[judgement][eqn].append((expr, value))
            # Both values are symbolic!
            else:
                self.deductions[expr] = current_val
                self.deduction_record[judgement][eqn].append((expr, current_val))
                if value != current_val:
                    self.deductions[value] = current_val
                    self.deduction_record[judgement][eqn].append((value, current_val))

    def get_var(self, var):
        ''' Return symbolic variable

            >>> system = EquationSolver()
            >>> x = system.get_var('x')
            >>> x
            x

            >>> x1 = system.get_var('x')
            >>> x1
            x

            >>> x1 is x
            True
        '''
        res = self.variables.get(var)
        if res is None:
            res = sympy.symbols(var, integer=True)
            self.variables[var] = res
        return res

    def set_to_max(self, expr):
        ''' Given an expression, update all terms so that it achieves it's maximum

            >>> system = EquationSolver()
            >>> expr = sympy.sympify('-2*x*y - 3 + 3*z23')
            >>> system.set_to_max(expr)
            >>> system.deductions
            {x*y: 0, z23: 1}
        '''
        coef_dict = expr.as_coefficients_dict()
        for var, coef in coef_dict.iteritems():
            if var == 1:
                continue
            if coef > 0:
                self.update_value(var, 1)
            elif coef < 0:
                self.update_value(var, 0)

    def set_to_min(self, expr):
        ''' Given an expression, update all terms so that it achieves it's minumum

            >>> system = EquationSolver()
            >>> expr = sympy.sympify('-2*x*y - 3 + 3*z23')
            >>> system.set_to_min(expr)
            >>> system.deductions
            {x*y: 1, z23: 0}

        '''
        coef_dict = expr.as_coefficients_dict()
        for var, coef in coef_dict.iteritems():
            if var == 1:
                continue
            if coef < 0:
                self.update_value(var, 1)
            elif coef > 0:
                self.update_value(var, 0)

    def apply_judgements(self, equations):
        ''' Apply judgements to a list of sympy equations and directly update
            self.deductions
        '''
        for eqn in equations:
            self.judgement_0(eqn)
            self.judgement_sq(eqn)
            self.judgement_mm(eqn)
            self.judgement_1(eqn)
            self.judgement_2(eqn)
            self.judgement_4(eqn)
            self.judgement_two_term(eqn)
            self.judgement_7i(eqn)
            self.judgement_7ii(eqn)
            self.judgement_9(eqn)
            self.judgement_10(eqn)
            #self.judgement_10i(eqn)

    def judgement_0(self, eqn):
        ''' Add x=y to deductions '''
        if len(eqn.lhs.atoms()) == len(eqn.rhs.atoms()) == 1:
            if eqn.lhs.is_constant():
                self.update_value(eqn.rhs, eqn.lhs)
            else:
                self.update_value(eqn.lhs, eqn.rhs)

    def judgement_sq(self, eqn):
        ''' If RHS is a non-zero constant and the LHS is a product of variables,
            then the variables must all be 1
            x*y*z=1 => x = y = 1

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x * y * z, 1)
            >>> system.judgement_sq(eqn)
            >>> system.deductions
            {x: 1, z: 1, y: 1}
        '''
        if eqn.rhs == 1:
            if len(eqn.lhs.as_ordered_terms()) == 1:
                for var in eqn.lhs.atoms(sympy.Symbol):
                    self.update_value(var, 1)

    def judgement_two_term(self, eqn):
        ''' If an expression has 2 variable terms, sub it in!
            This adds lots of complexity so the infrastructure, which is
            why we don't do it with 3 or 4 term sums

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y*z, 1)
            >>> system.judgement_two_term(eqn)
            >>> system.deductions
            {x + y*z: 1}
        '''
        if len(eqn.lhs.as_ordered_terms()) == 2 and eqn.rhs.is_constant():
            self.update_value(eqn.lhs, eqn.rhs)

    def clean_two_term_sub(self, eqns):
        ''' If an expression has 2 variable terms, sub it in!
            This adds lots of complexity so the infrastructure, which is
            why we don't do it with 3 or 4 term sums
            NOTE doesn't work fully, and will probably be deprecated

            >>> system = EquationSolver()
            >>> x, y, z, a, b = sympy.symbols('x y z a b')
            >>> system.deductions[x + y] = 1

            >>> eqn = sympy.Eq(x*z + y*z, a*b)
            >>> system.clean_two_term_sub([eqn])
            [z == a*b]

            >>> eqn = sympy.Eq(x*z + y*z + 3, a*b)
            >>> system.clean_two_term_sub([eqn])
            [x*z + y*z + 3 == a*b]

            >>> eqn = sympy.Eq(x*z + 2*y*z, a*b)
            >>> system.clean_two_term_sub([eqn])
            [x*z + 2*y*z == a*b]
        '''
        # Wild character we use for multiplication
        wild = sympy.Wild('w')
        to_replace = {}
        for expr, val in self.deductions.iteritems():
            if (num_add_terms(expr) == 2) and (val == 1):
                to_replace[(expr*wild).expand()] = wild

        def _helper_replace(eqn):
            # If we have too many variables, return so we don't take an age
            if len(eqn.atoms(sympy.Symbol)) > 10:
                return eqn
            for k, v in to_replace.iteritems():
                if eqn == True:
                    return True
                eqn = eqn.replace(k, v, exact=False)
            return eqn

        eqns = map(_helper_replace, eqns)
        eqns = filter(lambda x: x != True, eqns)
        return eqns


    def judgement_mm(self, eqn):
        ''' If min(rhs) == max(lhs), then we know what to do

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y + z, 3)
            >>> system.judgement_mm(eqn)
            >>> system.deductions
            {x: 1, z: 1, y: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 2*y, 5 - 2*z)
            >>> system.judgement_mm(eqn)
            >>> system.deductions
            {x: 1, z: 1, y: 1}
        '''

        def _helper(self, lhs, rhs):
            if min_value(lhs) == max_value(rhs):
                self.set_to_min(lhs)
                self.set_to_max(rhs)

        _helper(self, eqn.lhs, eqn.rhs)
        _helper(self, eqn.rhs, eqn.lhs)

    def judgement_1(self, eqn):
        ''' If x + y + z = 1 then xy = yz = zx = 0
            Generally true for any number of terms that = 1

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y + z, 1)
            >>> system.judgement_1(eqn)
            >>> system.deductions
            {x*z: 0, x*y: 0, y*z: 0}

            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + y + z + z2, 1)
            >>> system.judgement_1(eqn)
            >>> system.deductions
            {x*z: 0, z*z2: 0, y*z2: 0, x*z2: 0, y*z: 0, x*y: 0}
        '''
        if eqn.rhs != 1:
            return

        terms = eqn.lhs.as_ordered_terms()
        # We need more than 2 terms
        if len(terms) < 2:
            return

        # We also need every term to be positive
        for term in terms:
            if min_value(term) < 0:
                return

        pairs = itertools.combinations(terms, 2)
        for t1, t2 in pairs:
            self.update_value(t1*t2, 0)


#        if len(terms) == 3:
#            variables, coef = zip(*terms)
#            if all([c == 1 for c in coef]):
#                for v in itertools.permutations(variables, 2):
#                    self.update_value(v[0] * v[1], 0)


    def judgement_2(self, eqn):
        ''' If max(lhs) < max(rhs) and there is only 1 variable term
            in the RHS, then this term must be 1

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x+y, 2*z + 1)
            >>> system.judgement_2(eqn)
            >>> system.deductions
            {z: 0}
        '''
        lhs, rhs = eqn.lhs, eqn.rhs
        if (max_value(lhs) < max_value(rhs)):
            terms = rhs.as_ordered_terms()
            if (len(terms) == 2) and terms[-1].is_constant():
                self.set_to_min(rhs)

    def judgement_4(self, eqn):
        ''' If max(lhs) = max(rhs) and rhs is constant, then every term on the
            left is 1

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2)
            >>> system.judgement_4(eqn)
            >>> system.deductions
            {x: 1, y: 1}
        '''
        lhs, rhs = eqn.lhs, eqn.rhs
        if (max_value(lhs) == max_value(rhs)) and rhs.is_constant():
            self.set_to_max(lhs)

    def judgement_7i(self, eqn):
        ''' If a term being 1 would tip max(rhs) > max(lhs), then it must be 0

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2 + z)
            >>> system.judgement_7i(eqn)
            >>> system.deductions
            {z: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2 + z)
            >>> system.judgement_7i(eqn)
            >>> system.deductions
            {z: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 3 * y, 2 + z)
            >>> system.judgement_7i(eqn)
            >>> system.deductions
            {}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 5 * y, 2 + z)
            >>> system.judgement_7i(eqn)
            >>> system.deductions
            {y: 0}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 1 + z)
            >>> system.judgement_7i(eqn)
            >>> system.deductions
            {}
        '''
        eqn = balance_constant(eqn)

        def _helper(lhs, rhs):
            lhs_max = max_value(lhs)

            rhs_terms = rhs.as_ordered_terms()
            if rhs_terms[-1].is_constant():
                rhs_const = rhs_terms.pop()
            else:
                rhs_const = 0

            for term in rhs_terms:
                if max_value(term) + rhs_const > lhs_max:
                    self.set_to_min(term)

        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)

    def judgement_7ii(self, eqn):
        ''' x + 1 = 2y -> y = 1

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 1, 2*z)
            >>> system.judgement_7ii(eqn)
            >>> system.deductions
            {z: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x * y + 1, 2 * z)
            >>> system.judgement_7ii(eqn)
            >>> system.deductions
            {z: 1}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + 3 * y + 1, 2 * z + 1)
            >>> system.judgement_7ii(eqn)
            >>> system.deductions
            {}
        '''
        if ((eqn.lhs.as_coeff_add()[0] == 1) and
            (len(eqn.lhs.as_ordered_terms()) > 1) and
            (eqn.rhs.as_coeff_mul()[0] == 2) and
            (len(eqn.rhs.as_ordered_terms()) == 1)):
            self.set_to_max(eqn.rhs)

    def judgement_9(self, eqn):
        ''' x + y = 2z -> x = y = z

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y, 2*z)
            >>> system.judgement_9(eqn)
            >>> system.deductions
            {z: x, y: x}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x * y + x, 2 * z)
            >>> system.judgement_9(eqn)
            >>> system.deductions
            {x: x*y, z: x*y}
        '''
        lhs, rhs = eqn.lhs, eqn.rhs
        if len(lhs.as_ordered_terms()) == 2:
            t1, t2 = lhs.as_ordered_terms()
            if ((t1.as_coeff_mul()[0] == 1) and
                (not t1.is_constant()) and
                (t2.as_coeff_mul()[0] == 1) and
                (not t2.is_constant()) and
                (rhs.as_coeff_mul()[0] == 2)):
                self.update_value(t2, t1)
                self.update_value(rhs / 2, t1)

    def judgement_10(self, eqn):
        ''' x + 2y + z = 4z2 -> x = z
            Also works with any even RHS and any number of even variables on
            the LHS.

            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + 2*y + z, 2*z2)
            >>> system.judgement_10(eqn)
            >>> system.deductions
            {x: z}

            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + 2*y + 3*z, 2 * z2)
            >>> system.judgement_10(eqn)
            >>> system.deductions
            {x: 3*z}
        '''

        def _helper(lhs, rhs):
            for term in rhs.as_ordered_terms():
                if term.as_coeff_mul()[0] % 2:
                    return

            odd_terms = []
            for term in lhs.as_ordered_terms():
                term_coef = term.as_coeff_mul()[0]
                if (term_coef % 2):
                    if len(odd_terms) > 2:
                        return
                    else:
                        odd_terms.append(term)

            if len(odd_terms) == 2:
                if odd_terms[0].is_constant():
                    self.update_value(odd_terms[1], odd_terms[0])
                else:
                    self.update_value(*odd_terms)

        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)

    def judgement_10i(self, eqn):
        ''' x + y + z + 2a = 2b -> xy + xz = yx + yz = zx + zy = 0
            Also works with any even RHS and any number of even variables on
            the LHS.

            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + 2*y + z, 2*z2)
            >>> system.judgement_10(eqn)
            >>> system.deductions
            {x: z}

            >>> system = EquationSolver()
            >>> x, y, z, z2 = sympy.symbols('x y z z2')
            >>> eqn = sympy.Eq(x + 2*y + 3*z, 2 * z2)
            >>> system.judgement_10(eqn)
            >>> system.deductions
            {x: 3*z}

            >>> system = EquationSolver()
            >>> x, y, z = sympy.symbols('x y z')
            >>> eqn = sympy.Eq(x + y + 1, 2 * z)
            >>> system.judgement_10(eqn)
            >>> system.deductions
            {}
        '''

        def _helper(lhs, rhs):
            for term in rhs.as_ordered_terms():
                if term.as_coeff_mul()[0] % 2:
                    return

            odd_terms = []
            for term in lhs.as_ordered_terms():
                term_coef = term.as_coeff_mul()[0]
                if (term_coef % 2):
                    if len(odd_terms) > 3:
                        return
                    else:
                        if term.is_constant():
                            return
                        odd_terms.append(term)

            if len(odd_terms) == 3:
                print '\n*** 10i\n{}'.format(eqn)
                pairs = [a * b for a, b in itertools.combinations(odd_terms, 2)]
                expressions = [ab + bc for ab, bc in itertools.combinations(pairs, 2)]
                for expr in expressions:
                    print expr
                    self.update_value(expr, 0)


        _helper(eqn.lhs, eqn.rhs)
        _helper(eqn.rhs, eqn.lhs)


    ## Simulation


    def simplified_system(self):
        ''' Return a new EquationSolver with the same parameters, with new
            deductions
        '''
        new_equations = self.final_equations[:]
        new_variables = {str(v) : v for v in self.final_variables}
        return EquationSolver(new_equations, new_variables)

## Parsing and set up

def parse_equations(params):
    '''
        params[0][i] contains all leftsides
        params[1][i] contains all rightsides
    '''
    eqns = []
    variables = {}
    for lhs, rhs in itertools.izip(*params):
        lhs = _parse_expression(lhs, variables)
        rhs = _parse_expression(rhs, variables)
        eqn = sympy.Eq(lhs, rhs)
        eqns.append(eqn)

    return eqns, variables

def _parse_expression(expr, variables):
    ''' Take a list of terms, clean them up and return a sympy expression in terms
        of the global variables '''
    out = 0
    for term in expr:
        if term == '':
            continue
        out += _parse_term(term, variables)
    return out

def _parse_term(term, variables):
    ''' Take a term and clean it, replacing variables '''
    coef = int(ReHandler.get_coefficient(term))
    var = ReHandler.get_variables(term)
    var_instances = []
    for v in var:
        instance = variables.get(v)
        if instance is None:
            instance = sympy.symbols(v)
            variables[v] = instance
        var_instances.append(instance)
    return coef * sympy.prod(var_instances)

## Helper functions

def is_one_or_zero(val):
    ''' Self explanatory '''
    return (val == 0) or (val == 1)

def num_add_terms(expr, check=False):
    ''' Return the number of additive terms in an expression.
        Note doesn't work for multiplicative terms!

        >>> expr = 'x + 2*y'
        >>> num_add_terms(sympy.sympify(expr))
        2
        >>> expr = 'x + 2*y + 3*z + 5'
        >>> num_add_terms(sympy.sympify(expr))
        4
        >>> expr = 'x * 2*y'
        >>> num_add_terms(sympy.sympify(expr))
        1
        >>> expr = '(x + 2*y) * z'
        >>> num_add_terms(sympy.sympify(expr), check=False)
        1
    '''
    if check:
        assert expr.func == sympy.Add
    return len(expr.as_ordered_terms())

def is_constant(expr):
    ''' Determine whether an expression is constant
        >>> expr = 'x + 2*y'
        >>> is_constant(sympy.sympify(expr))
        False
        >>> expr = 'x + 5'
        >>> is_constant(sympy.sympify(expr))
        False
        >>> expr = '3'
        >>> is_constant(sympy.sympify(expr))
        True
        >>> expr = '2*x - 4'
        >>> is_constant(sympy.sympify(expr))
        False
    '''
    return len(expr.atoms(sympy.Symbol)) == 0

def divide_2(eqn):
    ''' Return equation divided by 2 if needed '''
    if eqn == True:
        return True
    for term in eqn.lhs.as_ordered_terms():
        if (term.as_coeff_mul()[0] % 2) == 1:
            return eqn
    return sympy.Eq(eqn.lhs / 2, eqn.rhs / 2)

def balance_constant(eqn):
    ''' Take an equation and tidy up the constant part

        >>> lhs = sympy.sympify('x + 2')
        >>> rhs = sympy.sympify('y + 3')
        >>> balance_constant(sympy.Eq(lhs, rhs))
        x == y + 1

        >>> lhs = sympy.sympify('x - 2')
        >>> rhs = sympy.sympify('y + 3')
        >>> balance_constant(sympy.Eq(lhs, rhs))
        x == y + 5

        >>> lhs = sympy.sympify('x + 2')
        >>> rhs = sympy.sympify('y - 3')
        >>> balance_constant(sympy.Eq(lhs, rhs))
        x + 5 == y
    '''
    if eqn == True:
        return True

    lhs_c = eqn.lhs.as_coeff_add()[0]
    rhs_c = eqn.rhs.as_coeff_add()[0]

    if (lhs_c < 0) or (rhs_c < 0):
        to_add = abs(min(lhs_c, rhs_c))
    else:
        to_add = - min(lhs_c, rhs_c)

    return sympy.Eq(eqn.lhs + to_add,
                    eqn.rhs + to_add)


def max_value(expr):
    ''' Return the max value of an expression

        >>> expr = sympy.sympify('x+1')
        >>> max_value(expr)
        2

        >>> expr = sympy.sympify('x-1')
        >>> max_value(expr)
        0

        >>> expr = sympy.sympify('x + y + 1')
        >>> max_value(expr)
        3

        >>> expr = sympy.sympify('2*x - y + 1')
        >>> max_value(expr)
        3

        >>> expr = sympy.sympify('x*y*z + a*b')
        >>> max_value(expr)
        2

        >>> expr = sympy.sympify('-2 * (x+1)')
        >>> max_value(expr)
        -2
    '''
    #expr = sympy.simplify(expr)
    coef = expr.as_coefficients_dict()
    max_ = 0
    for term, c in coef.iteritems():
        if c > 0 or term == 1:
            max_ += c
    return max_

def min_value(expr):
    ''' Return the min value of an expression

        >>> expr = sympy.sympify('x+1')
        >>> min_value(expr)
        1

        >>> expr = sympy.sympify('x-1')
        >>> min_value(expr)
        -1

        >>> expr = sympy.sympify('x + y +1')
        >>> min_value(expr)
        1

        >>> expr = sympy.sympify('2*x - y + 1')
        >>> min_value(expr)
        0

        >>> expr = sympy.sympify('x*y*z + a*b - 4*k + 2')
        >>> min_value(expr)
        -2

        >>> expr = sympy.sympify('-2 * (x+1)')
        >>> min_value(expr)
        -4
    '''
    #expr = sympy.simplify(expr)
    coef = expr.as_coefficients_dict()
    min_ = 0
    for term, c in coef.iteritems():
        if c < 0 or term == 1:
            min_ += c
    return min_

def remove_binary_squares(eqn):
    ''' Given an equation, remove all of the squares as any binary
        variable squared is itself.
        NOTE: This doesn't work at all for higher powers

        >>> expr = 'x**2 * y + z**3 + 2*z - 4'
        >>> remove_binary_squares(sympy.sympify(expr))
        x*y + z**3 + 2*z - 4
    '''
    if eqn == True:
        return True
    for var in eqn.atoms(sympy.Symbol):
        eqn = eqn.subs(var ** 2, var)
    return eqn

## Objective function help

def expression_to_coef_string(expr):
    ''' Write the objective function to a file, or printing if None.
        Also include the dictionary of variable number to original variable

        >>> inp = '40*s_1 + 30*s_1*s_2 + 100*s_1*s_2*s_3 - 15*s_2*s_3 - 20*s_3 + 4'
        >>> print expression_to_coef_string(inp)
        4
        1 2 30
        2 3 -15
        1 40
        3 -20
        1 2 3 100
        <BLANKLINE>
        {s_3: 3, s_2: 2, s_1: 1}
    '''

    if isinstance(expr, str):
        expr = sympy.sympify(expr)

    expr = expr.expand().simplify()
    atoms = expr.atoms(sympy.Symbol)
    atoms = sorted(list(atoms), key=lambda x: str(x))
    atom_map = {v : i + 1 for i, v in enumerate(atoms)}

    lines = []
    for term, coef in expr.as_coefficients_dict().iteritems():
        var_num = sorted([atom_map[atom] for atom in term.atoms(sympy.Symbol)])
        line = ' '.join(map(str, var_num))
        if line:
            line += ' ' + str(coef)
        else:
            line = str(coef)
        lines.append(line)

    coef_str = '\n'.join(lines)
    atom_str = str(atom_map)
    return '\n\n'.join([coef_str, atom_str])


# Inspection
def get_judgement():
    ''' Find the judgement calling update_value. Horrible hackery, but at least
        it's in 1 place
    '''
    up_again = ['set_to_min', 'set_to_max', 'update_value', '_helper']

    ind = 2
    frame = inspect.stack()[ind]
    caller_name = frame[-3]
    while caller_name in up_again:
        ind += 1
        frame = inspect.stack()[ind]
        caller_name = frame[-3]

    eqn = frame[0].f_locals.get('eqn')
    return caller_name, eqn


if __name__ == "__main__":
    import doctest
    doctest.testmod()
