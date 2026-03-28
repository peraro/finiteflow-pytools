'''FF2Sym: tools for using FiniteFlow with SymPy.

FF2Sym includes some utilities to facilitate the combined usage of
FiniteFlow and SymPy.

Note that this is NOT a full SymPy <-> FiniteFlow interface.

'''

import sympy as _s
from sympy.polys.matrices.linsolve import _linear_eq_to_dict as _leq2dict
import ff2cas as _ff2cas
import _ff2symse_common


# Implementation of CAS object

def Together(expr):
    '''Put expression under common denominator.'''
    return _ff2symse_common.Together(_s,expr)

def FastNumerDenom(expr):
    '''Return numerator and denominator, assuming the input is a
    product of polynomials raised to integer powers.'''
    return _ff2symse_common.FastNumerDenom(_s,expr)

def _Poly2FFCRules(variables, polydict, applyfun):
    return list((str(applyfun(val)),exps) for exps,val in polydict.items())

def _Sym2RatFun(variables, expr, applyfun):
    num,den = FastNumerDenom(expr)
    num = _s.Poly(num, *variables).as_dict()
    den = _s.Poly(den, *variables).as_dict()
    return (_Poly2FFCRules(variables,num,applyfun),
            _Poly2FFCRules(variables,den,applyfun))


class _CAS(_ff2cas.CAS):

    def __init__(self):
        self.Rational = _s.Rational
        self.Add = _s.Add
        self.Mul = _s.Mul

    def ToRatFunListData(self,variables,functions,applyfun=_ff2cas.identity):
        return list(_Sym2RatFun(variables, expr, applyfun) \
                    for expr in functions)

    def LCoefficientLists(self,variables, eqs):
        a, b = _leq2dict(eqs, variables)
        b = list(-el for el in b)
        return zip(a,b)

    def AsSeries(self,expr, x, x0, expmax):
        return expr + _s.O(x**(expmax+1),(x,x0))

    def VariablesIn(self,expr):
        return expr.free_symbols


_cas = _CAS()


# Definition of CAS functions

def ToRatFunList(variables, functions):
    '''Convert a list of rational functions from SymPy into a
    fflow.RatFunList'''
    return _cas.ToRatFunList(variables, functions)

def ToIdxRatFunList(variables, functions, indexes=None):
    '''Convert a list of rational functions from SymPy into a
    fflow.IdxRatFunList'''
    return _cas.ToIdxRatFunList(variables, functions, indexes)

def FromRatFunList(variables, functions):
    '''Convert a fflow.IdxRatFunList into a list of SymPy rational
    functions.'''
    return _cas.FromRatFunList(variables, functions)

def FromRatNumList(nums):
    '''Convert a list strings into SymPy rational numbers.'''
    return _cas.FromRatFunList(num)

def ToRatNumList(nums):
    '''Convert a list of SymPy rational numbers into strings.'''
    return _cas.ToRatNumList(nums)

def ToRatFunListFromCoeffs(coeffs, variables, functions):
    '''Convert a list of rational functions with symbolic coefficients
    from SymPy into a fflow.RatFunList which can be used in
    fflow.AlgRatFunEvalFromCoeffs.  The argument `coeffs` must be
    either a list or tuple of coefficients (in the same order they
    will be in the input node) or a function that converts each
    coefficient to the correct index.
    '''
    return _cas.ToRatFunListFromCoeffs(coeffs, variables, functions)

def ToAnalyticSparseLSolve(params, eqs, variables):
    '''ToAnalyticSparseLSolve(params, eqs, variables) returns
    (non_zero_els, non_zero_coeffs) defined as the inputs required by
    fflow.AlgAnalyticSparseLSolve.

    '''
    return _cas.ToAnalyticSparseLSolve(params, eqs, variables)

def ToNumericSparseLSolve(eqs, variables):
    '''ToNumericSparseLSolve(eqs, variables) returns
    (non_zero_els, non_zero_coeffs) defined as the inputs required by
    fflow.AlgNumericSparseLSolve.

    '''
    return _cas.ToNumericSparseLSolve(eqs, variables)

def LSolveDict(graph,node,variables,graphout):
    '''LSolveDict(graph,node,variables,graphout) where node is a
    linear solver node of the specified graph and graphout is its
    output converted to SymPy, returns a dictionary with the solution.
    The keys of the dictionary are the dependent variables, i.e. those
    for which a solution was found.

    '''
    return _cas.LSolveDict(graph,node,variables,graphout)

def SparseSolve(eqs, variables,
                parameters=None,
                indep_vars_only=False, mark_and_sweep=True,
                sparse_output = True, only_non_homogeneous = False,
                start_mod=0, max_primes=0, max_deg=0, dbginfo=0, polymethod=0,
                n_threads=0, needed_vars=None):
    '''SparseSolve(eqs,variables) returns a dictionary with the
    solution of the linear system.

    '''
    return _cas.SparseSolve(eqs, variables,
                            parameters=parameters,
                            indep_vars_only=indep_vars_only,
                            mark_and_sweep=mark_and_sweep,
                            sparse_output=sparse_output,
                            only_non_homogeneous=only_non_homogeneous,
                            start_mod=start_mod, max_primes=max_primes,
                            max_deg=max_deg, dbginfo=dbginfo,
                            polymethod=polymethod,
                            n_threads=n_threads, needed_vars=needed_vars)

def ToAnalyticDenseLSolve(params, eqs, variables):
    '''Returns coeffs defined as the inputs required by
    fflow.AlgAnalyticDenseLSolve.'''
    return _cas.ToAnalyticDenseLSolve(params, eqs, variables)

def ToNumericDenseLSolve(eqs, variables):
    '''Returns coeffs defined as the inputs required by
    fflow.AlgNumericDenseLSolve.'''
    return _cas.ToNumericDenseLSolve(eqs, variables)

def DenseSolve(eqs, variables, parameters=None,
               indep_vars_only=False, only_non_homogeneous = False,
               start_mod=0, max_primes=0, max_deg=0, dbginfo=0, polymethod=0,
               n_threads=0, needed_vars=None):
    '''Returns a dictionary with the solution of the linear system.'''
    return _cas.DenseSolve(eqs, variables, parameters=parameters,
                           indep_vars_only=indep_vars_only,
                           only_non_homogeneous = only_non_homogeneous,
                           start_mod=start_mod, max_primes=max_primes,
                           max_deg=max_deg, dbginfo=dbginfo,
                           polymethod=polymethod,
                           n_threads=n_threads, needed_vars=needed_vars)

def FromLaurentOutput(graph,node,ccs,expvar):
    '''Convert the output of a Laurent expansion node, i.e. the
    coefficients ccs, into SymPy.
    '''
    return _cas.FromLaurentOutput(graph,node,ccs,expvar)

def FromCoeffsAndExponents(coefficients,variables,exponents):
    '''Return a list of rational functions from a list of
    coefficients, variables and exponent data.'''
    return _cas.FromCoeffsAndExponents(coefficients,variables,exponents)

def AlgToExprEval(graph, in_node, expressions, variables):
    return _ff2symse_common.AlgToExprEval(_s, graph, in_node,
                                          expressions, variables)


# These are extra utilities to convert to/from Mathematica expressions

from sympy.parsing.mathematica import parse_mathematica as _parse_math

def MathToSym(exprstr):
    '''Convert a string containing a Mathematica expression to SymPy.'''
    return _parse_math(exprstr)

def SymToMath(expr):
    '''Convert a SymPy expression into a string containing a
    Mathematica expression.'''
    return _s.mathematica_code(expr)
