'''FF2SE: tools for using FiniteFlow with SymEngine.

FF2SE includes some utilities to facilitate the combined usage of
FiniteFlow and SymEngine.

Note that this is NOT a full SymEngine <-> FiniteFlow interface.

'''


import symengine as _s
from collections import defaultdict as _defdict
import ff2cas as _ff2cas
import _ff2symse_common
from fractions import Fraction as _Fraction


# Utilities

def CoefficientRules(poly,variables,varset=None):
    Integer = _s.Integer
    if varset is None:
        varset = set(variables)
    p = _s.sympify(poly).expand().as_coefficients_dict()
    crules = _defdict(lambda : Integer(0))
    for term,cc in p.items():
        if not cc == 0:
            termdict = term.as_powers_dict()
            ruleterm = Integer(1)
            rulecc = cc
            for b,e in termdict.items():
                if b in varset:
                    ruleterm *= b**e
                else:
                    rulecc *= b**e
            ruleterm = tuple(termdict[v] for v in variables)
            crules[ruleterm] += rulecc
    return dict(crules)

def FromCoefficientRules(polyrules,variables):
    Mul = _s.Mul
    Add = _s.Add
    return Add(*(cc * Mul(*(v**e for v,e in zip(variables,exps))) \
                 for exps,cc in polyrules.items()))

def LCoefficientList(expr,varset):
    Integer = _s.Integer
    if not type(varset) is set:
        raise TypeError("Variables must be passed as a set")
    if expr.is_Equality:
        if len(expr.args) != 2:
            raise ValueError("Equalities must have 2 arguments")
        exprl = zip(expr.args,(Integer(-1),Integer(1)))
    else:
        exprl = ((expr,Integer(1)),)
    crules = _defdict(lambda : Integer(0))
    constterm = Integer(0)
    for (arg,pref) in exprl:
        for term,cc in arg.as_coefficients_dict().items():
            if not cc == 0:
                termdict = term.as_powers_dict()
                thisvar = None
                rulecc = cc
                for b,e in termdict.items():
                    if b in varset:
                        if e==1 and thisvar is None:
                            thisvar = b
                        else:
                            raise ValueError("Input not linear in the variables")
                    else:
                        rulecc *= b**e
                if thisvar is None:
                    constterm += pref*rulecc
                else:
                    crules[thisvar] += pref*rulecc
    return (dict(crules),constterm)


def _ToPolyData(variables,poly,varset,applyfun):
    return list((str(applyfun(c)),e) for e,c in \
                CoefficientRules(poly,variables,varset=varset).items())

def Together(expr):
    '''Put expression under common denominator.'''
    return _ff2symse_common.Together(_s,expr)

def FastNumerDenom(expr):
    '''Return numerator and denominator, assuming the input is a
    product of polynomials raised to integer powers.'''
    return _ff2symse_common.FastNumerDenom(_s,expr)

def _ToFunData(variables,fun,varset,applyfun):
    #num,den = fun.as_numer_denom()
    num,den = FastNumerDenom(fun)
    return (_ToPolyData(variables,num,varset,applyfun),
            _ToPolyData(variables,den,varset,applyfun))


# Implementation of CAS object

class _CAS(_ff2cas.CAS):

    def __init__(self):
        self.Add = _s.Add
        self.Mul = _s.Mul

    def Rational(self,numstr):
        rat = _Fraction(numstr)
        return _s.Rational(rat.numerator, rat.denominator)

    def LCoefficientLists(self,variables,eqs):
        varset = set(variables)
        for eq in eqs:
            a,b = LCoefficientList(eq, varset)
            yield (a,-b)

    def ToRatFunListData(self,variables,functions,applyfun=_ff2cas.identity):
        varset = set(variables)
        return list(_ToFunData(variables,fun,varset,applyfun) \
                    for fun in functions)

    def VariablesIn(self,expr):
        return expr.free_symbols

    def IsZero(self,expr):
        return expr.is_zero


_cas = _CAS()



# Definition of CAS functions

def ToRatFunList(variables, functions):
    '''Convert a list of rational functions from SymEngine into a
    fflow.RatFunList

    '''
    return _cas.ToRatFunList(variables, functions)

def ToIdxRatFunList(variables, functions, indexes=None):
    '''Convert a list of rational functions from SymEngine into a
    fflow.IdxRatFunList

    '''
    return _cas.ToIdxRatFunList(variables, functions, indexes)

def FromRatFunList(variables, functions):
    '''Convert a fflow.IdxRatFunList into a list of SymEngine rational
    functions.'''
    return _cas.FromRatFunList(variables, functions)

def FromRatNumList(nums):
    '''Convert a list strings into SymEngine rational numbers.'''
    return _cas.FromRatFunList(num)

def ToRatNumList(nums):
    '''Convert a list of SymEngine rational numbers into strings.'''
    return _cas.ToRatNumList(nums)

def ToRatFunListFromCoeffs(coeffs, variables, functions):
    '''Convert a list of rational functions with symbolic coefficients
    from SymEngine into a fflow.RatFunList which can be used in
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
    '''ToNumericSparseLSolve(eqs, variables) returns (non_zero_els,
    non_zero_coeffs) defined as the inputs required by
    fflow.AlgNumericSparseLSolve.

    '''
    return _cas.ToNumericSparseLSolve(eqs, variables)

def LSolveDict(graph,node,variables,graphout):
    '''LSolveDict(graph,node,variables,graphout) where node is a
    linear solver node of the specified graph and graphout is its
    output converted to SymEngine, returns a dictionary with the
    solution.  The keys of the dictionary are the dependent variables,
    i.e. those for which a solution was found.

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
    coefficients ccs, into SymEngine.
    '''
    return _cas.FromLaurentOutput(graph,node,ccs,expvar)

def FromCoeffsAndExponents(coefficients,variables,exponents):
    '''Return a list of rational functions from a list of
    coefficients, variables and exponent data.'''
    return _cas.FromCoeffsAndExponents(coefficients,variables,exponents)

def AlgToExprEval(graph, in_node, expressions, variables):
    return _ff2symse_common.AlgToExprEval(_s, graph, in_node,
                                          expressions, variables)
