'''Generic implementation of tools for communicating between FiniteFlow
and a Computer Algebra System (CAS) in Python.

Specific implementations must inherit from the CAS class.
'''

import fflow as _ff
from itertools import islice as _islice, chain as _chain

identity = lambda x : x

class CAS:
    '''To implement conversions between FiniteFlow and a CAS, derive
    from this class and provide the required methods.  On top of
    these, expressions returned by the CAS must overload arithmetic
    operators.

    '''


    # Required methods to be implemented by derived classes

    def Rational(self,numstr):
        '''Rational numbers from strings'''
        raise NotImplementedError()

    def Add(self,*args):
        raise NotImplementedError()

    def Mul(self,*args):
        raise NotImplementedError()

    def ToRatFunListData(self,variables,expr,applyfun=identity):
        '''Convert a list of rational functions to Finiteflow monomial
        data.  The last argument `applyfun`, if provided, is a
        function to be applied to the polynomial coefficients before
        converting them to strings.
        '''
        raise NotImplementedError()

    def LCoefficientLists(self,varslist,expr):
        '''Returns an iterable of tuples containing a dictionary of
           linear coefficients and a constant term for each expression
           (which may be a linear combination or an equation).
        '''
        raise NotImplementedError()

    def VariablesIn(self,expr):
        '''Returns a list of variables appearing in the expression'''
        raise NotImplementedError()


    # Optional methods

    def AsSeries(self,expr,x,x0,expmax):
        '''Turns a Laurent polynomial into a series.  Can be optionally
        overridden by derived classes.'''
        return expr


    # Provided methods

    def ToRatFunList(self, variables, functions):
        '''Convert a list of rational functions from the CAS into a
        fflow.RatFunList

        '''
        fundata = self.ToRatFunListData(variables, functions)
        return _ff.NewRatFunList(len(variables), fundata)

    def ToIdxRatFunList(self, variables, functions, indexes=None):
        '''Convert a list of rational functions from the CAS into a
        fflow.IdxRatFunList

        '''
        if indexes is None:
            (unique, indexes) = _ff._getUniqueAndIdx(functions)
        else:
            unique = functions
        ret0 = self.ToRatFunList(variables, unique)
        return _ff.MoveRatFunToIdx(ret0, indexes)

    def ToRatFunListFromCoeffs(self, coeffs, variables, functions):
        '''Convert a list of rational functions with symbolic
        coefficients from the CAS into a fflow.RatFunList which can be
        used in fflow.AlgRatFunEvalFromCoeffs.  The argument `coeffs`
        must be either a list or tuple of coefficients (in the same
        order they will be in the input node) or a function that
        converts each coefficient to the correct index.
        '''
        if type(coeffs) is list or type(coeffs) is tuple:
            ccdict = dict(zip(coeffs,range(len(coeffs))))
            applyfun = lambda x : ccdict[x]
        else:
            applyfun = coeffs
        fundata = self.ToRatFunListData(variables, functions, applyfun)
        return _ff.NewRatFunList(len(variables), fundata)

    def _PolyData2Sym(self, variables, polydata):
        Mul = self.Mul
        Add = self.Add
        Rational = self.Rational
        return Add(*(Mul(Rational(cc),*(v**e for v,e in zip(variables,exps))) \
                     for cc,exps in polydata))

    def _RatFun2Sym(self, variables, num, den):
        return (self._PolyData2Sym(variables,num)/
                self._PolyData2Sym(variables,den))

    def FromRatFunList(self, variables, functions):
        '''Convert a fflow.IdxRatFunList into a list of rational
        functions in the CAS.

        '''
        fundata = functions.monomials()
        return list(self._RatFun2Sym(variables,num,den) for num,den in fundata)

    def FromRatNumList(self, nums):
        '''Convert a list strings into rational numbers of the CAS.'''
        return list(self.Rational(el) for el in nums)

    def ToRatNumList(self, nums):
        '''Convert a list of rational numbers from the CAS into strings.'''
        return list(str(el) for el in nums)


    def _ToLSolveInput(self, getclist, eqs, variables):
        var2col = dict(zip(variables,range(len(variables))))
        coeffs = []
        cols = []
        for lhs,rhs in self.LCoefficientLists(variables, eqs):
            thisrowcols = []
            for x,cc in lhs.items():
                thisrowcols.append((var2col[x],cc))
            if not rhs.is_zero:
                thisrowcols.append((len(variables),rhs))
            thisrowcols.sort(key = lambda x : x[0])
            cols.append(list(cc for cc,_ in thisrowcols))
            coeffs.extend(cc for _,cc in thisrowcols)
        coeffs  = getclist(coeffs)
        return (cols, coeffs)

    def ToAnalyticSparseLSolve(self, params, eqs, variables):
        '''Returns (non_zero_els, non_zero_coeffs) defined as the inputs
        required by fflow.AlgAnalyticSparseLSolve.

        '''
        def getclist(ccs):
            return self.ToIdxRatFunList(params, ccs)
        return self._ToLSolveInput(getclist, eqs, variables)

    def ToNumericSparseLSolve(self, eqs, variables):
        '''Returns (non_zero_els, non_zero_coeffs) defined as the inputs
        required by fflow.AlgNumericSparseLSolve.

        '''
        return self._ToLSolveInput(self.ToRatNumList, eqs, variables)


    def LSolveDict(self,graph,node,variables,graphout):
        '''LSolveDict(graph,node,variables,graphout) where node is a
        sparse solver node of the specified graph and graphout is its
        output converted to the CAS, returns a dictionary with the
        solution.  The keys of the dictionary are the dependent
        variables, i.e. those for which a solution was found.

        '''
        sol = dict()

        depvars = _ff.LSolveDepVars(graph,node)
        depvars = list(variables[idx] for idx in depvars)
        idx=0
        extvars = variables + [1]
        if _ff.LSolveOutputIsSparse(graph,node):
            for i in range(len(depvars)):
                indepvars = _ff.LSolveIndepVars(graph,node,i)
                indeplen = len(indepvars)
                rhsccs = _islice(graphout,idx,idx+indeplen)
                sol[depvars[i]] = self.Add(*(c * extvars[v] \
                                             for c,v in zip(rhsccs,indepvars)))
                idx+=indeplen

        else:
            indepvars = _ff.LSolveIndepVars(graph,node)
            if len(graphout) == len(depvars)*len(indepvars):
                pass
            elif len(graphout) == len(depvars)*(len(indepvars) + 1):
                indepvars = list(_chain(indepvars,(len(variables),)))
            else:
                raise ValueError("Graph output has invalid length")
            indeplen = len(indepvars)
            for i in range(len(depvars)):
                rhsccs = _islice(graphout,idx,idx+indeplen)
                sol[depvars[i]] = self.Add(*(c * extvars[v] \
                                             for c,v in zip(rhsccs,indepvars)))
                idx+=indeplen

        return sol


    def SparseSolve(self,
                    eqs, variables,
                    parameters=None,
                    indep_vars_only=False, mark_and_sweep=True,
                    sparse_output = True, only_non_homogeneous = False,
                    start_mod=0, max_primes=0, max_deg=0, dbginfo=0, polymethod=0,
                    n_threads=0, needed_vars=None):
        '''Returns a dictionary with the solution of the linear system.

        '''

        if parameters is None:
            parameters = set()
            varset = set(variables)
            for eq in eqs:
                parameters.update(self.VariablesIn(eq)-varset)
            parameters = list(parameters)

        with _ff.GraphContextWithInput(len(parameters)) as (g,inp):

            if len(parameters) > 0:
                (cols,ccs) = self.ToAnalyticSparseLSolve(parameters,
                                                         eqs, variables)
                ls = _ff.AlgAnalyticSparseLSolve(g, inp,
                                                 len(variables), cols, ccs)
            else:
                (cols,ccs) = self.ToNumericSparseLSolve(eqs, variables)
                ls = _ff.AlgNumericSparseLSolve(g, len(variables), cols, ccs,
                                                needed_vars=needed_vars)
            _ff.SetOutputNode(g,ls)
            _ff.LSolveOptimizeZeroVars(g,ls)

            if (not indep_vars_only) and sparse_output:
                _ff.LSolveSparseOutput(g,ls)
            if only_non_homogeneous:
                _ff.LSolveOnlyNonHomogeneous(g,ls)

            _ff.LearnEx(g,prime_no=start_mod)

            if indep_vars_only:
                indepvars = _ff.LSolveIndepVars(g,ls)
                indepvars = list(variables[idx] for idx in indepvars)
                return indepvars

            if mark_and_sweep:
                _ff.LSolveMarkAndSweepEqs(g, ls)
                _ff.LSolveDeleteUnneededEqs(g, ls)

            if len(parameters) > 0:
                rec = _ff.ReconstructFunction(g,
                                              start_mod=start_mod,
                                              max_primes=max_primes,
                                              max_deg=max_deg,
                                              dbginfo=dbginfo,
                                              polymethod=polymethod,
                                              n_threads=n_threads)
                if not type(rec) is _ff.RatFunList:
                    return rec
                rec = self.FromRatFunList(parameters, rec)
            else:
                rec = _ff.ReconstructNumeric(g,
                                             start_mod=start_mod,
                                             max_primes=max_primes,
                                             dbginfo=dbginfo,
                                             n_threads=n_threads)
                if not type(rec) is list:
                    return rec
                rec = self.FromRatNumList(rec)

            return self.LSolveDict(g,ls,variables,rec)


    def FromLaurentOutput(self,graph,node,ccs,expvar):
        '''Convert the output of a Laurent expansion node, i.e. the
        coefficients ccs, into the CAS.

        '''
        ret = _ff.LaurentOutput(graph,node,ccs)
        max_orders = _ff.LaurentMaxOrders(graph,node)
        AsSeries = self.AsSeries
        Add = self.Add
        return list(AsSeries(Add(*(cc*expvar**o for (o,cc) in el)),
                             expvar,0,maxo) \
                    for el,maxo in zip(ret,max_orders))

    def _PolyCCsExps2Sym(self, cciter, variables, exponents):
        Mul = self.Mul
        Add = self.Add
        Rational = self.Rational
        return Add(*(Mul(next(cciter),*(v**e for v,e in zip(variables,exps))) \
                     for exps in exponents))

    def _RatCCsExps2Sym(self, cciter, variables, num, den):
        return (self._PolyCCsExps2Sym(cciter,variables,num)/
                self._PolyCCsExps2Sym(cciter,variables,den))

    def FromCoeffsAndExponents(self,coefficients,variables,exponents):
        '''Return a list of rational functions from a list of
        coefficients, variables and exponent data.'''
        flat = _ff._Flattened
        RatCCsExps2Sym = self._RatCCsExps2Sym
        return list(RatCCsExps2Sym(flat(coefficients),variables,num,den) \
                    for num,den in exponents)
