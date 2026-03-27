from fflow import *
import os
USE_SYMENGINE = os.getenv("FF_SYMENGINE") == "1"
if not USE_SYMENGINE:
    from sympy import *
    from ff2sym import *
    print("Using SymPy")
else:
    from symengine import *
    from ff2se import *
    print("Using SymEngine")

x,y,z,a,b,c = symbols("x y z a b c")

eqs = [
    Eq(x + 2*y + 3*z, a),
    Eq(7*x + 3*y - 4*z, b),
    Eq(x + 3*y + 13*z, (c+1)/(c-1))
]

print("Test SparseSolve")
sol = SparseSolve(eqs,[x,y,z])
print("-> ",sol)
print("-> ", list(Together((eq.args[0]-eq.args[1]).xreplace(sol))
                  .as_numer_denom()[0].expand() \
                  for eq in eqs))

params = [a,b,c]
xs = [x,y,z]
with GraphContextWithInput(len(params)) as (g,inp):
    print("Test LSolve in a graph")
    cols,ccs = ToAnalyticSparseLSolve(params,eqs,xs)
    ls = AlgAnalyticSparseLSolve(g,inp,len(xs),cols,ccs)
    SetOutputNode(g,ls)
    LSolveOptimizeZeroVars(g,ls)
    Learn(g)
    print("->  indep eqs: ",LSolveIndepEqs(g,ls))
    LSolveMarkAndSweepEqs(g,ls)
    LSolveDeleteUnneededEqs(g,ls)
    #print("zeroes : ",list(xs[idx] for idx in LSolveZeroVars(g,ls)))
    rec = ReconstructFunction(g)
    rec = FromRatFunList(params,rec)
    rec = LSolveDict(g,ls,xs,rec)
    print("-> ",rec)

with GraphContextWithInput(len(xs)) as (g,inp):
    print("Test Laurent in a graph")
    rf = AlgRatFunEval(g,inp,ToRatFunList(xs,[(x+y+z)/(x**2 - x*y - x*z)]))
    SetOutputNode(g,rf)
    with GraphContextWithInput(len(xs)-1) as (gl,inpl):
        laur = AlgLaurent(gl,inpl,g,2)
        SetOutputNode(gl,laur)
        Learn(gl)
        rec = FromRatFunList(xs[2:],ReconstructFunction(gl))
        print("-> ",FromLaurentOutput(gl,laur,rec,xs[0]))


with GraphContextWithInput(len(xs)) as (g,inp):
    print("Test subgraph reconstruct and EvalFromCoeffs")
    c = Function("c")
    func = [
        (c(0)*x + c(1)*y)/(c(2) + c(3)*x*y)
    ]
    conv = {
        c(0) : z**3,
        c(1) : z**2,
        c(2) : Integer(1),
        c(3) : z
    }
    cc_syms = list(c(i) for i in range(4))
    cc_vals = list(conv[c(i)] for i in range(4))
    print("-> ",func[0].xreplace(conv))

    ff_func0 = ToRatFunList(xs,[func[0].xreplace(conv)])
    f0 = AlgRatFunEval(g,inp,ff_func0)

    zonly = AlgSlice(g,inp,2,3)
    xyonly = AlgSlice(g,inp,0,2)
    ci = AlgRatFunEval(g,inp, ToRatFunList(xs,cc_vals))

    # from list of coeffs
    ff_func1 = ToRatFunListFromCoeffs(cc_syms, [x,y], func)
    f1 = AlgRatFunEvalFromCoeffs(g,ci,xyonly,ff_func1)

    # from a function
    ff_func2 = ToRatFunListFromCoeffs(lambda cc : int(cc.args[0]), [x,y], func)
    f2 = AlgRatFunEvalFromCoeffs(g,ci,xyonly,ff_func2)

    pt = [1234567890,567890123456,34567890123]
    SetOutputNode(g,f0)
    ev0 = EvaluateGraph(g,pt,0)
    SetOutputNode(g,f1)
    ev1 = EvaluateGraph(g,pt,0)
    SetOutputNode(g,f2)
    ev2 = EvaluateGraph(g,pt,0)

    print("-> ", ev0 == ev1 == ev2)

    with GraphContextWithInput(1) as (g2,in2):
        sub = AlgSubgraphRec(g2,in2,g,2)
        SetOutputNode(g2,sub)
        Learn(g2)

        rc = Function("rc")
        rcs = list(rc(i) for i in range(GraphNParsOut(g2)))

        out = FromCoeffsAndExponents(rcs,[x,y],SubgraphRecExponents(g2,sub))
        print("-> ",out)

        rec = FromRatFunList([z],ReconstructFunction(g2))
        rec = dict(zip(rcs,rec))
        print("-> ", out[0].xreplace(rec) == func[0].xreplace(conv))
