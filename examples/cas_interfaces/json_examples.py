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


test_list = [42,1,2,3,42]
U32ListToJSON("mylist.fflow.json", test_list)
if not (test_list == U32ListFromJSON("mylist.fflow.json")):
    print("ERROR: JSON list test failed")
    exit(1)


from def_linear_system import *

half_neqs = len(eqs) // 2

eqs1 = eqs[:half_neqs]
eqs2 = eqs[half_neqs:]

def export(filename, equations):
    (nz_els, nz_ccs) = ToAnalyticSparseLSolve([x,y],equations,unknowns)
    SparseEqsToJSON(filename, nz_els, nz_ccs)

export("eqs1.fflow.json",eqs1)
export("eqs2.fflow.json",eqs2)
needevars = unknowns

SparseSystemToJSON("system.fflow.json", len(eqs), len(unknowns), 2,
                   PositionsInList(needevars,unknowns),
                   ["eqs1.fflow.json", "eqs2.fflow.json"])

with GraphContextWithInput(2) as (g,inp):

    ls = AlgJSONSparseLSolve(g,inp,"system.fflow.json")
    SetOutputNode(g,ls)
    LSolveOptimizeZeroVars(g,ls)
    LSolveSparseOutput(g,ls)

    Learn(g)
    LSolveMarkAndSweepEqs(g, ls)
    LSolveDeleteUnneededEqs(g, ls)

    rec = FromRatFunList([x,y], ReconstructFunction(g))
    sol = LSolveDict(g,ls,unknowns,rec)

    check_solution(sol)
