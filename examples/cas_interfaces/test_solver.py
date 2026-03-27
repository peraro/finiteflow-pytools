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

from def_linear_system import *

sol = SparseSolve(eqs,unknowns)
check_solution(sol)
