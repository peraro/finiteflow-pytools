'''Definition of a linear system for some tests.'''

import os
USE_SYMENGINE = os.getenv("FF_SYMENGINE") == "1"
if not USE_SYMENGINE:
    from sympy import *
    from ff2sym import *
else:
    from symengine import *
    from ff2se import *

x,y,z = symbols("x y z")
za = Function("za")
zb = Function("zb")
zc = Function("zc")
zd = Function("zd")
ze = Function("ze")
zf = Function("zf")

eqs = [Eq(y*zb(0) + y*zc(0) + y*zd(0) + y*ze(0), 0), Eq(y*zb(5) + y*zc(5) + y*zd(5) + y*ze(5), 0), Eq(y*zc(0) + y*zd(0) + y*zf(0), 0), Eq(za(0) + y*zb(0) + y*zb(4) + y*zc(0) + y*zc(4) + y*zc(5) + y*zd(4) + y*zd(5) + y*ze(4) + y*zf(5), 0), Eq(za(5) + y*zb(5) + y*zc(5), 0), Eq(y*zc(4) + y*zd(4) + y*zf(4), 0), Eq(za(4) + y*zb(4) + y*zc(4), 0), Eq(y*zb(0) + y*ze(0) + y*zf(0), 0), Eq(za(0) + y*zb(0) + y*zb(3) + y*zb(5) + y*zc(0) + y*zc(3) + y*zd(3) + y*ze(3) + y*ze(5) + y*zf(5), 0), Eq(za(5) + y*zb(5) + y*zc(5), 0), Eq(za(0) + y*zb(0) + y*zb(4) + y*zc(0) + y*zc(3) + y*zd(3) + y*ze(4) + y*zf(3) + y*zf(4), 0), Eq(za(3) + za(4) + za(5) + y*zb(3) + y*zb(4) + y*zb(5) + y*zc(3) + y*zc(4) + y*zc(5), 0), Eq(za(4) + y*zb(4) + y*zc(4), 0), Eq(y*zb(3) + y*ze(3) + y*zf(3), 0), Eq(za(3) + y*zb(3) + y*zc(3), 0), Eq(za(3) + y*zb(3) + y*zc(3), 0), Eq(y*zb(0) + y*ze(0) + y*zf(0), 0), Eq(za(0) + y*zb(2) + y*zb(5) + y*zc(2) + y*zd(0) + y*zd(2) + y*ze(0) + y*ze(2) + y*ze(5) + y*zf(5), 0), Eq(za(5) + y*zd(5) + y*ze(5), 0), Eq(za(0) + y*zb(0) + y*zb(4) + y*zc(2) + y*zd(0) + y*zd(2) + y*ze(4) + y*zf(0) + y*zf(2) + y*zf(4), 0), Eq(za(0) + za(2) + za(4) + za(5) + y*zb(2) + y*zb(5) + y*zc(2) + y*zd(4) + y*zd(5) + y*ze(4) + y*zf(5), 0), Eq(za(5), 0), Eq(za(4) + y*zb(4) + y*zd(4) + y*zf(4), 0), Eq(za(4), 0), Eq(y*zb(0) + y*zb(2) + y*zb(3) + y*ze(0) + y*ze(2) + y*ze(3) + y*zf(0) + y*zf(2) + y*zf(3), 0), Eq(za(0) + za(2) + za(3) + y*zb(2) + y*zb(5) + y*zc(2) + y*zd(3) + y*ze(3) + y*ze(5) + y*zf(5), 0), Eq(za(5), 0), Eq(za(0) + za(2) + za(3) + y*zb(2) + y*zb(3) + y*zb(4) + y*zc(2) + y*zd(3) + y*ze(4) + y*zf(3) + y*zf(4), 0), Eq(za(3) + za(4) + za(5), 0), Eq(za(4), 0), Eq(y*zb(3) + y*ze(3) + y*zf(3), 0), Eq(za(3), 0), Eq(za(3), 0), Eq(y*zb(2) + y*ze(2) + y*zf(2), 0), Eq(za(2) + y*zd(2) + y*ze(2), 0), Eq(za(2) + y*zb(2) + y*zd(2) + y*zf(2), 0), Eq(za(2), 0), Eq(y*zb(2) + y*ze(2) + y*zf(2), 0), Eq(za(2), 0), Eq(za(2), 0), Eq(y*zc(0) + y*zd(0) + y*zf(0), 0), Eq(za(0) + y*zb(1) + y*zc(1) + y*zc(5) + y*zd(0) + y*zd(1) + y*zd(5) + y*ze(0) + y*ze(1) + y*zf(5), 0), Eq(za(5) + y*zd(5) + y*ze(5), 0), Eq(y*zc(0) + y*zc(1) + y*zc(4) + y*zd(0) + y*zd(1) + y*zd(4) + y*zf(0) + y*zf(1) + y*zf(4), 0), Eq(za(0) + za(1) + za(4) + y*zb(1) + y*zc(1) + y*zc(5) + y*zd(4) + y*zd(5) + y*ze(4) + y*zf(5), 0), Eq(za(5), 0), Eq(y*zc(4) + y*zd(4) + y*zf(4), 0), Eq(za(4), 0), Eq(za(0) + y*zb(1) + y*zc(0) + y*zc(3) + y*zd(3) + y*ze(0) + y*ze(1) + y*zf(0) + y*zf(1) + y*zf(3), 0), Eq(za(0) + za(1) + za(3) + za(5) + y*zb(1) + y*zc(1) + y*zc(5) + y*zd(3) + y*ze(3) + y*ze(5) + y*zf(5), 0), Eq(za(5), 0), Eq(za(0) + za(1) + za(4) + y*zb(1) + y*zc(1) + y*zc(3) + y*zc(4) + y*zd(3) + y*ze(4) + y*zf(3) + y*zf(4), 0), Eq(za(3) + za(4) + za(5), 0), Eq(za(4), 0), Eq(za(3) + y*zc(3) + y*ze(3) + y*zf(3), 0), Eq(za(3), 0), Eq(za(3), 0), Eq(za(0) + y*zb(1) + y*zc(2) + y*zd(0) + y*zd(2) + y*ze(0) + y*ze(1) + y*zf(1) + y*zf(2), 0), Eq(za(1) + za(2) + za(5) + y*zd(1) + y*zd(2) + y*zd(5) + y*ze(1) + y*ze(2) + y*ze(5), 0), Eq(za(0) + za(1) + za(4) + y*zb(1) + y*zc(2) + y*zd(1) + y*zd(2) + y*zd(4) + y*ze(4) + y*zf(1) + y*zf(2), 0), Eq(za(1) + za(2) + za(5), 0), Eq(za(4), 0), Eq(za(0) + za(2) + za(3) + y*zb(1) + y*zc(2) + y*zd(3) + y*ze(1) + y*ze(2) + y*ze(3) + y*zf(1) + y*zf(2), 0), Eq(za(1) + za(2) + za(5), 0), Eq(za(1) + za(2) + za(3) + za(4), 0), Eq(za(3), 0), Eq(za(2) + y*zd(2) + y*ze(2), 0), Eq(za(2), 0), Eq(za(2), 0), Eq(y*zc(1) + y*zd(1) + y*zf(1), 0), Eq(za(1) + y*zd(1) + y*ze(1), 0), Eq(y*zc(1) + y*zd(1) + y*zf(1), 0), Eq(za(1), 0), Eq(za(1) + y*zc(1) + y*ze(1) + y*zf(1), 0), Eq(za(1), 0), Eq(za(1), 0), Eq(za(1) + y*zd(1) + y*ze(1), 0), Eq(za(1), 0), Eq(za(1), 0)]
unknowns = [za(1), za(2), za(3), za(4), za(5), za(0), zb(1), zb(2), zb(3), zb(4),zb(5), zb(0), zc(1), zc(2), zc(3), zc(4), zc(5), zc(0), zd(1), zd(2), zd(3), zd(4), zd(5), zd(0), ze(1), ze(2), ze(3), ze(4), ze(5),ze(0), zf(1), zf(2), zf(3), zf(4), zf(5), zf(0)]


def check_solution(sol):
    if USE_SYMENGINE:
        check = set(Eq((eq.args[0]-eq.args[1]).xreplace(sol).expand(),0) \
                    for eq in eqs)
    else:
        check = set(eq.xreplace(sol).expand() for eq in eqs)

    if not check ==  {Eq(0,0)}:
        print("ERROR: SparseSolve test failed")
        exit(1)
    else:
        print("Test passed!")
