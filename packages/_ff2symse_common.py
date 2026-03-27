'''Common functions used in FF2Sym and FF2SE
'''

import fflow as _ff


def Together(mod,expr):
    '''Put expression under common denominator.'''
    num,den = expr.as_numer_denom()
    return num/den

def FastNumerDenom(mod,expr):
    '''Return numerator and denominator, assuming the input is a
    product of polynomials raised to integer powers.'''
    pows = expr.as_powers_dict()
    num = set()
    den = set()
    for k,v in pows.items():
        if v < 0:
            den.add(k**-v)
        else:
            num.add(k**v)
    return mod.Mul(*num), mod.Mul(*den)

def AlgToExprEval(mod, graph, in_node, expressions, variables):
    dummies = tuple(mod.Dummy("t"+str(i)) for i in range(len(variables)))
    todummies = dict(zip(variables,dummies))
    strfuns = [str(fun.xreplace(todummies)) for fun in expressions]
    return _ff.AlgRatExprEval(graph, in_node, strfuns,
                              variable_prefix='_t',
                              n_vars=len(variables))
