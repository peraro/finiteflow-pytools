'''This example illustrates the usage of FiniteFlow for reduction to
master integrals and reconstructing differential equations.

The system of IBP identities and the derivatives of the master
integrals need to be provided externally (we provide Mathematica files
which generate them using the FFIntRed package).  Check out the
README.md file in this directory for more information.

'''

from fflow import *
from functools import reduce
import json # to import/export information


# We reconstruct w.r.t. these variables
variables = "d,s12,s23".split(",")

# NOTE: s23 is set to 1 in the reconstructed result but we still
# compute the DEs w.r.t. it
invariants = "s12,s23".split(",")

graph,inp = NewGraphWithInput(len(variables))

ibps = AlgJSONSparseLSolve(graph,inp,"ibps/system.json")
LSolveOnlyHomogeneous(graph,ibps)
LSolveSparseOutput(graph,ibps)
SetOutputNode(graph,ibps)

Learn(graph)
nonmis = LSolveDepVars(graph,ibps)

# We use a sparse output, hence we get a separate list of MIs per
# unknown
mis = [LSolveIndepVars(graph,ibps,i) for i in range(len(nonmis))]

# This takes the union of all MIs
union_mis = sorted(list(reduce(lambda x,y : x.union(y), (set(x) for x in mis))))

try:

    # see if we know the MIs already
    with open("ibps/mis.json") as f:
        known_mis = json.load(f)

except FileNotFoundError:

    # If we don't know the MIs, we save them and quit for now.
    # Then we will compute derivatives in a separate program.
    with open("ibps/mis.json","w") as f:
        json.dump(union_mis,f)
    print("The list of MIs has been saved.")
    exit(0)


if known_mis != union_mis:
    print("New list of MIs differs from the saved one.")
    print("{} vs {}".format(union_mis,known_mis))
    exit(1)


# Optimize system
LSolveMarkAndSweepEqs(graph,ibps)
LSolveDeleteUnneededEqs(graph,ibps)
# and print some info
indep_eqs = LSolveNIndepEqs(graph,ibps)
tot_eqs = LSolveNEqsNVars(graph,ibps)[0]
print("System trimmed down to {} equations out of {}".format(indep_eqs,tot_eqs))


# MIs could appear in the derivatives => augment the system with the
# reduction of MIs to MIs
misred = AlgRatNumEval(graph,["1" for _ in range(len(union_mis))])
ibpsfull = AlgChain(graph,[ibps,misred])


# indexes of integrals in the list of nonmis
nonmisidx = dict(zip(nonmis + union_mis, range(len(nonmis)+len(union_mis))))

# indexes of integrals in the list of mis
misidx = dict(zip(union_mis,range(len(mis))))

# at this stage we expect the derivatives to have been computed
derivs = AlgJSONRatFunEval(graph,inp,"ibps/derivatives_coefficients.json")
with open("ibps/derivatives.json") as f:
    ccidxs = json.load(f)
ccidxs = [[sorted(e,key=lambda x:nonmisidx[x[0]]) for e in deriv] \
          for deriv in ccidxs]

def take_pattern():
    for deriv in ccidxs:
        for e in deriv:
            for cc in e:
                yield (0,cc[1])
unreduced = AlgTake(graph,[derivs],list(take_pattern()))


def mat1_rows():
    for deriv in ccidxs:
        for e in deriv:
            yield [nonmisidx[cc[0]] for cc in e]

def mat2_rows():
    for mislist in mis:
        yield [misidx[el] for el in mislist]
    for el in union_mis:
        yield [misidx[el]]

de = AlgSparseMatMul(graph,unreduced,ibpsfull,
                     len(invariants)*len(union_mis),
                     len(nonmis)+len(union_mis), len(union_mis),
                     list(mat1_rows()), list(mat2_rows()))
SetOutputNode(graph,de)


# To simplify the reconstruction we restrict it to a subset of unique
# elements
print("Taking unique elements of output")
unique,indexes = TakeUnique(graph,de)
SetOutputNode(graph,unique)

print("Reconstruct analytic result")
rec = ReconstructFunction(graph,dbginfo=1)
if not type(rec) is RatFunList:
    print("Reconstruction failed with status: {}".format(repr(rec)))
    exit(1)


# There are many different ways to save the result.  The most portable
# one would be building a list of strings containing the analytic
# expressions.

print("Format and save the output")

unique_analytic = rec.to_string(variables)
full_analytic = list(unique_analytic[i] for i in indexes)

# Format it as an array of len(invariants) matrices
n_mis = len(union_mis)
def de_matrix(input_array):
    return [
        input_array[i*n_mis:(i+1)*n_mis] for i in range(n_mis)
    ]

n_invs = len(invariants)
full_out = [
    de_matrix(full_analytic[i*n_mis*n_mis:(i+1)*n_mis*n_mis]) \
    for i in range(n_invs)
]
with open("ibps/de.json","w") as f:
    json.dump(full_out,f)
