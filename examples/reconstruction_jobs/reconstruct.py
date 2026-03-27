from fflow import *
import glob

(nparsin,nparsout) = NParsFromDegreeFile("degrees.fflow")

# No need to use the actual graph, since we just want to attempt a
# reconstruction from saved evaluations
graph =  NewGraphDummy(nparsin, nparsout)

print(LoadDegrees(graph,"degrees.fflow"))
print(LoadEvaluations(graph, glob.glob("evaluations_*.fflow")))

ret = ReconstructFromCurrentEvaluations(graph)

if type(ret) is RatFunList:
    print("Function successfully reconstructed")
    zvars = ["z1", "z2", "z3", "z4"]
    print(ret.to_string(zvars))
else:
    print("Reconstruction failed with return value: {}".format(repr(ret)))
