from fflow import *
import sys, glob

# set number of primes
nprimes = 1

# or read it from the command line
for i in range(len(sys.argv)):
    if sys.argv[i] == 'nprimes':
        nprimes = int(sys.argv[i+1])
        break

(nparsin,nparsout) = NParsFromDegreeFile("degrees.fflow")

# No need to use the actual graph to generate points, since we have
# the degrees already
graph =  NewGraphDummy(nparsin, nparsout)

print(LoadDegrees(graph,"degrees.fflow"))
print(LoadEvaluations(graph, glob.glob("evaluations_*.fflow")))

savefile = "points_" + str(nprimes) + ".fflow"

print(DumpSamplePoints(graph,savefile,max_primes=nprimes))

print(("Generated {} new points to be " +
       "evaluated in '{}'").format(NSamplePointsInFile(savefile), savefile))
