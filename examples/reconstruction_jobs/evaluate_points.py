from fflow import *
import sys

# set inputs
filename="points_1.fflow"
start=0
npoints=226

# or read them from the command line
for i in range(len(sys.argv)):
    if sys.argv[i] == 'file':
        filename = sys.argv[i+1]
    elif sys.argv[i] == 'start':
        start = int(sys.argv[i+1])
    elif sys.argv[i] == 'npoints':
        npoints = int(sys.argv[i+1])

# load the graph
from define_graph import *

#print(LoadDegrees(graph,"degrees.fflow"))
print(EvaluatePointsInFile(graph,filename,start,npoints))

outfile = "evaluations_" + str(start) + "_" + str(npoints) + "_" + filename

if DumpEvaluations(graph,outfile) == Success():
    print("Saved evaluations in {}".format(outfile))
