from fflow import *
from define_graph import *

# NOTE: AllDegrees returns just the total degrees but it computes and
# stores internally all the degrees w.r.t. all the
# variables. DumpDegrees saves all of them.
print(AllDegrees(graph))
print(DumpDegrees(graph,"degrees.fflow"))
