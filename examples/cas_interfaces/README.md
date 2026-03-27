# Tests with SymPy and SymEngine

A few tests and examples using FiniteFlow with [SymPy](https://www.sympy.org) and the Python interface of [SymEngine](https://symengine.org) via `ff2sym` or `ff2se` (see their description [here](../../packages/README.md)).

By default the examples use SymPy.  To select SymEngine instead of SymPy, set the `FF_SYMENGINE` environment variable to `1`.  For instance
```
python3 basic_examples.py # will use SymPy
FF_SYMENGINE=1 python3 basic_examples.py # will use SymEngine
```
