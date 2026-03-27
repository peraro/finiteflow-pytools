# Packages

Python packages using FiniteFlow's Python interface:

- [ff2cas](ff2cas.py): tools for using FiniteFlow with Computer Algebra Systems (CAS) in Python.  It is not meant to be used directly, but it contains common code for implementing an interface that facilitates usage of FiniteFlow with a CAS.  See `ff2sym` and `ff2se` for specific implementations.

- [ff2sym](ff2sym.py): tools for using FiniteFlow with [SymPy](https://www.sympy.org).  SymPy is a CAS implemented in Python.  It is not very performant, but it has a significant number of features.  You may install SymPy with the command `pip3 install sympy --user` (or see its official documentation for more options).

- [ff2se](ff2se.py): tools for using FiniteFlow with the Python interface of [SymEngine](https://symengine.org).  SymEngine is a C++ CAS, which is much more efficient than SymPy despite having a similar interface, but it has more basic features.  You may install SymEngine's Python interface with the command `pip3 install symengine --user` (or see its official documentation for more options).

In order to use these packages, we recommend adding
```
export PYTHONPATH=$PYTHONPATH:/path/to/finiteflow-pytools/packages
```
to your `~/.bashrc`, `~/.zshrc` or `~/.bash_profile`.
