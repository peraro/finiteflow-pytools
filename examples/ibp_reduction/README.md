IBP reduction and differential equations
========================================

In this example, we use FiniteFlow's Python API for Integration By
Parts (IBP) reduction.  In particular, we reconstruct the differential
equations for an integral family.  While the example is simple by
modern standards, it can be used as a template for more complex use
cases.

FiniteFlow relies on external programs for the generation of IBPs and
other identities.  For this example we have used the
[CALICO](https://github.com/fontana-g/calico) Mathematica package to
generate the input.  We provide the source of the scripts
`generate_ibps.wl` and `compute_derivatives.py` for generating the IBP
system and computing the derivatives of the Master Integrals (MIs)
respectively, both of which import `common_defs.m` which contains some
common definitions.  However, any other program capable of providing
such input can, in principle, be used.  In many applications,
generating such information is not a bottleneck, while the IBP
reduction and the functional reconstruction are.  The latter often
require the use of high-performance machines and computing clusters.
This example thus showcases how to achieve the IBP reduction and
functional reconstruction using the Python API of FiniteFlow.

We need to first generate information about the IBP system and the
integral family.  Here, we do that via the `generate_ibps.wl` script
```
  math -script generate_ibps.wl
```
which creates the folder `ibps/` containing such information, using
the CALICO package.

For this application, we assume we don't know in advance the list of
MIs.  Hence we reduce a selection of integrals (a.k.a. needed
integrals or target integrals) which is likely to be a superset of
those we need for computing differential equations.  The script
`generate_ibps.wl` saves this and other information in the JSON file
`ibps/system.json` and other auxiliary files.

Next, we run the script `ibps_and_des.py`
```
  python3 ibps_and_des.py
```
The script will look for files containing the list of MIs and their
derivatives.  If these are not present, it will assume that the MIs are
not known and therefore it will only solve the system numerically and
save a list of indexes corresponding to the MIs into a file.

Once the list of MIs is known, we compute their derivatives.  This is
done via the `compute_derivatives.wl` script
```
  math -script compute_derivatives.wl
```
which saves the information into the file `ibps/de.json`.  It also
updates the list of needed integrals in `ibps/system.json`.

Finally, we run again the script `ibps_and_des.py`, which now will
find the derivatives of the MIs and build a graph which reconstructs
the DEs of the MIs.  The information can thus be saved in various
formats.  We format it as analytic expressions written into strings
and then put into a list of matrices for the DEs, which we save into
`ibps/de.json`.  The file can be loaded into e.g. Mathematica or other
programs and languages for further use.

The example can be reset by deleting the `ibps/` directory.
