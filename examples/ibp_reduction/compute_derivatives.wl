(* ::Package:: *)

<<CALICO`


Quiet[SetDirectory[NotebookDirectory[]]];


(* ::Text:: *)
(*Load common definition for the integral family.*)


Get["common_defs.m"];


(* ::Text:: *)
(*Load information about the list of integrals in the IBP system.*)


sysdata = Import["ibps/sys_data.m"];


integrals = "Integrals"/.sysdata;


(* ::Text:: *)
(*Get the list of master integrals (MIs), obtained by solving the system in Python and exported as a list of indexes in JSON format.*)


mastersidxs = Import["ibps/mis.json"]+1;
masters = integrals[[mastersidxs]];


masters


masters//Length


(* ::Text:: *)
(*Now we compute the derivates of the MIs.*)


twist = CATMultiBaikov[{extpoly,baikovpoly},{extexponent,baikovexponent},zs];


diffop = CATDiffOperator[twist,zs,{s12,s23},2,1];


derivs = CATInvMonIdsFromDiffOperators["fam",diffop,zs];


derivatives = Transpose[ simplifyeqs[ Collect[ derivs@@#[[2]]&/@masters, _CATInt, Factor] ] ];


(* ::Text:: *)
(*The new list of needed integrals contains those appearing in the derivatives.  We update "ibps/system.json" with the new list.*)


newneededints = integralsin[derivatives];


(* check *)
If[!SubsetQ["Needed"/.sysdata,newneededints],
  Print["ERROR: new set of needed integrals must be "
        <>"a subset of the previous one!"];
  Quit[1];
];


FFSparseSystemToJSON["ibps/system.json","NEqs"/.sysdata,integrals,{d,s12,s23},{"ibps/sids_calicoids.json"},"NeededVars"->newneededints];


(* ::Text:: *)
(*Now we format the IBP data and export it to JSON format.  One file contains the list of (unique) rational functions multiplying the loop integrals in the derivatives.  The other one contains, for each external invariant and each master integral, a list of pairs {intidx,ccidx} with the indexes of the integral into the list of the unknowns and the function the multiplies it into the previous list.*)


derivdata = Map[FFSparseRowRules[#,integrals]&,derivatives,{2}];
uniquederivcc = Flatten[derivdata[[;;,;;,;;,2]]]//DeleteDuplicates;
uniquetopos = Association[Thread[uniquederivcc->Range[0,Length[uniquederivcc]-1]]];
derivdata = Map[{#[[1]]-1,uniquetopos[#[[2]]]}&,derivdata,{3}];


FFRatFunToJSON["ibps/derivatives_coefficients.json",{d,s12,s23},uniquederivcc]


Export["ibps/derivatives.json",derivdata,"RawJSON","Compact"->True]
