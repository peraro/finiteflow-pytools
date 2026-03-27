(* ::Package:: *)

<<CALICO`


Quiet[SetDirectory[NotebookDirectory[]]];


Get["common_defs.m"];


twist = CATBaikov[baikovpoly,d,2,3,zs];


ann = CATAnnihilator[twist,zs,2,1];


Length/@ann[[1]]


tmpids = CATInvMonIdsFromAnnihilators["fam",ann,zs];


maxrank=2;
maxdots=1;
seeds = Join@@(CATGenerateSeeds[#,{0,maxdots},{0,maxrank}]&/@nonzerosects);


equations = Join[simplifyeqs[ Flatten[tmpids@@#&/@seeds] ], Import["simple_maps.m"]];


intS = -Plus@@Select[#[[2]],#<0&]&;
intT = Plus@@sectorof[#]&;
intU = ((Plus@@Select[#[[2]],#>0&]) - intT[#])&;


integrals  = integralsin[equations];
needed = Select[integrals,intU[#]<=1 && intS[#]<=2 &];


(* serialize the identities in JSON format. *)
position = Association[{}];
Table[position[integrals[[ii]]]=ii-1;,{ii,Length[integrals]}];
Quiet[CreateDirectory["ibps"],CreateDirectory::eexist];
FFSparseEqsToJSON["ibps/sids_calicoids.json",{d,s12,s23},Thread[equations==0],integrals,CATUnorderedIntCases,position];
FFSparseSystemToJSON["ibps/system.json",Length[equations],integrals,{d,s12,s23},{"ibps/sids_calicoids.json"},"NeededVars"->needed];


(* export additional useful info *)
Export["ibps/sys_data.m",
{
"NEqs"->Length[equations],
"Integrals"->integrals,
"Needed"->needed
}
]
