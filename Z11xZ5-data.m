(* 
    Group data on Z11 \rtimes Z4
    -Alan Tran
*)

(* ::Package:: *)

BeginPackage["Z11xZ5data`"];

MapTo::usage = "MapTo[g] returns g in [\!\(\*SuperscriptBox[\(a\), \(j\)]\),\!\(\*SuperscriptBox[\(b\), \(k\)]\)] format";
Inv::usage = "Inv[g] returns inverse of g";
Mult::usage = "Mult[g1,g2,...gn] returns g1*g2*...gn";
Conj::usage = "Conj[g,h] returns Mult[g,h,Inv[g]]";

GetConjClass::usage = "%[g] returns conjugacy class of g";
GetCentralizer::usage = "%[g] returns centralizer of g";
GetRepOfConj::usage = "%[g] returns conj class rep associated to conj class of g";
GetCosetReps::usage = "%[g] returns coset reps of G/cent(g)";

ConjClassReps::usage = "Variable holdidng {1,3,6,2,4,7,11}";


(* ::Input::Initialization:: *)
Begin["Global`"]
Z11xZ5MultTable={{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55},{2,4,5,7,8,9,11,12,13,14,1,16,17,18,19,3,21,22,23,24,6,26,27,28,29,10,31,32,33,34,15,36,37,38,39,20,41,42,43,44,25,46,47,48,49,30,50,51,52,35,53,54,40,55,45},{3,19,6,28,24,10,51,33,29,15,26,54,38,34,20,31,7,43,39,25,36,12,48,44,30,41,17,52,49,35,46,22,4,2,40,50,27,8,5,45,53,32,13,9,1,55,37,18,14,11,42,23,16,47,21},{4,7,8,11,12,13,1,16,17,18,2,3,21,22,23,5,6,26,27,28,9,10,31,32,33,14,15,36,37,38,19,20,41,42,43,24,25,46,47,48,29,30,50,51,52,34,35,53,54,39,40,55,44,45,49},{5,23,9,32,28,14,53,37,33,19,10,55,42,38,24,15,11,47,43,29,20,16,51,48,34,25,21,54,52,39,30,26,7,4,44,35,31,12,8,49,40,36,17,13,2,45,41,22,18,1,46,27,3,50,6},{6,39,10,52,44,15,42,4,49,20,41,47,8,2,25,46,51,13,5,30,50,54,18,9,35,53,7,23,14,40,55,12,28,19,45,11,17,33,24,1,16,22,38,29,3,21,27,43,34,26,32,48,31,37,36},{7,11,12,1,16,17,2,3,21,22,4,5,6,26,27,8,9,10,31,32,13,14,15,36,37,18,19,20,41,42,23,24,25,46,47,28,29,30,50,51,33,34,35,53,54,38,39,40,55,43,44,45,48,49,52},{8,27,13,36,32,18,40,41,37,23,14,45,46,42,28,19,1,50,47,33,24,3,53,51,38,29,6,55,54,43,34,10,11,7,48,39,15,16,12,52,44,20,21,17,4,49,25,26,22,2,30,31,5,35,9},{9,43,14,54,48,19,46,7,52,24,25,50,12,4,29,30,53,17,8,34,35,55,22,13,39,40,11,27,18,44,45,16,32,23,49,1,21,37,28,2,3,26,42,33,5,6,31,47,38,10,36,51,15,41,20},{10,5,15,23,9,20,32,28,14,25,53,37,33,19,30,55,42,38,24,35,11,47,43,29,40,16,51,48,34,45,21,54,52,39,1,26,7,4,44,3,31,12,8,49,6,36,17,13,2,41,22,18,46,27,50},{11,1,16,2,3,21,4,5,6,26,7,8,9,10,31,12,13,14,15,36,17,18,19,20,41,22,23,24,25,46,27,28,29,30,50,32,33,34,35,53,37,38,39,40,55,42,43,44,45,47,48,49,51,52,54},{12,31,17,20,36,22,44,25,41,27,18,49,30,46,32,23,2,35,50,37,28,5,40,53,42,33,9,45,55,47,38,14,1,11,51,43,19,3,16,54,48,24,6,21,7,52,29,10,26,4,34,15,8,39,13},{13,47,18,55,51,23,30,11,54,28,29,35,16,7,33,34,40,21,12,38,39,45,26,17,43,44,1,31,22,48,49,3,36,27,52,2,6,41,32,4,5,10,46,37,8,9,15,50,42,14,20,53,19,25,24},{14,8,19,27,13,24,36,32,18,29,40,41,37,23,34,45,46,42,28,39,1,50,47,33,44,3,53,51,38,49,6,55,54,43,2,10,11,7,48,5,15,16,12,52,9,20,21,17,4,25,26,22,30,31,35},{15,24,20,48,29,25,22,52,34,30,16,27,4,39,35,21,32,8,44,40,26,37,13,49,45,31,42,18,2,1,36,47,23,5,3,41,51,28,9,6,46,54,33,14,10,50,7,38,19,53,12,43,55,17,11},{16,15,21,24,20,26,48,29,25,31,22,52,34,30,36,27,4,39,35,41,32,8,44,40,46,37,13,49,45,50,42,18,2,1,53,47,23,5,3,55,51,28,9,6,11,54,33,14,10,7,38,19,12,43,17},{17,50,22,45,53,27,34,1,55,32,33,39,3,11,37,38,44,6,16,42,43,49,10,21,47,48,2,15,26,51,52,5,20,31,54,4,9,25,36,7,8,14,30,41,12,13,19,35,46,18,24,40,23,29,28},{18,12,23,31,17,28,20,36,22,33,44,25,41,27,38,49,30,46,32,43,2,35,50,37,48,5,40,53,42,52,9,45,55,47,4,14,1,11,51,8,19,3,16,54,13,24,6,21,7,29,10,26,34,15,39},{19,28,24,51,33,29,26,54,38,34,3,31,7,43,39,6,36,12,48,44,10,41,17,52,49,15,46,22,4,2,20,50,27,8,5,25,53,32,13,9,30,55,37,18,14,35,11,42,23,40,16,47,45,21,1},{20,44,25,18,49,30,12,23,2,35,31,17,28,5,40,36,22,33,9,45,41,27,38,14,1,46,32,43,19,3,50,37,48,24,6,53,42,52,29,10,55,47,4,34,15,11,51,8,39,16,54,13,21,7,26},{21,35,26,49,40,31,38,2,45,36,37,43,5,1,41,42,48,9,3,46,47,52,14,6,50,51,4,19,10,53,54,8,24,15,55,7,13,29,20,11,12,18,34,25,16,17,23,39,30,22,28,44,27,33,32},{22,16,27,15,21,32,24,20,26,37,48,29,25,31,42,52,34,30,36,47,4,39,35,41,51,8,44,40,46,54,13,49,45,50,7,18,2,1,53,12,23,5,3,55,17,28,9,6,11,33,14,10,38,19,43},{23,32,28,53,37,33,10,55,42,38,5,15,11,47,43,9,20,16,51,48,14,25,21,54,52,19,30,26,7,4,24,35,31,12,8,29,40,36,17,13,34,45,41,22,18,39,1,46,27,44,3,50,49,6,2},{24,48,29,22,52,34,16,27,4,39,15,21,32,8,44,20,26,37,13,49,25,31,42,18,2,30,36,47,23,5,35,41,51,28,9,40,46,54,33,14,45,50,7,38,19,1,53,12,43,3,55,17,6,11,10},{25,9,30,43,14,35,54,48,19,40,46,7,52,24,45,50,12,4,29,1,53,17,8,34,3,55,22,13,39,6,11,27,18,44,10,16,32,23,49,15,21,37,28,2,20,26,42,33,5,31,47,38,36,51,41},{26,3,31,19,6,36,28,24,10,41,51,33,29,15,46,54,38,34,20,50,7,43,39,25,53,12,48,44,30,55,17,52,49,35,11,22,4,2,40,16,27,8,5,45,21,32,13,9,1,37,18,14,42,23,47},{27,36,32,40,41,37,14,45,46,42,8,19,1,50,47,13,24,3,53,51,18,29,6,55,54,23,34,10,11,7,28,39,15,16,12,33,44,20,21,17,38,49,25,26,22,43,2,30,31,48,5,35,52,9,4},{28,51,33,26,54,38,3,31,7,43,19,6,36,12,48,24,10,41,17,52,29,15,46,22,4,34,20,50,27,8,39,25,53,32,13,44,30,55,37,18,49,35,11,42,23,2,40,16,47,5,45,21,9,1,14},{29,13,34,47,18,39,55,51,23,44,30,11,54,28,49,35,16,7,33,2,40,21,12,38,5,45,26,17,43,9,1,31,22,48,14,3,36,27,52,19,6,41,32,4,24,10,46,37,8,15,50,42,20,53,25},{30,29,35,13,34,40,47,18,39,45,55,51,23,44,1,11,54,28,49,3,16,7,33,2,6,21,12,38,5,10,26,17,43,9,15,31,22,48,14,20,36,27,52,19,25,41,32,4,24,46,37,8,50,42,53},{31,20,36,44,25,41,18,49,30,46,12,23,2,35,50,17,28,5,40,53,22,33,9,45,55,27,38,14,1,11,32,43,19,3,16,37,48,24,6,21,42,52,29,10,26,47,4,34,15,51,8,39,54,13,7},{32,53,37,10,55,42,5,15,11,47,23,9,20,16,51,28,14,25,21,54,33,19,30,26,7,38,24,35,31,12,43,29,40,36,17,48,34,45,41,22,52,39,1,46,27,4,44,3,50,8,49,6,13,2,18},{33,17,38,50,22,43,45,53,27,48,34,1,55,32,52,39,3,11,37,4,44,6,16,42,8,49,10,21,47,13,2,15,26,51,18,5,20,31,54,23,9,25,36,7,28,14,30,41,12,19,35,46,24,40,29},{34,33,39,17,38,44,50,22,43,49,45,53,27,48,2,1,55,32,52,5,3,11,37,4,9,6,16,42,8,14,10,21,47,13,19,15,26,51,18,24,20,31,54,23,29,25,36,7,28,30,41,12,35,46,40},{35,49,40,38,2,45,37,43,5,1,21,42,48,9,3,26,47,52,14,6,31,51,4,19,10,36,54,8,24,15,41,7,13,29,20,46,12,18,34,25,50,17,23,39,30,53,22,28,44,55,27,33,11,32,16},{36,40,41,14,45,46,8,19,1,50,27,13,24,3,53,32,18,29,6,55,37,23,34,10,11,42,28,39,15,16,47,33,44,20,21,51,38,49,25,26,54,43,2,30,31,7,48,5,35,12,52,9,17,4,22},{37,21,42,35,26,47,49,40,31,51,38,2,45,36,54,43,5,1,41,7,48,9,3,46,12,52,14,6,50,17,4,19,10,53,22,8,24,15,55,27,13,29,20,11,32,18,34,25,16,23,39,30,28,44,33},{38,37,43,21,42,48,35,26,47,52,49,40,31,51,4,2,45,36,54,8,5,1,41,7,13,9,3,46,12,18,14,6,50,17,23,19,10,53,22,28,24,15,55,27,33,29,20,11,32,34,25,16,39,30,44},{39,52,44,42,4,49,41,47,8,2,6,46,51,13,5,10,50,54,18,9,15,53,7,23,14,20,55,12,28,19,25,11,17,33,24,30,16,22,38,29,35,21,27,43,34,40,26,32,48,45,31,37,1,36,3},{40,14,45,8,19,1,27,13,24,3,36,32,18,29,6,41,37,23,34,10,46,42,28,39,15,50,47,33,44,20,53,51,38,49,25,55,54,43,2,30,11,7,48,5,35,16,12,52,9,21,17,4,26,22,31},{41,6,46,39,10,50,52,44,15,53,42,4,49,20,55,47,8,2,25,11,51,13,5,30,16,54,18,9,35,21,7,23,14,40,26,12,28,19,45,31,17,33,24,1,36,22,38,29,3,27,43,34,32,48,37},{42,41,47,6,46,51,39,10,50,54,52,44,15,53,7,4,49,20,55,12,8,2,25,11,17,13,5,30,16,22,18,9,35,21,27,23,14,40,26,32,28,19,45,31,37,33,24,1,36,38,29,3,43,34,48},{43,54,48,46,7,52,25,50,12,4,9,30,53,17,8,14,35,55,22,13,19,40,11,27,18,24,45,16,32,23,29,1,21,37,28,34,3,26,42,33,39,6,31,47,38,44,10,36,51,49,15,41,2,20,5},{44,18,49,12,23,2,31,17,28,5,20,36,22,33,9,25,41,27,38,14,30,46,32,43,19,35,50,37,48,24,40,53,42,52,29,45,55,47,4,34,1,11,51,8,39,3,16,54,13,6,21,7,10,26,15},{45,34,1,33,39,3,17,38,44,6,50,22,43,49,10,53,27,48,2,15,55,32,52,5,20,11,37,4,9,25,16,42,8,14,30,21,47,13,19,35,26,51,18,24,40,31,54,23,29,36,7,28,41,12,46},{46,25,50,9,30,53,43,14,35,55,54,48,19,40,11,7,52,24,45,16,12,4,29,1,21,17,8,34,3,26,22,13,39,6,31,27,18,44,10,36,32,23,49,15,41,37,28,2,20,42,33,5,47,38,51},{47,55,51,30,11,54,29,35,16,7,13,34,40,21,12,18,39,45,26,17,23,44,1,31,22,28,49,3,36,27,33,2,6,41,32,38,5,10,46,37,43,9,15,50,42,48,14,20,53,52,19,25,4,24,8},{48,22,52,16,27,4,15,21,32,8,24,20,26,37,13,29,25,31,42,18,34,30,36,47,23,39,35,41,51,28,44,40,46,54,33,49,45,50,7,38,2,1,53,12,43,5,3,55,17,9,6,11,14,10,19},{49,38,2,37,43,5,21,42,48,9,35,26,47,52,14,40,31,51,4,19,45,36,54,8,24,1,41,7,13,29,3,46,12,18,34,6,50,17,23,39,10,53,22,28,44,15,55,27,33,20,11,32,25,16,30},{50,45,53,34,1,55,33,39,3,11,17,38,44,6,16,22,43,49,10,21,27,48,2,15,26,32,52,5,20,31,37,4,9,25,36,42,8,14,30,41,47,13,19,35,46,51,18,24,40,54,23,29,7,28,12},{51,26,54,3,31,7,19,6,36,12,28,24,10,41,17,33,29,15,46,22,38,34,20,50,27,43,39,25,53,32,48,44,30,55,37,52,49,35,11,42,4,2,40,16,47,8,5,45,21,13,9,1,18,14,23},{52,42,4,41,47,8,6,46,51,13,39,10,50,54,18,44,15,53,7,23,49,20,55,12,28,2,25,11,17,33,5,30,16,22,38,9,35,21,27,43,14,40,26,32,48,19,45,31,37,24,1,36,29,3,34},{53,10,55,5,15,11,23,9,20,16,32,28,14,25,21,37,33,19,30,26,42,38,24,35,31,47,43,29,40,36,51,48,34,45,41,54,52,39,1,46,7,4,44,3,50,12,8,49,6,17,13,2,22,18,27},{54,46,7,25,50,12,9,30,53,17,43,14,35,55,22,48,19,40,11,27,52,24,45,16,32,4,29,1,21,37,8,34,3,26,42,13,39,6,31,47,18,44,10,36,51,23,49,15,41,28,2,20,33,5,38},{55,30,11,29,35,16,13,34,40,21,47,18,39,45,26,51,23,44,1,31,54,28,49,3,36,7,33,2,6,41,12,38,5,10,46,17,43,9,15,50,22,48,14,20,53,27,52,19,25,32,4,24,37,8,42}};

Group=Z11xZ5MultTable;
GroupInGenerators={{0,0},{0,1},{1,0},{0,2},{1,1},{2,0},{0,3},{1,2},{2,1},{3,0},{0,4},{1,3},{2,2},{3,1},{4,0},{1,4},{2,3},{3,2},{4,1},{5,0},{2,4},{3,3},{4,2},{5,1},{6,0},{3,4},{4,3},{5,2},{6,1},{7,0},{4,4},{5,3},{6,2},{7,1},{8,0},{5,4},{6,3},{7,2},{8,1},{9,0},{6,4},{7,3},{8,2},{9,1},{10,0},{7,4},{8,3},{9,2},{10,1},{8,4},{9,3},{10,2},{9,4},{10,3},{10,4}};
MapTo[g_]:=GroupInGenerators[[g]];


(* ::Input::Initialization:: *)
Inv[g_]:=Position[Group[[g]],1][[1,1]];
Mul[g1_,g2_]:=Group[[g1,g2]];
Mult[g1_,g2__]:=Fold[Mul,g1,{g2}];
Conj[h_,g_]:=Mult[h,g,Inv[h]];


(* ::Input::Initialization:: *)
Elements=Z11xZ5MultTable[[1]];
order=Group //Length;

GetConjClass[g_]:=Sort[DeleteDuplicates[Table[Mult[Elements[[i]],g,Inv[Elements[[i]]]],{i,1,order}]]];

GetCentralizer[t_]:=Module[{left,right},
    left=Table[Mul[g,t],{g,1,order}];
    right=Table[Mul[t,g],{g,1,order}];
    Return[Position[left-right, 0] //Flatten];
];

GetRepOfConj[g_]:=Sort[GetConjClass[g]][[1]];

GetCosetReps[k_]:=Module[{cent,cosets},
	cent=GetCentralizer[k];
	cosets=Table[Mult[i,cent]//Sort,{i,1,55}]//DeleteDuplicates;
	cosets[[All,1]]
];
ConjClassReps={1,3,6,2,4,7,11};


End[]
EndPackage[]



