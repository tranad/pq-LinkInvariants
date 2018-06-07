(*
    Compute link invariants for Z(Vec_G)
    -Alan Tran
*)

(* ::Package:: *)

BeginPackage["DZ11xZ5`"];

ComputeLinkInvariant::usage = "Main function %[u,{{braid},numStrand}]"
<<"Z11xZ5-data.m"
Begin["Global`"]
DistributeDefinitions[MapTo,Inv,Mult,Conj,GetConjClass,GetCentralizer,GetRepOfConj,GetCosetReps,ConjClassReps]

(* The cocyle *)
\[Omega][gg_,hh_,kk_,u_]:=Module[{g,h,k},
    g=MapTo[gg][[2]];
    h=MapTo[hh][[2]];
    k=MapTo[kk][[2]];
    Return[Exp[(I 2 \[Pi])/25 u k(Mod[g,5]+Mod[h,5]-Mod[g+h,5])]]
];

(* Appears when computing the action of a braid *)
\[Theta][g_,x_,y_,u_]:=\[Omega][g,x,y,u] Conjugate[\[Omega][x,Conj[Inv[x],g],y,u]] \[Omega][x,y,Conj[Inv[Mult[x,y]],g],u];

(* Number of projective irreps corresponding to centralizer of g *)
NumProjIrreps[g_]:=Module[{ll=Length[GetCentralizer[g]]},If[ll==55,7,ll]];

(* Character table for the group G=(11,5,4) *)
CharTableG[g_,n_]:=Module[{l,cc,s},
    cc=GetConjClass[g];
    l=GetConjClass[g]//Length;
    s=Sum[Exp[(I 2 \[Pi])/11 MapTo[cc[[iiii]]][[1]]],{iiii,1,l}];
    Which[l==1, Which[n<=4,1,n==5,5,n==6,5],
        l==11, Which[n<=4,Exp[(I 2 \[Pi])/5 n MapTo[g][[2]]],n==5,0,n==6,0],
        l==5, Which[n<=4,1,n==5,s,n==6,Conjugate[s]]
    ]
];

(* Projective characters for centralizer of t, evaluated on g with the n'th irrep for the u'th cocycle *)
\[CapitalPi][t_,g_,n_,u_]:=Module[{CENTL},
    CENTL=GetCentralizer[t]//Length;
    Which[CENTL==55,CharTableG[g,n],
        CENTL==11,Return[Exp[(I 2 \[Pi])/11 n MapTo[Mult[1,g]][[1]] ],Module],
        CENTL==5,Return[Exp[(I 2 \[Pi])/5 n MapTo[Mult[1,g]][[2]]]Exp[(I 2 \[Pi])/5^2 u MapTo[g][[2]]MapTo[t][[2]] ], Module]
    ]
];

(* Find the permutation of the under strand during a braid *)
FindRpY[t2_,g_,r2_]:=Module[{R,Y,RxY,coords,rp,y,numSol},
    R=GetCosetReps[t2]; 
    Y=GetCentralizer[t2];
    RxY=Table[Mult[R[[iii]],Y[[jjj]]],{iii,1, Length[R]},{jjj,1,Length[Y]}];
    coords=Position[RxY,Mult[g,r2]]//Flatten;
    rp=R[[coords[[1]]]];
    y=Y[[coords[[2]]]];
    Return[{rp,y},Module]
];

(* (inverse of the) associator function *) 
Associator[index_,vecG_,u_]:= Which[
		index == 1, Return[1],
		index==2, Return[\[Omega][vecG[[1]],vecG[[2]],vecG[[3]],u]],
		index>=3, Return[ \[Omega][vecG[[index-1]],vecG[[index]],vecG[[index+1]],u]*
            \[Omega][vecG[[index-2]],Mult[vecG[[index-1]],vecG[[index]]],vecG[[index+1]],u]*
            Associator[index-1, vecG, u]]
];
	
(* Compute the action for a braid. 1 is the over strand and 2 is the understrand *)
Action[sign_,t1_,t2_,r1_,r2_,n1_,n2_,u_]:=Module[{g, h, hp, rpy, multfactor},
    g=Mult[r1,t1,Inv[r1]];
    If[sign>0,
        rpy = FindRpY[t2,g,r2];
        h = Conj[rpy[[1]], t2];
        multfactor = \[Theta][h, g, r2, u]*Conjugate[\[Theta][h, rpy[[1]], rpy[[2]], u]]
        ,
        rpy = FindRpY[t2,Inv[g],r2];
        h = Conj[rpy[[1]], t2];
        multfactor = Conjugate[\[Theta][Conj[g,h], g, Inv[g], u]]* \[Theta][h, Inv[g], r2, u] * Conjugate[\[Theta][h, rpy[[1]], rpy[[2]], u]]
    ];
    Return[{multfactor, rpy[[1]], rpy[[2]]}, Module]
];

(* Finds the strands belonging to the same link component for a given ending configuration of strands *)
TracePath[configE_]:=ConnectedComponents[Graph[MapThread[DirectedEdge,{Range[1,Length[configE]],configE}]]];


(* ::Input::Initialization:: *)
(* 
    Computes a braid sequence with input strands labeled by [t, r, n] where t is a conj class rep
    r is an element of the coset representative G/Z(t) and n is a labeling for a projective irrep.
    At each step of the braidword the action is computed. 
*)
ComputeBraidSequence[braidword_,vecTin_,vecRin_,vecPin_,u_]:=
Module[{actseq,vecTout,vecRout,vecPout,tAbove,tBelow,rAbove,rBelow,nAbove,nBelow,
    index,sign,act,pathStart,pathEnd,yList,yClosedList,yAct,revBraid,yindex,componentsIn,componentsOut},
    (* initialize starting configuration of strands (strands evolve bottom to top) *)
	pathStart=Range[1,Length[vecTin]];
	pathEnd=pathStart;
	vecTout=vecTin;
	vecRout=vecRin;
	vecPout=vecPin;
	actseq=1;
	yList=Table[{1},{i,1,Length[vecPin]}];
    (* Loop through the braid word *)
	For[iter=1,iter<=Length[braidword],iter++,
		sign=Sign[braidword[[iter]]];
		index=Abs[braidword[[iter]]];
		If[index==0,Continue[]];
        actseq*=Conjugate[Associator[index, Table[Conj[vecRout[[i]],vecTout[[i]]],{i,1,Length[vecTin]}],u]];
        tAbove=vecTout[[index+KroneckerDelta[sign,-1]]];
        tBelow=vecTout[[index+KroneckerDelta[sign,1]]];
        rAbove=vecRout[[index+KroneckerDelta[sign,-1]]];
        rBelow=vecRout[[index+KroneckerDelta[sign,1]]];
        nAbove=vecPout[[index+KroneckerDelta[sign,-1]]];
        nBelow=vecPout[[index+KroneckerDelta[sign,1]]];
        act=Action[sign,tAbove,tBelow,rAbove,rBelow,nAbove,nBelow,u];
        actseq*=act[[1]];
        If[actseq==0,Return[0,Module]];
        yList[[index+KroneckerDelta[sign,1]]]=Append[yList[[index+KroneckerDelta[sign,1]]],act[[3]]];
        vecRout=ReplacePart[vecRout,(index+KroneckerDelta[sign,1])->act[[2]]];
        $HistoryLength=0;
        vecTout[[{index,index+1}]]=vecTout[[{index+1,index}]];
        vecRout[[{index,index+1}]]=vecRout[[{index+1,index}]];
        vecPout[[{index,index+1}]]=vecPout[[{index+1,index}]];
        pathEnd[[{index,index+1}]]=pathEnd[[{index+1,index}]];
        yList[[{index,index+1}]]=yList[[{index+1,index}]];
        (*after swapping, go back to original association pattern*)
        actseq*=Associator[index, Table[Conj[vecRout[[i]],vecTout[[i]]],{i,1,Length[vecTin]}],u];
    ];
    If[vecTout==vecTin&&vecRout==vecRin&&vecPout==vecPin,
        componentsIn=TracePath[pathEnd];
        componentsOut=Table[Table[Position[pathEnd,componentsIn[[i,j]]],{j,1,Length[componentsIn[[i]]]}]//Flatten,{i,1,Length[componentsIn]}];
        yClosedList=Table[If[Length[componentsOut[[i]]]>1,Join[yList[[componentsOut[[i]]]]]//Flatten,yList[[componentsOut[[i,1]]]]],{i,1,Length[componentsOut]}];
        yAct=Product[\[CapitalPi][vecTout[[componentsOut[[i,1]]]],Mult[1,1,Sequence@@yClosedList[[i]]],vecPout[[componentsOut[[i,1]]]],u] Product[\[Theta][vecTout[[componentsOut[[i,1]]]],yClosedList[[i,j]],Mult[1,1,Sequence@@yClosedList[[i]][[j+1;;]]],u],{j,1,Length[yClosedList[[i]]]-1}],{i,1,Length[componentsOut[[All,1]]]}];
        Return[actseq*yAct],
        Return[0]
    ]
];

(* generate labels for the link *)
ConjClassReps[[{1,2,3,4,5,6,7}]]
MakeLoopedTuples[braid_,numStrands_,u_,fluxList_:{1,2,3,4,5,6,7},chargeStart_:0,chargeOffset_:0]:=
Module[{tups,sTP,strandConj,strandIrrep,temp,\[Sigma],t,tp,pathEnd,comps,numComps,TP,loopTuples},
	pathEnd=Range[1,numStrands];
	For[iter=1,iter<=Length[braid],iter++,
		\[Sigma]=Abs[braid[[iter]]];
		$HistoryLength=0;
		pathEnd[[{\[Sigma],\[Sigma]+1}]]=pathEnd[[{\[Sigma]+1,\[Sigma]}]];
	];
	comps=TracePath[pathEnd];
	numComps=comps//Length;
	Print["Number of link components/number of strands: ",numComps,"/",numStrands];
	t=Tuples[ConjClassReps[[fluxList]],numComps];
	tp=ParallelTable[{t[[i]],Tuples[Table[Range[chargeStart,NumProjIrreps[t[[i,iter]]]-1-chargeOffset],{iter,1,numComps}]]},{i,1,Length[t]}];
	TP=Flatten[ParallelTable[{tp[[i,1]],tp[[i,2,j]]},{i,1,Length[tp]},{j,1,Length[tp[[i,2]]]}],1];
	sTP=SortBy[TP,Sum[(Riffle[Table[Position[ConjClassReps,#[[1,i]]],{i,1,Length[#[[1]]]}],#[[2]]][[i]]+1)*100^(Length[#[[1]]]*2-i),{i,1,Length[#[[1]]]*2}]&];
	strandConj=Table[sTP[[All,1]][[All,i]],{i,1,numComps}];
	strandIrrep=Table[sTP[[All,2]][[All,i]],{i,1,numComps}];
	temp=Range[1,2*numStrands];
	For[s=1,s<=numComps,s++,
		For[ss=1,ss<=Length[comps[[s]]],ss++,
			temp=ReplacePart[temp,comps[[s,ss]]->strandConj[[s]]];
			temp=ReplacePart[temp,numStrands+comps[[s,ss]]->strandIrrep[[s]]]
			]
		];
	temp=Transpose[temp];
	temp=Partition[#,numStrands]&/@temp;
	tups=Prepend[{#},{braid,numStrands,u}]&/@temp;
	{numComps,Flatten[#,1]&/@tups}
	];

(* do the trace over the braidword *)
TraceActSeq[tuple_]:=Module[{braidword,nStrands,vecTin,vecPin,u,vecRin,cosetR,intC},
	braidword=tuple[[1]];
	vecTin=tuple[[4]];
	vecPin=tuple[[5]];
	nStrands=tuple[[2]];
	u=tuple[[3]];
	cosetR=GetCosetReps[#]&/@vecTin;
	vecRin=Tuples[cosetR];
	Return[Sum[ComputeBraidSequence[braidword,vecTin,vecRin[[i]],vecPin,u],{i,1,Length[vecRin]}]];
];


(* Main function to compute a link invariant *)
ComputeLinkInvariant[uu_,braidAndNumStrands_,fluxList_:{1,2,3,4,5,6,7},chargeStart_:0,chargeEndOffset_:0]:=
Module[{tuples,dimList,partitionCount,x},
	Print["+++++-----+++++-----+++++-----+++++-----+++++-----+++++"]
	Print["Computing braid = ", braidAndNumStrands, " for u=",uu,"."];
	Print["For conjugacy classes (reps): ", ConjClassReps[[fluxList]]];
	Print["and for projective irreps: ", chargeStart," to ",(NumProjIrreps[ConjClassReps[[#]]]&/@fluxList)-1-chargeEndOffset,"."];
	Print["Using kernels: ", Kernels[]];
	tuples=MakeLoopedTuples[braidAndNumStrands[[1]],braidAndNumStrands[[2]],uu,fluxList,chargeStart,chargeEndOffset];
	dimList=NumProjIrreps[#]&/@ConjClassReps;
	partitionCount=Sum[(dimList[[fluxList[[i]]]]-chargeStart-chargeEndOffset),{i,1,Length[fluxList]}];
	If[tuples[[1]]==2,
		Print["Reshape to ", partitionCount,"x",partitionCount, " matrix"];
		x=(ParallelMap[TraceActSeq,tuples[[2]]]~Partition~partitionCount)//AbsoluteTiming,
		x=ParallelMap[TraceActSeq,tuples[[2]]]//AbsoluteTiming
	];
	Print["Elapsed time: ",x[[1]],"s = ",x[[1]]/60,"min = ",x[[1]]/60/60,"hr."];
	Return[x[[2]]];
];


End[]
EndPackage[]
