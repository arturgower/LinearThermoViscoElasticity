(* ::Package:: *)

(* ::Code::Bold:: *)
(**)


(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
BeginPackage["LoadThermalWavesAlegbra`"];

ClearAll @@ Names["LoadThermalWavesAlegbra`*"];

LoadTWEigensystem::usage = "LoadTWEigensystem returns the formulas for {{\[Lambda][\[Xi]],\[Lambda][\[Zeta]],\[Lambda][\[Psi]]},{V[\[Xi]],V[\[Zeta]]}}  for general thermo-acoustic waves"
CalculateTWEigensystem::usage = "CalculateTWEigensystem calculates the formulas for {{\[Lambda][\[Xi]],\[Lambda][\[Zeta]],\[Lambda][\[Psi]]},{V[\[Xi]],V[\[Zeta]]}}  for general thermo-acoustic waves"

Begin["`Private`"]

LoadTWEigensystem["ECAH"]= Module[{c2,\[Alpha],kp,kt2,B},
  c2 =( 4"\[Mu]"+3 "\[Gamma]" "\[Beta]")/(3"\[Rho]");
  \[Alpha] = \[FormalW]^2/(2 c2^(3/2) "\[Rho]") ("\[Eta]\[Beta]"+4/3 "\[Eta]\[Mu]"+(("\[Gamma]"-1)/("Cv""\[Gamma]")) "\[Kappa]"(1- (4"\[Mu]")/(3 c2 "\[Rho]")));
  kp = Sqrt[\[FormalW]^2/c2] + I \[Alpha];
  kt2 = I "Cp" "\[Rho]" \[FormalW]/"\[Kappa]";
  B = "\[Beta]"-I \[FormalW] "\[Eta]\[Beta]" +4"\[Mu]"/3 - 4 I "\[Eta]\[Mu]" \[FormalW]/3;
 Return@{
    { kp^2, kt2, ("\[Rho]" \[FormalW]^2)/("\[Mu]"-I "\[Eta]\[Mu]" \[FormalW])}
  , {
      {1,("\[Rho]" \[FormalW]^2 - B kp^2)/("T0""\[Alpha]""\[Beta]")},
      {1,("\[Rho]" \[FormalW]^2 - B kt2)/("T0""\[Alpha]""\[Beta]")} 
    }
}
];

LoadTWEigensystem["ArtQ"]= Module[{Q,\[Delta],kp2,kt2,B},
  B = "\[Beta]"-I \[FormalW] "\[Eta]\[Beta]" +4"\[Mu]"/3 - 4 I "\[Eta]\[Mu]" \[FormalW]/3;
  Q =("\[Gamma]"-1)(("\[Beta]")/B-1);
  \[Delta] = - ((I \[FormalW] "\[Kappa]")/("Cv" B));
  kp2 = ("\[Rho]"\[FormalW]^2)/B 1/(Q+"\[Gamma]") (1 -( Q+"\[Gamma]"-1)/(Q+"\[Gamma]")^2 \[Delta]);
  kt2 = (I "\[Rho]"\[FormalW] "Cv")/("\[Kappa]") (Q+"\[Gamma]")(1 +( Q+"\[Gamma]"-1)/(Q+"\[Gamma]")^2 \[Delta]);
 Return@{
    { kp2, kt2, ("\[Rho]" \[FormalW]^2)/("\[Mu]"-I "\[Eta]\[Mu]" \[FormalW])}
  , {
      {1,("\[Rho]" \[FormalW]^2 - B kp2)/("T0""\[Alpha]""\[Beta]")}
    , {1,("\[Rho]" \[FormalW]^2 - B kt2)/("T0""\[Alpha]""\[Beta]")} }
}
];

LoadTWEigensystem["Exact"]= {{
    (\[FormalW] (-3 "T0" ("\[Alpha]")^2 ("\[Beta]")^2-"Cv" (3 "\[Beta]"+4 "\[Mu]") "\[Rho]"+I (3 "Cv" "\[Eta]\[Beta]"+4 "Cv" "\[Eta]\[Mu]"+3 "\[Kappa]") "\[Rho]" \[FormalW]))/(2 "\[Kappa]" (3 I "\[Beta]"+4 I "\[Mu]"+3 "\[Eta]\[Beta]" \[FormalW]+4 "\[Eta]\[Mu]" \[FormalW]))+Sqrt[(3 "Cv" ("\[Rho]")^2 \[FormalW]^3)/(3 I "\[Beta]" "\[Kappa]"+4 I "\[Kappa]" "\[Mu]"+3 "\[Eta]\[Beta]" "\[Kappa]" \[FormalW]+4 "\[Eta]\[Mu]" "\[Kappa]" \[FormalW])+(\[FormalW]^2 (3 I "T0" ("\[Alpha]")^2 ("\[Beta]")^2+"\[Rho]" (3 "\[Kappa]" \[FormalW]+"Cv" (3 I "\[Beta]"+4 I "\[Mu]"+3 "\[Eta]\[Beta]" \[FormalW]+4 "\[Eta]\[Mu]" \[FormalW])))^2)/(4 ("\[Kappa]")^2 (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW])^2)]
  ,
    1/2 ((I "Cv" "\[Rho]" \[FormalW])/("\[Kappa]")+(3 \[FormalW] (I "T0" ("\[Alpha]")^2 ("\[Beta]")^2+"\[Kappa]" "\[Rho]" \[FormalW]))/("\[Kappa]" (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW]))-Sqrt[-((9 ("T0")^2 ("\[Alpha]")^4 ("\[Beta]")^4+6 "T0" ("\[Alpha]")^2 ("\[Beta]")^2 "\[Rho]" ("Cv" (3 "\[Beta]"+4 "\[Mu]")-I (3 "Cv" "\[Eta]\[Beta]"+4 "Cv" "\[Eta]\[Mu]"+3 "\[Kappa]") \[FormalW])+(3 I "\[Kappa]" "\[Rho]" \[FormalW]+"Cv" "\[Rho]" (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW]))^2)/(("\[Kappa]")^2 (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW])^2))] Abs[\[FormalW]])
  ,
    ("\[Rho]" \[FormalW]^2)/("\[Mu]"-I "\[Eta]\[Mu]" \[FormalW])
  }
, {
  {1,-((6 I "\[Alpha]" "\[Beta]" "\[Rho]" \[FormalW]^3)/("\[Kappa]" (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW]) (-((I "Cv" "\[Rho]" \[FormalW])/("\[Kappa]"))+(3 \[FormalW] ("T0" ("\[Alpha]")^2 ("\[Beta]")^2+I "\[Kappa]" "\[Rho]" \[FormalW]))/("\[Kappa]" (3 I "\[Beta]"+4 I "\[Mu]"+3 "\[Eta]\[Beta]" \[FormalW]+4 "\[Eta]\[Mu]" \[FormalW]))+Sqrt[-((9 ("T0")^2 ("\[Alpha]")^4 ("\[Beta]")^4+6 "T0" ("\[Alpha]")^2 ("\[Beta]")^2 "\[Rho]" ("Cv" (3 "\[Beta]"+4 "\[Mu]")-I (3 "Cv" "\[Eta]\[Beta]"+4 "Cv" "\[Eta]\[Mu]"+3 "\[Kappa]") \[FormalW])+(3 I "\[Kappa]" "\[Rho]" \[FormalW]+"Cv" "\[Rho]" (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW]))^2)/(("\[Kappa]")^2 (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW])^2))] Abs[\[FormalW]])))}
, {-(1/(6 "\[Alpha]" "\[Beta]" "\[Rho]" \[FormalW]^3))I "\[Kappa]" (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW]) ((I "Cv" "\[Rho]" \[FormalW])/("\[Kappa]")-(3 \[FormalW] ("T0" ("\[Alpha]")^2 ("\[Beta]")^2+I "\[Kappa]" "\[Rho]" \[FormalW]))/("\[Kappa]" (3 I "\[Beta]"+4 I "\[Mu]"+3 "\[Eta]\[Beta]" \[FormalW]+4 "\[Eta]\[Mu]" \[FormalW]))+Sqrt[-((9 ("T0")^2 ("\[Alpha]")^4 ("\[Beta]")^4+6 "T0" ("\[Alpha]")^2 ("\[Beta]")^2 "\[Rho]" ("Cv" (3 "\[Beta]"+4 "\[Mu]")-I (3 "Cv" "\[Eta]\[Beta]"+4 "Cv" "\[Eta]\[Mu]"+3 "\[Kappa]") \[FormalW])+(3 I "\[Kappa]" "\[Rho]" \[FormalW]+"Cv" "\[Rho]" (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW]))^2)/(("\[Kappa]")^2 (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW])^2))] Abs[\[FormalW]]),1} }
};
LoadTWEigensystem["LeadingOrder"]= {{
   (3 "Cv" "T0" ("\[Rho]")^2 \[FormalW]^2)/(3 ("pT")^2+"Cv" "T0" "\[Rho]" (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW]))
  ,(I \[FormalW] )/("\[Kappa]") ("Cv" "\[Rho]"+(3 ("pT")^2)/("T0" (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW])))
  ,("\[Rho]" \[FormalW]^2)/("\[Mu]"-I "\[Eta]\[Mu]" \[FormalW])
  }
,{
  {1,(3 "pT" "\[Rho]" \[FormalW]^2)/(3 ("pT")^2+"Cv" "T0" "\[Rho]" (3 "\[Beta]"+4 "\[Mu]"-I (3 "\[Eta]\[Beta]"+4 "\[Eta]\[Mu]") \[FormalW]))}
 ,{1,-((\[FormalW] (3 I ("pT")^2+"Cv" "T0" "\[Rho]" (3 I "\[Beta]"+4 I "\[Mu]"+3 "\[Eta]\[Beta]" \[FormalW]+4 "\[Eta]\[Mu]" \[FormalW])))/(3 "pT" "T0" "\[Kappa]"))}
 }
};

LoadTWEigensystem["EC"]= {{
   -(("\[Rho]" \[FormalW]^2 (3 (-1+"\[Gamma]") "\[Kappa]" \[FormalW]+"Cp" (-6 I "\[Beta]" "\[Gamma]"+3 "\[Eta]\[Beta]" \[FormalW]+4 "\[Eta]\[Mu]" \[FormalW]))^2)/(36 ("Cp")^2 ("\[Beta]")^3 ("\[Gamma]")^3))
  , I "Cp" "\[Rho]" \[FormalW]/"\[Kappa]"
  , ("\[Rho]" \[FormalW]^2)/("\[Mu]"-I "\[Eta]\[Mu]" \[FormalW])
  }
,{
  {1,((-1+"\[Gamma]") "\[Rho]" \[FormalW]^2)/("T0" "\[Alpha]" "\[Beta]" "\[Gamma]")}
 ,{1,-((I "Cp" "\[Rho]" \[FormalW])/("T0" "\[Alpha]" "\[Kappa]"))}
 }
};


(* ::Code::Bold:: *)
(*Code need to generate the formulas for the wavenumbers *)


(* ::Code::Initialization::Bold:: *)
CalculateTWEigensystem := Block[{$Assumptions={\[Alpha]\[Element]Reals, \[FormalW] \[Element] Reals,pT \[Element] Reals,L\[Alpha]>0,L\[Beta]>0,Cv>0,k\[Theta]>0,L>0,L\[Eta]\[Mu]>0,\[Beta]>0,T0>0,\[Kappa]>0,\[Lambda]>0,\[Rho]>0,L\[Mu]>=0,\[Eta]\[Mu]>=0,\[Eta]\[Lambda]>=0,\[Eta]\[Beta]>= 0,\[Mu]>=0,L\[Eta]!=0,\[Alpha]!=0,3 I+3 L\[Eta]\[Lambda]+6 L\[Eta]\[Mu]+4 I L\[Mu]^2!= 0}}, 
  (*When we reach the final equations we will turn these parameters into Strings*)
  parametersToString = (#->ToString@SymbolName[#])&/@{\[Mu],\[Kappa],\[Eta]\[Mu],\[Eta]\[Beta],\[Rho],Cv,T0,\[Alpha],\[Beta],L};
  eqs = {
    I \[FormalW] \[Rho] Cv \[Theta][x/L]+\[Kappa] D[\[Theta][x/L],{x,2}] + I \[FormalW] pT/T0 D[ \[Phi][x/L],{x,2}]==0 ,
    \[FormalW]^2 \[Rho]  \[Phi][x/L]+(\[Lambda] - I \[FormalW] \[Eta]\[Lambda] +2 \[Mu]-2I \[FormalW] \[Eta]\[Mu] ) D[\[Phi][x/L],{x,2}]-  pT  \[Theta][x/L]==0,
    \[FormalW]^2 \[Rho] \[Psi][x/L]+(\[Mu] - I \[FormalW] \[Eta]\[Mu]) D[\[Psi][x/L],{x,2}] == 0
  }/.{ \[Alpha] -> pT/(T0 \[Beta]),\[Lambda]-> \[Beta]-2 \[Mu]/3, \[Eta]\[Lambda]-> \[Eta]\[Beta]-2\[Eta]\[Mu]/3,n_[x/L]:> n[x] }/.{\[Phi][x]-> L^2 \[Phi][x],\[Psi][x]->  L^2 \[Psi][x],\[Phi]''[x]->  L^2 \[Phi]''[x],\[Psi]''[x]-> L^2 \[Psi]''[x]};

  (*Introduce non-dimensional parameters*)

  listStateVars ={\[Theta][x],\[Phi][x], \[Psi][x]};
  listDiv={\[Theta]''[x],\[Phi]''[x],\[Psi]''[x]};
  listWaveNumbers= {k\[Theta]^2,k\[Phi]^2,k\[Psi]^2};

  factors=SeriesCoefficient[eqs[[#,1]]-eqs[[#,2]],{listDiv[[#]],0,1}]&/@{1,2,3};
  eqs=(Simplify[eqs[[#,1]]/factors[[#]]]-Simplify[eqs[[#,2]]/factors[[#]]]  ==0 )&/@{1,2,3};

  subWaveNumbers= Thread[listWaveNumbers->
    ( SeriesCoefficient[#[[1,1]],{#[[2]],0,1}]&/@Transpose[{eqs,listStateVars}])
  ];
  subWaveNumbers= (subWaveNumbers ~ Join~Map[ 1/#[[1]] -> 1/#[[2]]&, subWaveNumbers/.Rule->List]);

  listNonDim= {L\[Phi],L\[Theta]};
  listDiv2={\[Phi]''[x],\[Theta][x]};
  subNonDim= Thread[listNonDim->
    (SeriesCoefficient[#[[1,1]],{#[[2]],0,1}]&/@Transpose[{eqs[[1;;2]],listDiv2}])
  ]//Simplify;
  (*subNonDim[[2]]=subNonDim[[2,1]]\[Rule]- subNonDim[[2,2]];*)

  eqsW=( FullSimplify[#[[1,1]]/.#[[3]]->0] +#[[2]] #[[3]]==0 )&/@Transpose[{eqs,listWaveNumbers,listStateVars}];
  eqsW[[1;;2]]=( FullSimplify[#[[1,1]]/.#[[3]]->0] +#[[2]] #[[3]]==0 )&/@Transpose[{eqsW[[1;;2]],listNonDim,listDiv2}];

  If[Length@FullSimplify[And@@(eqsW~Join~eqs/.subNonDim/.subWaveNumbers)]=!=  3,Print["False, Equations are NOT the same for the given assumptions!"]];

  (*Diagonalize system*)
  A = {Coefficient[#[[1]],\[Phi][x]],Coefficient[#[[1]],\[Theta][x]]}&/@eqsW[[1;;2]];
  B = {Coefficient[#[[1]],\[Phi]''[x]],Coefficient[#[[1]],\[Theta]''[x]]}&/@eqsW[[1;;2]];
  If[ True=!= Simplify[A.{\[Phi][x],\[Theta][x]} + B.{\[Phi]''[x],\[Theta]''[x]}== eqsW[[1;;2]]/.Equal[n_,0]->n],Print["Warning: system is not the same!"]];
  eqa = {a==Simplify[Plus@@Eigensystem[ Inverse[B].A ][[1]]]/2};
  suba = eqa/.Equal-> Rule;
  SubToa = First@Solve[eqa,listNonDim[[2]]];
  subTob = {Sqrt[a^2-k\[Theta]^2 k\[Phi]^2]-> b};
  subb = {b-> Sqrt[a^2-k\[Theta]^2 k\[Phi]^2]};
  {{l[\[Zeta]],l[\[Xi]]},{v[\[Zeta]],v[\[Xi]]}}=Simplify[Simplify[Eigensystem[Inverse[B].A] /.SubToa]/.subTob];
  v[\[Xi]] = v[\[Xi]]/(v[\[Xi]][[1]])//Simplify; (* most common convention*)
  v[\[Zeta]] = v[\[Zeta]]/(v[\[Zeta]][[2]])//Simplify; (* NOT common convention*)
 
(*To check it works uncomment below*)
 (*And[Inverse[B].A.v[\[Xi]] \[Equal]  v[\[Xi]]l[\[Xi]] /.subb/.suba //FullSimplify,
  Inverse[B].A.v[\[Zeta]]  \[Equal] v[\[Zeta]]l[\[Zeta]]/.subb /.suba //FullSimplify]*)
  l[\[Psi]] = \[Rho] \[FormalW]^2/(\[Mu] - I \[Eta]\[Mu] \[FormalW]);
  {{l[\[Xi]],l[\[Zeta]],l[\[Psi]]},{v[\[Xi]],v[\[Zeta]]}}/.subb//.suba/.subWaveNumbers/.{pT-> T0 \[Alpha] \[Beta]}~Join~(subNonDim/.pT-> (T0 \[Alpha] \[Beta]))/.parametersToString
];

LoadTWScattering := Module[{TWIn,TW\[Psi],TW\[Xi],TW\[Zeta],R,T},
 ClearAll["height",TWIn,TW\[Psi],TW\[Xi],TW\[Zeta],R,T,subNtoT];
"ClampedBC"; "FreeBC"; "FixedTempBC";"NoHeatFluxBC";
parameters = {"\[Mu]","\[Beta]","\[Eta]\[Mu]","\[Eta]\[Beta]","\[Rho]","L","Cv","T0","\[Alpha]","\[Kappa]"};
subNto[T] = #->(#<>"-")&/@parameters;
subNto[R] = #->(#<>"+")&/@parameters;

(*Incident wave, chosen to be able to represent any plane wave, i.e. thermalwaves or shear wave*)
TWIn = <|"Name"-> "Incident Wave","rng\[Omega]"->{\[Omega]},"WaveNumbers"->{{k[1],k[2],k[3]}},"WaveVectors"-> {{{V[1],V[2]},{\[CapitalPsi][1],\[CapitalPsi][2],\[CapitalPsi][3]}}},"parameters"->{}|>;

(*Reflected and Transmitted Waves*)
  TW\[Xi][R_] := <|"Name"-> "\[Xi] Wave","rng\[Omega]"->{\[Omega]},"WaveNumbers"->{{k[1],k[2],k\[Xi][R]}},"WaveVectors"-> {{\[Xi]0[R]{V[R,\[Xi],\[Phi]],V[R,\[Xi],\[Theta]]},{0,0,0}}},"parameters"->subNto[R]|>;
  TW\[Zeta][R_] := <|"Name"-> "\[Zeta] Wave","rng\[Omega]"->{\[Omega]},"WaveNumbers"->{{k[1],k[2],k\[Zeta][R]}},"WaveVectors"-> {{\[Zeta]0[R]{V[R,\[Zeta],\[Phi]],V[R,\[Zeta],\[Theta]]},{0,0,0}}},"parameters"->subNto[R]|>;
  TW\[Psi][R_] := <|"Name"-> "Shear Wave","rng\[Omega]"->{\[Omega]},"WaveNumbers"->{{k[1],k[2],k\[Psi][R]}},"WaveVectors"-> {{{0,0},{\[CapitalPsi][R,1],\[CapitalPsi][R,2],\[CapitalPsi][R,3]}}},"parameters"->subNto[R]|>;

  unknowns = {\[CapitalPsi][R,1],\[CapitalPsi][R,2],\[CapitalPsi][R,3],\[Xi]0[R],\[Zeta]0[R]};
  unknowns = unknowns~Join~(unknowns/.R-> T);
  listTWRs = {TW\[Xi]@#,TW\[Zeta]@#,TW\[Psi]@#}&@R;
  listTWTs = {TW\[Xi]@#,TW\[Zeta]@#,TW\[Psi]@#}&@T;

(*Continuity boundary conditions*)
  uS =DisplacementTW[ TWIn, {0,0,-"height"}] + Plus@@ Map[DisplacementTW[#,{0,0,0}]&,listTWRs];
  uS -=(1-"ClampedBC") Plus@@ Map[DisplacementTW[#,{0,0,0}]&,listTWTs];
  uS = (1- "FreeBC")uS ;

  \[Sigma] = StressTW[TWIn,{0,0,-"height"}] + Plus@@ Map[StressTW[#,{0,0,0}]&,listTWRs];
  \[Sigma]-= (1- "FreeBC")Plus@@ Map[StressTW[#,{0,0,0}]&,listTWTs];
  \[Sigma]S = -#.{0,0,1}&/@\[Sigma];
  \[Sigma]S= (1-"ClampedBC")\[Sigma]S;

  \[Theta]RS =\[Theta]TW[ TWIn, {0,0,z-"height"}] + Plus@@Map[\[Theta]TW[#,{0,0,z}]&,listTWRs];
  \[Theta]TS = Plus@@ Map[\[Theta]TW[#,{0,0,z}]&,listTWTs];
  D3\[Theta]S =D[\[Theta]RS,z] - (1 - "NoHeatFluxBC")D[\[Theta]TS,z]/.z-> 0;
  D3\[Theta]S= (1 -"FixedTempBC" )D3\[Theta]S;
  \[Theta]S = \[Theta]RS - (1 -  "FixedTempBC")\[Theta]TS/.z->0;
  \[Theta]S = (1 - "NoHeatFluxBC")\[Theta]S;

  (*constrain the shear waves to have zero divergence*)
  \[Psi]Constraint = #["WaveNumbers"][[1]].#["WaveVectors"][[1,2]]==0 & /@{TW\[Psi][R],TW\[Psi][T]};
  BC1 =  Thread[uS[[1]] ==0 ]~Join~Thread[\[Sigma]S[[1]] ==0 ];
  BC2 ={\[Theta]S[[1]] ==0,D3\[Theta]S[[1]] ==0 };
  eqns =Simplify[ BC1~Join~BC2~Join~\[Psi]Constraint , TimeConstraint -> 2]
];

End[];
EndPackage[];



