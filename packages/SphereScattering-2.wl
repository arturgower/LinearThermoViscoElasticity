(* ::Package:: *)

(* ::Code::Bold:: *)
(**)


(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
BeginPackage["SphereScattering`"];

ClearAll @@ Names["SphereScattering`*"];

$Path = $Path~Join~{NotebookDirectory[],NotebookDirectory[]<>"packages/"};
Get["ThermalPotentials.wl"];
Get["SumsAlgebra.wl"];
(*Get["SumsAlgebra.wl", Path\[Rule] {NotebookDirectory[],NotebookDirectory[]<>"packages/"}];
Get["ThermalPotentials.wl", Path\[Rule] {NotebookDirectory[],NotebookDirectory[]<>"packages/"}];*)
Vorticity\[CapitalPhi]::usage = "takes a 2D vector field and returns a function of (r,\[Theta]) that gives the vorticity."

OutgoingDisplacements::usage = "Returns three functions {u\[CurlyPhi],u\[CurlyTheta],u\[Chi]}, each is a function u\[CurlyPhi][r,\[Theta]] = {\!\(\*SubscriptBox[\(u\[CurlyPhi]\), \(r\)]\),\!\(\*SubscriptBox[\(u\[CurlyPhi]\), \(\[Theta]\)]\)}"
InternalDisplacements::usage = "Returns three functions {u\[CurlyPhi],u\[CurlyTheta],u\[Chi]}, each is a function u\[CurlyPhi][r,\[Theta]] = {\!\(\*SubscriptBox[\(u\[CurlyPhi]\), \(r\)]\),\!\(\*SubscriptBox[\(u\[CurlyPhi]\), \(\[Theta]\)]\)}"
IncidentPressureWaveDisplacement::usage = "";

IncidentPressure::usage = "IncidentPressure[subNo_,\[Omega]_, opts:OptionsPattern[]] returns a function of (r,\[Theta])"

InternalPotentials::usage = ""
ScatteredPotentials::usage = ""
OutgoingPressure::usage = "OutgoingPressure[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]] returns a function of (r,\[Theta])"
OutgoingTemperatures::usage = "";
InternalTemperatures::usage = "";
IncidentPressureWaveTemperature::usage = "";
TotalInternalTemperature::usage = "";
TotalDisplacement::usage="";
TotalInternalDisplacement::usage="";

SolveSphericalBC::usage = "Solves the boundary conditions for a set of external waves and internal waves";
emptyTW::usage = "an empty ThermalWave or ShearWave";
internalSeries::usage = "use internalSeries[\[Phi]] gives a general series expansion in terms of BesselJ functions";
externalSeries::usage = "use externalSeries[\[Phi]] gives a general series expansion in terms of Hankel functions";



(* ::Code::Bold:: *)
(**)


(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
Begin["`Private`"]

standardOptions = {
   "ThermalSystem" -> "ECAH", 
   "radius"-> 1/4 10^-6, 
   "LegendreOrder"->3,
   "IncidentLegendreOrder"->16,
   WorkingPrecision->40,
   "IncidentScatteringCoefficientFunction" -> (1 &), 
   "Verbose"-> False
 }; 

Options[TotalInternalDisplacement]= standardOptions;
Options[TotalDisplacement]= standardOptions;

Options[OutgoingTemperatures]= standardOptions;

Options[OutgoingDisplacements]= standardOptions;
Options[OutgoingPressure]= standardOptions;
Options[ScatteredPotentials]= standardOptions;
Options[InternalPotentials]= standardOptions;

Options[InternalTemperatures]= standardOptions;
Options[TotalInternalTemperature]= standardOptions;
Options[InternalDisplacements]= standardOptions;

Options[IncidentPressure]= standardOptions;
Options[IncidentPressureWaveDisplacement]= standardOptions;
Options[IncidentPressureWaveTemperature]= standardOptions;
Options[SolveSphericalBC]= standardOptions;

   
$Assumptions = {\[FormalR]>= 0, \[FormalN] \[Element] Integers,0<= \[FormalT]<= \[Pi]};
subNto[o_] := #->(#<>o)&/@AllParameters;
subNoPressure = {\[FormalCapitalA]["\[CurlyPhi]",_,_] -> 0};
subNoThermal = {\[FormalCapitalA]["\[CurlyTheta]",_,_] -> 0};
subNoShear = {\[FormalCapitalA]["\[Chi]",_,_] -> 0};

emptyTW = <|"Name" -> "the emptyness", "rng\[Omega]"-> {0.}, 
             "WaveNumbers"->{0.}, "WaveModes"->{0.,0.}, "parameters"->{}|>;


(* ::Code::Bold:: *)
(*NOTE: looks like Valerie's vorticity is negative the below*)


(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
Vorticity\[CapitalPhi][us_]:= Module[{f}, f[\[FormalR]_,\[FormalT]_]= -1/\[FormalR] (us[\[FormalR],\[FormalT]][[2]] + \[FormalR] D[us[\[FormalR],\[FormalT]][[2]],\[FormalR]]- D[us[\[FormalR],\[FormalT]][[1]],\[FormalT]]); Return[f]]; 

IncidentPressureWaveTemperature[subN_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subn, precision = Max[OptionValue@WorkingPrecision, Precision@subN],
NN = OptionValue@"IncidentLegendreOrder", subWaveCoefficients,TW\[CurlyPhi],TW\[CurlyTheta],N\[Tau]0,r,\[Theta]},
  
  subn = SetPrecision[AddToNsub[subN,"ThermalSystem" -> OptionValue["ThermalSystem"]],precision];
  {TW\[CurlyPhi],TW\[CurlyTheta]} = ThermalPotentials[subn,SetPrecision[{\[Omega]},precision], "ThermalSystem" -> OptionValue["ThermalSystem"]];
  subWaveCoefficients = {
    "T\[CurlyPhi]o" -> TW\[CurlyPhi]["WaveModes"][[1,2]],"U\[CurlyPhi]o" -> TW\[CurlyPhi]["WaveModes"][[1,1]], 
    \[FormalK]["\[CurlyPhi]","o"]->TW\[CurlyPhi]["WaveNumbers"][[1]], \[FormalCapitalA]["\[CurlyPhi]0",\[FormalN]_,"o"] -> OptionValue["IncidentScatteringCoefficientFunction"][\[FormalN]]
  };
  N\[Tau]0[r_,\[Theta]_] = Activate[\[Tau]0/.{\[Infinity] ->NN}/.subWaveCoefficients/.(subn/.subNto["o"])/.{\[FormalR]->r,\[FormalT]->\[Theta]}];
 Return[N\[Tau]0]
];

IncidentPressure[subN_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subn, precision = Max[OptionValue@WorkingPrecision, Precision@subN],
NN = OptionValue@"IncidentLegendreOrder", subWaveCoefficients,TW\[CurlyPhi],TW\[CurlyTheta],r,\[Theta]},
  subn = SetPrecision[AddToNsub[subN,"ThermalSystem" -> OptionValue["ThermalSystem"]],precision];
  {TW\[CurlyPhi],TW\[CurlyTheta]} = ThermalPotentials[subn,SetPrecision[{\[Omega]},precision],"ThermalSystem" -> OptionValue["ThermalSystem"]];
  subWaveCoefficients= {
    "T\[CurlyPhi]o" -> TW\[CurlyPhi]["WaveModes"][[1,2]],"U\[CurlyPhi]o" -> TW\[CurlyPhi]["WaveModes"][[1,1]], 
    \[FormalK]["\[CurlyPhi]","o"]->TW\[CurlyPhi]["WaveNumbers"][[1]], 
    \[FormalCapitalA]["\[CurlyPhi]0",\[FormalN]_,"o"] -> OptionValue["IncidentScatteringCoefficientFunction"][\[FormalN]]
  };
  Np[r_,\[Theta]_] = Activate[ -("\[Beta]o"-I "\[Eta]\[Beta]o" \[Omega])(D[\[FormalR]^2 ur0,\[FormalR]]/(\[FormalR]^2) 
        +(u\[Theta]0/.subLegendreDerivativeIdentity)/\[FormalR])+"\[Alpha]o""\[Beta]o"\[Tau]0/.{\[Infinity] -> NN}/.subWaveCoefficients/.(subn/.subNto["o"])/.{\[FormalT]->\[Theta],\[FormalR]->r}];
 Return[Np]
];

IncidentPressureWaveDisplacement[subN_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subn, precision = Max[OptionValue@WorkingPrecision, Precision@subN],
NN = OptionValue@"IncidentLegendreOrder", subWaveCoefficients,TW\[CurlyPhi],TW\[CurlyTheta],Nu\[CurlyPhi], r, \[Theta]},
  
  subn = SetPrecision[AddToNsub[subN,"ThermalSystem" -> OptionValue["ThermalSystem"]], precision];
  {TW\[CurlyPhi],TW\[CurlyTheta]} = ThermalPotentials[subn, SetPrecision[{\[Omega]},precision],"ThermalSystem" -> OptionValue["ThermalSystem"]];
  subWaveCoefficients= {\[FormalK]["\[CurlyPhi]","o"]->TW\[CurlyPhi]["WaveNumbers"][[1]], "U\[CurlyPhi]o" -> TW\[CurlyPhi]["WaveModes"][[1,1]], 
  \[FormalCapitalA]["\[CurlyPhi]0",\[FormalN]_,"o"] -> OptionValue["IncidentScatteringCoefficientFunction"][\[FormalN]]};
  Nu\[CurlyPhi][r_,\[Theta]_] = Activate[{ur0,u\[Theta]0}/.{\[Infinity] -> NN}/.subWaveCoefficients/.(subn/.subNto["o"])/.{\[FormalR]->r,\[FormalT]->\[Theta]}];
  If[OptionValue@"IncidentLegendreOrder", Print["LegendreOrder =",NN]];
 Return[Nu\[CurlyPhi]]
];

TotalInternalTemperature[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subN,subsols, N = OptionValue@"LegendreOrder",NTi, r, \[Theta]},
  {subsols,subN} = SolveSphericalBC[subNo,subNi,\[Omega], FilterRules[{opts}, Options[SolveSphericalBC]]];
  NTi[r_,\[Theta]_] = Activate[Ti/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.Flatten[subsols];
 Return[NTi]
];

InternalTemperatures[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subN,subsols, N = OptionValue@"LegendreOrder",NT\[CurlyPhi],NT\[CurlyTheta],NT\[Chi], r, \[Theta]},
  {subsols,subN} = SolveSphericalBC[subNo,subNi,\[Omega], FilterRules[{opts}, Options[SolveSphericalBC]]];

  NT\[CurlyPhi][r_,\[Theta]_] = Activate[Ti/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoThermal/.subNoShear/.Flatten[subsols];
  NT\[CurlyTheta][r_,\[Theta]_] = Activate[Ti/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoPressure/.subNoShear/.Flatten[subsols];
  NT\[Chi][r_,\[Theta]_] = Activate[Ti/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoPressure/.subNoThermal/.Flatten[subsols];
 Return[{NT\[CurlyPhi],NT\[CurlyTheta],NT\[Chi]}]
];

TotalInternalDisplacement[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subN,subsols, N = OptionValue@"LegendreOrder",Nu, r, \[Theta]},
  {subsols,subN} = SolveSphericalBC[subNo,subNi,\[Omega], FilterRules[{opts}, Options[SolveSphericalBC]]];
  Nu[r_,\[Theta]_] = Activate[{uri,u\[Theta]i}/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.Flatten[subsols];
 Return[Nu]
];

InternalDisplacements[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subN,subsols, N = OptionValue@"LegendreOrder",Nu\[CurlyPhi],Nu\[CurlyTheta],Nu\[Chi],r,\[Theta]},
  {subsols,subN} = SolveSphericalBC[subNo,subNi,\[Omega], FilterRules[{opts}, Options[SolveSphericalBC]]];

  Nu\[CurlyPhi][r_,\[Theta]_] = Activate[{uri,u\[Theta]i}/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoThermal/.subNoShear/.Flatten[subsols];
  Nu\[CurlyTheta][r_,\[Theta]_] = Activate[{uri,u\[Theta]i}/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoPressure/.subNoShear/.Flatten[subsols];
  Nu\[Chi][r_,\[Theta]_] = Activate[{uri,u\[Theta]i}/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoPressure/.subNoThermal/.Flatten[subsols];
 Return[{Nu\[CurlyPhi],Nu\[CurlyTheta],Nu\[Chi]}]
];

OutgoingTemperatures[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subN,subsols, N = OptionValue@"LegendreOrder",NT\[CurlyPhi],NT\[CurlyTheta],NT\[Chi], r, \[Theta]},
  {subsols,subN} = SolveSphericalBC[subNo,subNi,\[Omega], FilterRules[{opts}, Options[SolveSphericalBC]]];

  NT\[CurlyPhi][r_,\[Theta]_] = Activate[To/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoThermal/.subNoShear/.Flatten[subsols];
  NT\[CurlyTheta][r_,\[Theta]_] = Activate[To/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoPressure/.subNoShear/.Flatten[subsols];
  NT\[Chi][r_,\[Theta]_] = Activate[To/.{\[Infinity] -> N}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoPressure/.subNoThermal/.Flatten[subsols];
 Return[{NT\[CurlyPhi],NT\[CurlyTheta],NT\[Chi]}]
];

OutgoingDisplacements[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subN,subsols, NN = OptionValue@"LegendreOrder",Nu\[CurlyPhi],Nu\[CurlyTheta],Nu\[Chi], r, \[Theta]},
  {subsols,subN} = SolveSphericalBC[subNo,subNi,\[Omega], FilterRules[{opts}, Options[SolveSphericalBC]]];

  Nu\[CurlyPhi][r_,\[Theta]_] = Activate[{uro,u\[Theta]o}/.{\[Infinity] -> NN}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoThermal/.subNoShear/.Flatten[subsols];
  Nu\[CurlyTheta][r_,\[Theta]_] = Activate[{uro,u\[Theta]o}/.{\[Infinity] -> NN}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoPressure/.subNoShear/.Flatten[subsols];
  Nu\[Chi][r_,\[Theta]_] = Activate[{uro,u\[Theta]o}/.{\[Infinity] -> NN}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.subNoPressure/.subNoThermal/.Flatten[subsols];
  If[OptionValue@"LegendreOrder", Print["LegendreOrder =",NN]];
 Return[{Nu\[CurlyPhi],Nu\[CurlyTheta],Nu\[Chi]}]
];

TotalDisplacement[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subN,subsols, NN = OptionValue@"LegendreOrder", Nu, Nuo, Nu0,r, \[Theta]},
  {subsols,subN} = SolveSphericalBC[subNo,subNi,\[Omega], FilterRules[{opts}, Options[SolveSphericalBC]]];
  Nu0 = IncidentPressureWaveDisplacement[subNo,\[Omega], opts];
  Nuo[r_,\[Theta]_] = Activate[{uro,u\[Theta]o}/.{\[Infinity] ->NN}/.subN/.{\[FormalR]->r,\[FormalT]->\[Theta]}]/.Flatten[subsols];
  Nu[r_,\[Theta]_] = Nu0[r,\[Theta]] + Nuo[r,\[Theta]];
 Return[Nu]
];

OutgoingPressure[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subN,subsols, N = OptionValue@"LegendreOrder",Np,Nu\[CurlyTheta],Nu\[Chi], \[Theta], r},
  {subsols,subN} = SolveSphericalBC[subNo,subNi,\[Omega], FilterRules[{opts}, Options[SolveSphericalBC]]];
  Np[r_,\[Theta]_] = Activate[ -("\[Beta]o"-I "\[Eta]\[Beta]o" \[Omega])(D[\[FormalR]^2 uro,\[FormalR]]/\[FormalR]^2 
      +(u\[Theta]o/.subLegendreDerivativeIdentity)/\[FormalR])+"\[Alpha]o""\[Beta]o"To/.{\[Infinity] -> N}/.subN/.{\[FormalT]->\[Theta],\[FormalR]->r}]/.subNoShear/.Flatten[subsols];
 Return[Np]
];

ScatteredPotentials[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subN,subsols, N = OptionValue@"LegendreOrder",waves, \[Theta], r},
  {subsols,subN} = SolveSphericalBC[subNo,subNi,\[Omega], FilterRules[{opts}, Options[SolveSphericalBC]]];
  waves[r_,\[Theta]_] = Activate[{outgoingSeries["\[CurlyPhi]"],outgoingSeries["\[CurlyTheta]"],outgoingSeries["\[Chi]"]}/.{\[Infinity] -> N}/.subN/.{\[FormalT]->\[Theta],\[FormalR]->r}]/.Flatten[subsols];
 Return[waves]
];
InternalPotentials[subNo_,subNi_,\[Omega]_, opts:OptionsPattern[]]:= 
Module[{subN,subsols, N = OptionValue@"LegendreOrder",waves, \[Theta], r},
  {subsols,subN} = SolveSphericalBC[subNo,subNi,\[Omega], FilterRules[{opts}, Options[SolveSphericalBC]]];
  waves[r_,\[Theta]_] = Activate[{internalSeries["\[CurlyPhi]"],internalSeries["\[CurlyTheta]"],internalSeries["\[Chi]"]}/.{\[Infinity] -> N}/.subN/.{\[FormalT]->\[Theta],\[FormalR]->r}]/.Flatten[subsols];
 Return[waves]
];


(* ::Code::Bold:: *)
(**)


(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
subNQ = (Head@# == List && Head[#[[1]]] == Rule)&;
SolveSphericalBC[subNo_?subNQ, subNi_?subNQ, \[Omega]_, opts:OptionsPattern[]]:= 
Module[{precision = Max[OptionValue@WorkingPrecision, Precision@subNo], 
subno, subni, rng\[Omega], TW\[CurlyPhi],TW\[CurlyTheta],TW\[Psi],TWi\[CurlyPhi],TWi\[CurlyTheta],TWi\[Psi]},

  subno = SetPrecision[AddToNsub[subNo,"ThermalSystem" -> OptionValue["ThermalSystem"]],precision];
  subni = SetPrecision[AddToNsub[subNi,"ThermalSystem" -> OptionValue["ThermalSystem"]],precision];

  {TW\[CurlyPhi],TW\[CurlyTheta]} = ThermalPotentials[subno,{\[Omega]},"ThermalSystem" -> OptionValue["ThermalSystem"]];
  TW\[Psi] = ShearPotential[subno,{\[Omega]},"ThermalSystem" -> OptionValue["ThermalSystem"]];
 
  {TWi\[CurlyPhi],TWi\[CurlyTheta]} = ThermalPotentials[subni,{\[Omega]},"ThermalSystem" -> OptionValue["ThermalSystem"]];
  TWi\[Psi] = ShearPotential[subni,{\[Omega]},"ThermalSystem" -> OptionValue["ThermalSystem"]];
  
  Return[SolveSphericalBC[{TW\[CurlyPhi],TW\[CurlyTheta],TW\[Psi]},{TWi\[CurlyPhi],TWi\[CurlyTheta],TWi\[Psi]},opts]]
]; 

SolveSphericalBC[{TW\[CurlyPhi]:(_?AssociationQ):emptyTW,TW\[CurlyTheta]_:emptyTW,TW\[Psi]_:emptyTW},
                 {TWi\[CurlyPhi]:(_?AssociationQ):emptyTW,TWi\[CurlyTheta]_:emptyTW,TWi\[Psi]_:emptyTW}
                 ,opts:OptionsPattern[]]:= 
Module[{subN,\[Omega],subNo,subNi,subsols, subincidentPlanewave, n, residuals, Nb, NM,
        NN = OptionValue@"LegendreOrder" , a0, Ncoefs, undeterminedAs, subAs,
        precision = OptionValue@WorkingPrecision},
 
  If[Precision["\[Mu]"/.TW\[CurlyPhi]["parameters"]] > precision, precision = Precision["\[Mu]"/.TW\[CurlyPhi]["parameters"]]];
  a0 = SetPrecision[OptionValue@"radius", precision];
  subincidentPlanewave= {\[FormalCapitalA]["\[CurlyPhi]0",n_,"o"] :> OptionValue["IncidentScatteringCoefficientFunction"][n]};

   \[Omega] = TW\[CurlyPhi]["rng\[Omega]"][[1]];
   subN = {"radius"->a0}~Join~TWsTosubParameters[{TW\[CurlyPhi],TW\[CurlyTheta],TW\[Psi]},{TWi\[CurlyPhi],TWi\[CurlyTheta],TWi\[Psi]}];

   NM = M/.subN;
   Nb = b/.subN/.subincidentPlanewave;

   subsols = With[{cs = coefs/.\[FormalM]->#, eqs =Thread[NM.coefs + Nb == 0]/.\[FormalM]->#}, 
      NSolve[eqs, cs, WorkingPrecision -> precision]
   ]&/@Range[0,NN];

   undeterminedAs = Union@Cases[Transpose[Flatten[subsols]/.Rule-> List][[2]],\[FormalCapitalA][n__]:>\[FormalCapitalA][n],\[Infinity]];
   If[Length@undeterminedAs > 0, 
     Ncoefs = coefs/.\[FormalM]->#&/@Range[0,NN];
     undeterminedAs = Union[undeterminedAs,Flatten@Cases[Ncoefs/.Flatten@subsols,\[FormalCapitalA][n__]:>\[FormalCapitalA][n],\[Infinity]]];   
     If[ OptionValue["Verbose"],
         Print["Will take the coefficients below to be zero, as they could not determined from the boundary conditions: \n", undeterminedAs]
     ];    
     subAs = Thread[undeterminedAs -> 0];
     subsols = (subsols/.subAs)~Join~subAs;
   ];

   (*residuals = (Norm[NM.coefs + Nb/.m->#/.Flatten[subsols[[#+1]]]/.subAs]/Norm[Nb/.m->#]) &/@Range[0,NN];*)
   residuals = (Norm[NM.coefs + Nb/.\[FormalM]->#/.Flatten[subsols]]/Norm[Nb/.\[FormalM]->#]) &/@Range[0,NN];
   If[ OptionValue["Verbose"] || Norm[residuals]>10^-6,
       Print["Relative residuals of the boundary condition system ( each residual corresponds to one of the Legendre Polynomial)  = \n",residuals]
       (*Print["The condition numbers of the system cond(M)= \n",LinearAlgebra`MatrixConditionNumber[NM/.m->#]&/@Range[0,NN]];*)
   ];
   
 Return[{subsols,subN}]
];
TWsTosubParameters[TWos_,TWis_] := (
    {
    \[FormalK]["\[CurlyPhi]","o"]-> TWos[[1]]["WaveNumbers"][[1]], \[FormalK]["\[CurlyTheta]","o"]-> TWos[[2]]["WaveNumbers"][[1]], \[FormalK]["\[Chi]","o"]-> TWos[[3]]["WaveNumbers"][[1]],
    "T\[CurlyPhi]o" -> TWos[[1]]["WaveModes"][[1,2]], "T\[CurlyTheta]o" -> TWos[[2]]["WaveModes"][[1,2]], 
    "U\[CurlyPhi]o" -> TWos[[1]]["WaveModes"][[1,1]], "U\[CurlyTheta]o" -> TWos[[2]]["WaveModes"][[1,1]],
    \[FormalK]["\[CurlyPhi]","i"]-> TWis[[1]]["WaveNumbers"][[1]], \[FormalK]["\[CurlyTheta]","i"]-> TWis[[2]]["WaveNumbers"][[1]], \[FormalK]["\[Chi]","i"]-> TWis[[3]]["WaveNumbers"][[1]],
    "T\[CurlyPhi]i" -> TWis[[1]]["WaveModes"][[1,2]], "T\[CurlyTheta]i" -> TWis[[2]]["WaveModes"][[1,2]],
    "U\[CurlyPhi]i" -> TWis[[1]]["WaveModes"][[1,1]], "U\[CurlyTheta]i" -> TWis[[2]]["WaveModes"][[1,1]]
  }~Join~(TWos[[1]]["parameters"]/.subNto["o"])~Join~(TWis[[1]]["parameters"]/.subNto["i"])
);


(* ::Code::Bold:: *)
(* To check identities, see section Test identities at the end of this document*)


(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
subLegendreDerivativeIdentity = {LegendreP[\[FormalM]_,1,Cos[\[FormalT]]]-> -\[FormalM](\[FormalM]+1)LegendreP[\[FormalM],Cos[\[FormalT]]]} ;
(*Applies the derivative and then subsitutes the identity Eq (1.28) in SphereScattering.pdf*) 

DifferentiateTheta::\[Theta]Dependance="The argument `1` should only depend on \[Theta] through terms of the form LegendreP[m,Cos[\[Theta]]].";
DifferentiateTheta[s_] := (
    If[ (s/.LegendreP[\[FormalM]_,Cos[\[FormalT]]]->\[FormalL][\[FormalM]]/.\[FormalT]-> RandomReal[\[Pi]]/.\[FormalL][\[FormalM]_]->LegendreP[\[FormalM],Cos[\[FormalT]]] )=!= s, Message[DifferentiateTheta::\[Theta]Dependance,s]];
    s/.{LegendreP[\[FormalM]_,Cos[\[FormalT]]]-> LegendreP[\[FormalM],1,Cos[\[FormalT]]]}
);
(* Assumes the only \[Theta] dependance is through LegendreP[m,Cos[\[Theta]]], and substitutes LegendreP[m,Cos[\[Theta]]] \[Rule] D[LegendreP[m,Cos[\[Theta]]],\[Theta]] = LegendreP[m,1,Cos[\[Theta]]]*) 

DivideByLegendreP::PDependance="The argument `1` be LegendreP[n_,Cos[\[Theta]]] times something.";
DivideByLegendreP[s_] := Module[{m,i,n,j,legendre,list},
  list = Cases[s, (LegendreP[n_,Cos[\[FormalT]]])->{n,0},\[Infinity]]\[Union]Cases[s, LegendreP[n_,j_,Cos[\[FormalT]]]-> {n,j},\[Infinity]] ;
  If[Length[list]>1,  Message[DivideByLegendreP::PDependance,s]];
  legendre = LegendreP[ list[[1,1]],list[[1,2]],Cos[\[FormalT]]];

  If[ Simplify[Coefficient[s,legendre] legendre == s] =!= True, Message[DivideByLegendreP::PDependance,s]];
  Coefficient[s,legendre]
];

displacement[waveSeries_]:=Module[{ur,u\[Theta]},
  ur = - "U\[CurlyPhi]"D[waveSeries["\[CurlyPhi]"],\[FormalR]]- "U\[CurlyTheta]"D[waveSeries["\[CurlyTheta]"],\[FormalR]] - 1/\[FormalR](DifferentiateTheta[waveSeries["\[Chi]"]] /.subLegendreDerivativeIdentity);
  u\[Theta] = - (1/\[FormalR])"U\[CurlyPhi]" DifferentiateTheta[waveSeries["\[CurlyPhi]"]] - (1/\[FormalR])"U\[CurlyTheta]" DifferentiateTheta[waveSeries["\[CurlyTheta]"]] + (1/\[FormalR]) D[\[FormalR] DifferentiateTheta[waveSeries["\[Chi]"]], \[FormalR]] ;
  Return[SimplifySum/@{ur,u\[Theta]}]
];
traction[u_,T_]:=Module[{\[Sigma]rr,\[Sigma]r\[Theta],ur,u\[Theta]},
  ur=u[[1]]; u\[Theta]=u[[2]];
  \[Sigma]rr = 2 ("\[Lambda]")ur/\[FormalR] + ("\[Lambda]"+2"\[Mu]")\[FormalR] D[ur,\[FormalR]] + ("\[Lambda]")(1/\[FormalR])(u\[Theta]/.subLegendreDerivativeIdentity) - "pT" T;
  \[Sigma]r\[Theta] = "\[Mu]" D[u\[Theta],\[FormalR]] + ("\[Mu]")/\[FormalR] (DifferentiateTheta[ur] - u\[Theta]) ;
  Return[SimplifySum/@{\[Sigma]rr,\[Sigma]r\[Theta]}]
];


(* ::Code::Bold:: *)
(*Incident wave*)


(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
 Unprotect[b,M,coef0s,coefs]//Quiet;
 Unprotect[\[Sigma]rri,\[Sigma]r\[Theta]i,\[Sigma]rr,\[Sigma]r\[Theta],Ti,To,Temp,\[Tau]0,uri,u\[Theta]i,ur,u\[Theta],uro,u\[Theta]o,ur0,u\[Theta]0]//Quiet;
 Unprotect[incidentSeries,outgoingSeries,internalSeries,subincident]//Quiet;
  
   
  incidentSeries[\[Phi]_] := Inactive[Sum][ I^\[FormalN] (2\[FormalN]+1) \[FormalCapitalA][\[Phi],\[FormalN],"o"]Inactive[SphericalBesselJ][\[FormalN],\[FormalK][\[Phi],"o"]\[FormalR]] LegendreP[\[FormalN],Cos[\[FormalT]]],{\[FormalN],0,\[Infinity]}];
(*Outside waves*)
  outgoingSeries[\[Phi]_] := Inactive[Sum][ I^\[FormalN] (2\[FormalN]+1) \[FormalCapitalA][\[Phi],\[FormalN],"o"]Inactive[SphericalHankelH1][\[FormalN],\[FormalK][\[Phi],"o"]\[FormalR]]/SphericalHankelH1[\[FormalN],\[FormalK][\[Phi],"o"]"radius"] LegendreP[\[FormalN],Cos[\[FormalT]]],{\[FormalN],0,\[Infinity]}];
(*Inside waves*)
  internalSeries[\[Phi]_] := Inactive[Sum][ I^\[FormalN] (2\[FormalN]+1) \[FormalCapitalA][\[Phi],\[FormalN],"i"]Inactive[SphericalBesselJ][\[FormalN],\[FormalK][\[Phi],"i"]\[FormalR]]/SphericalBesselJ[\[FormalN],\[FormalK][\[Phi],"i"]"radius"] LegendreP[\[FormalN],Cos[\[FormalT]]],{\[FormalN],0,\[Infinity]}];

(*Incident pressure dominated waves*)
  subincident= {\[FormalCapitalA]["\[Chi]",_,"o"]->0, \[FormalCapitalA]["\[CurlyTheta]",_,"o"]-> 0, \[FormalCapitalA]["\[CurlyPhi]",\[FormalN]_,"o"]:> \[FormalCapitalA]["\[CurlyPhi]0",\[FormalN],"o"]};
  {ur0,u\[Theta]0} = displacement[incidentSeries]/.subNto["o"]/.subincident;
  (*Incident temperature field*) 
  \[Tau]0 = - SimplifySum["T0""T\[CurlyPhi]" incidentSeries["\[CurlyPhi]"]/.subNto["o"]]/.subincident;

(*Outgoing waves*)
  {uro,u\[Theta]o} = displacement[outgoingSeries]/.subNto["o"];
  ur = SimplifySum[uro + ur0];
  u\[Theta] = SimplifySum[u\[Theta]o + u\[Theta]0];

(*Outgoing temperature field*)
  To = -SimplifySum["T0""T\[CurlyPhi]" outgoingSeries["\[CurlyPhi]"] + "T0""T\[CurlyTheta]" outgoingSeries["\[CurlyTheta]"]/.subNto["o"]];
  Temp = SimplifySum[\[Tau]0 + To];

  {\[Sigma]rr,\[Sigma]r\[Theta]} = traction[{ur,u\[Theta]},Temp]/.subNto["o"];

(*Internal waves*)
  Ti = -SimplifySum["T0""T\[CurlyPhi]" internalSeries["\[CurlyPhi]"] + "T0""T\[CurlyTheta]" internalSeries["\[CurlyTheta]"]/.subNto["i"]];

  {uri,u\[Theta]i} = displacement[internalSeries]/.subNto["i"];
  {\[Sigma]rri,\[Sigma]r\[Theta]i} = traction[{uri,u\[Theta]i},Ti]/.subNto["i"];

(*Apply the boundary conditions*)
roots = {};
AppendTo[roots,DivideByLegendreP[ur-uri/.subpickOut[\[FormalN],\[FormalM]]] ];
AppendTo[roots,DivideByLegendreP[u\[Theta]-u\[Theta]i/.subpickOut[\[FormalN],\[FormalM]]] ];
AppendTo[roots,DivideByLegendreP[\[Sigma]rr-\[Sigma]rri/.subpickOut[\[FormalN],\[FormalM]]] ];
AppendTo[roots,DivideByLegendreP[\[Sigma]r\[Theta]-\[Sigma]r\[Theta]i/.subpickOut[\[FormalN],\[FormalM]]] ];
AppendTo[roots,DivideByLegendreP[Temp-Ti/.subpickOut[\[FormalN],\[FormalM]]] ];
AppendTo[roots,DivideByLegendreP["\[Kappa]o"D[Temp,\[FormalR]]-"\[Kappa]i"D[Ti,\[FormalR]]/.subpickOut[\[FormalN],\[FormalM]]] ];

(*just to simplify*)
terms = {1/\[FormalR] I^\[FormalM] (1+2\[FormalM]), 1/\[FormalR] I^\[FormalM] (1+2\[FormalM]),I^\[FormalM](1+2\[FormalM]), 1/\[FormalR]^2 I^\[FormalM] (1+2\[FormalM]), I^\[FormalM] (1+2\[FormalM]), I^\[FormalM] (1+2\[FormalM])};
roots =(roots/terms /. \[FormalR]->"radius") //Simplify;

(*unknowns to solve for*)
coefs = {\[FormalCapitalA]["\[CurlyPhi]",\[FormalM],"o"],\[FormalCapitalA]["\[CurlyPhi]",\[FormalM],"i"],\[FormalCapitalA]["\[CurlyTheta]",\[FormalM],"o"],\[FormalCapitalA]["\[CurlyTheta]",\[FormalM],"i"],\[FormalCapitalA]["\[Chi]",\[FormalM],"o"],\[FormalCapitalA]["\[Chi]",\[FormalM],"i"]};
(*coefficients of the incident wave*)
coef0s = {\[FormalCapitalA]["\[CurlyPhi]0",\[FormalM],"o"]};
M = Outer[Coefficient[#1,#2]&,roots,coefs]//Simplify;
b0 = Outer[Coefficient[#1,#2]&,roots,coef0s]//Simplify;
Print["Extracted coefficients correctly: ", M.coefs + b0.coef0s == roots //Simplify];

b = Simplify[Activate[b0.coef0s]];
M = Simplify[Activate[M]];

Protect[b,M,coef0s,coefs];
Protect[\[Sigma]rri,\[Sigma]r\[Theta]i,\[Sigma]rr,\[Sigma]r\[Theta],Ti,To,Temp,\[Tau]0,uri,u\[Theta]i,ur,u\[Theta],uro,u\[Theta]o,ur0,u\[Theta]0];
Protect[incidentSeries,outgoingSeries,internalSeries,subincident];

End[];
EndPackage[];


