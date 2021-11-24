(* ::Package:: *)

(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
BeginPackage["ThermalPotentials`"];

Get["LoadThermalWavesAlgebra.wl"]

(*ThermalWavesAlegbra*)
ThermalPotentials::usage = "returns to association, first for the pressure dominated and second for the thermal dominated."
ShearPotential::usage = "return an association for the shear wave, which includes the wavenumbers, and parameters used."
AddToNsub::usage = "adds missing parameters which are calculated from other parameters"
AllParameters::usage = "variable with all possible parameters"
parameterEquations::usage = "All "

(*ParaToExp::usage = "sub all parameters from a string to symbol"
ExpToPara::usage = "sub all parameter symbols to strings"*)
subParaToExp::usage = "sub all parameters from a string to symbol"
subExpToPara::usage = "sub all parameter symbols to strings"
TWEigensystem::usage = "export LoadTWEigensystem from LoadThermalWavesAlgebra.wl"


(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
Begin["`Private`"]

TWEigensystem = LoadTWEigensystem["ECAH"];

standardOptions = {"ThermalSystem"-> "ArtQ", "Verbose"-> False}; (*"Exact", "LeadingOrder" or "EC" *)
Options[ThermalPotentials] = standardOptions;
Options[ShearPotential] = standardOptions;
Options[AddToNsub] = standardOptions;


(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
parameterEquations = {
    0== - "\[Lambda]" + "\[Beta]"-2/3 "\[Mu]",
    0== "\[Gamma]" "\[Beta]" +4 "\[Mu]"/3 -  "\[Rho]""cp2" ,(*approximate compressional low frequency wave speed squared*)
    0== -"\[Eta]\[Lambda]" + "\[Eta]\[Beta]"-2/3 "\[Eta]\[Mu]", 0== -"pT" + "\[Alpha]" "\[Beta]" "T0",
    0==-"\[Gamma]" + "Cp"/"Cv", 0 == -"\[Rho]""Cp" + "\[Rho]""Cv" +"\[Alpha]""pT"
};

AllParameters= {"\[Mu]", "\[Lambda]","\[Beta]","cp2","\[Kappa]","\[Eta]\[Mu]","\[Eta]\[Beta]","\[Rho]","Cv","Cp","U\[CurlyPhi]","U\[CurlyTheta]","T\[CurlyPhi]","T\[CurlyTheta]","T0","\[Alpha]","L","\[Gamma]","pT","\[Eta]\[Lambda]"} ;
subParaToExp := Thread[AllParameters -> ToExpression/@AllParameters] ;
subExpToPara := Thread[ ToExpression/@AllParameters -> AllParameters];

(*ParaToExp[subn_] := Thread[Transpose[subn/.Rule-> List][[1]]->ToExpression/@Transpose[subn/.Rule-> List][[1]]];
ExpToPara[subn_] := Thread[ToExpression/@Transpose[subn/.Rule-> List][[1]]\[Rule] Transpose[subn/.Rule-> List][[1]]];
*)


(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
(* where \[Xi] \[Equal] \[ExponentialE]^(\[PlusMinus] \[ImaginaryI] Sqrt[\[Lambda][\[Xi]]]x.n), \[Zeta] \[Equal] \[ExponentialE]^(\[PlusMinus] \[ImaginaryI] Sqrt[\[Lambda][\[Zeta]]]x.n), \[VerticalSeparator]n\[VerticalSeparator] \[Equal] 1, 
and 
{\[Phi], \[Theta]} \[Equal]  \[Xi] V[\[Xi]] + \[Zeta] V[\[Zeta]] *)

OutgoingWaveNumber = 
  If[ Im@#<0, 
       Print["Wave amplitude increases in time: are you using negative viscosity or numerical precision too low."];
       Print["Will chop imaginary part of wavenumber to remedy"];
       Re[#],
       #
  ] &@ If[Re[#]>0, #,-#] &;

ThermalPotentials[subN_,rng\[Omega]_, opts:OptionsPattern[]]:=Module[{NSqrt\[Lambda]\[CurlyPhi]\[Omega],NSqrt\[Lambda]\[CurlyTheta]\[Omega],NV\[CurlyPhi],NV\[CurlyTheta],\[Lambda],V,subn = AddToNsub@subN},
  {{\[Lambda]["\[CurlyPhi]"],\[Lambda]["\[CurlyTheta]"],\[Lambda]["\[Psi]"]},{V["\[CurlyPhi]"],V["\[CurlyTheta]"]}} = LoadTWEigensystem[OptionValue@"ThermalSystem"];
  NSqrt\[Lambda]\[CurlyPhi]\[Omega][\[Omega]_] = OutgoingWaveNumber[Sqrt[\[Lambda]["\[CurlyPhi]"]]/.subn/.\[FormalW]-> \[Omega]];
  NSqrt\[Lambda]\[CurlyTheta]\[Omega][\[Omega]_] = OutgoingWaveNumber[Sqrt[\[Lambda]["\[CurlyTheta]"]]/.subn/.\[FormalW]-> \[Omega]];
  NV\[CurlyPhi][\[Omega]_] = V["\[CurlyPhi]"]/.subn/.\[FormalW]-> \[Omega];
  NV\[CurlyTheta][\[Omega]_] = V["\[CurlyTheta]"]/.subn/.\[FormalW]-> \[Omega];
  Return@{<|"Name" -> "Pressure dominated Wave \[CurlyPhi]", "rng\[Omega]"-> rng\[Omega], 
             "WaveNumbers"-> NSqrt\[Lambda]\[CurlyPhi]\[Omega]/@rng\[Omega], "WaveModes"->NV\[CurlyPhi]/@rng\[Omega], "parameters"->subn|>,
           <|"Name" -> "Thermal dominated Wave \[CurlyTheta]", "rng\[Omega]"-> rng\[Omega], 
             "WaveNumbers"->NSqrt\[Lambda]\[CurlyTheta]\[Omega]/@rng\[Omega], "WaveModes"-> NV\[CurlyTheta]/@rng\[Omega], "parameters"->subn|>
         }
];

ShearPotential[subN_,rng\[Omega]_, opts:OptionsPattern[]]:=Module[{NSqrt\[Lambda]\[Psi]\[Omega], \[Lambda]},
  \[Lambda]["\[Psi]"] = LoadTWEigensystem[OptionValue@"ThermalSystem"][[1,3]];
  NSqrt\[Lambda]\[Psi]\[Omega][\[Omega]_] = OutgoingWaveNumber[Sqrt[\[Lambda]["\[Psi]"]]/.subN]/.\[FormalW]-> \[Omega];
  Return@<|"Name" -> "Shear Wave", "rng\[Omega]"-> rng\[Omega], 
             "WaveNumbers"->NSqrt\[Lambda]\[Psi]\[Omega]/@rng\[Omega], "parameters"->AddToNsub@subN|> 
];

TWQ[TW_]:= TrueQ[Length@TW["rng\[Omega]"] == Length@TW["WaveNumbers"] == Length@TW["WaveVectors"]];

AddToNsub[subN_,opts:OptionsPattern[]]:= Module[
    {
     subn = subN,addparameters,subsols, 
     precision = Precision@subN, i, \[Lambda],eqResidues, maxeq
    },
  precision = If[Precision@subn < \[Infinity], Precision@subn, 15];
  subn = (#[[1]]-> (#[[2]]/.subn))&/@(subn/.Rule-> List);
  \[Lambda]["\[CurlyPhi]"] = LoadTWEigensystem[OptionValue@"ThermalSystem"][[1,1]]; 
  (*addparameters = Flatten@Reap[If[Head[#]===String, Sow[#]]&/@(AllParameters/.subn)][[2]];*)
  addparameters = Union@Cases[parameterEquations/.subn,_String,\[Infinity]];
  If[Length@addparameters == 0,
    eqResidues= Abs[#/#[[1]]]&/@(parameterEquations/.Equal[a_,b_]:>a-b)/.subn;
    maxeq = Extract[parameterEquations, Position[eqResidues,Max[eqResidues]]];
    If[Max[eqResidues]>0.005 || OptionValue@"Verbose", 
      Print["All parameters specified. The largest relative error was ",
          Round[ 1000 Max[eqResidues]]/10, "% from equation: ", maxeq
      ];
      Print["Here are the equations which need to be satisfied: \n",parameterEquations]
    ];  
     
    Return@subn
  ];
  subn = subn~Join~Thread[addparameters->ToExpression/@addparameters];
  subsols = Quiet@Solve[parameterEquations /.subn];
  
  If[Length@subsols ==0, subsols = NSolve[parameterEquations /.subn]];
  
  If[Length@subsols >1,
    negativeParameter[sol_] := Or@@(Negative/@ Transpose[sol/.Rule-> List][[2]]);
    i=1; While[negativeParameter[subsols[[i]]] && i<= Length[subsols],i++];
    If[i> Length[subsols], 
      Print["no viable solution found for remaining physical constants:", addparameters]; 
      Print["See proposed solutions:", subsols, "\n to the equations: ", parameterEquations/.subn];
    ];
    subn = subn/.subsols[[i]]/.subExpToPara;
    , 
    subn = 
      If[Length@subsols == 1, 
        If [!NumericQ@Total[(ToExpression/@addparameters)/.subsols[[1]]],
           Print["Not enough physical constants were given, see remaining equations: ", parameterEquations/.subn];
        ];
        subn = subn/.subsols[[1]]/.subExpToPara;
        subn~Join~{"LowFrequencyPressureWaveSpeed"-> \[FormalW]/Re@OutgoingWaveNumber[Sqrt[\[Lambda]["\[CurlyPhi]"]]/.subn/.\[FormalW]->10^(-precision/100)]/.\[FormalW]->N[10^(-precision/100)]}
      ,
        Print["Possibly too many physical constants were given because no solution was found for remaining physical constants:", addparameters]; 
        Print["See the equations that need to be satisfied: ", parameterEquations]; 
        subn/.subExpToPara
     ];
  ];
  subn = SetPrecision[subn,precision];
Return@subn
]; 




(* ::Code::Bold:: *)
(**)


(* ::Code::Initialization::Bold:: *)
End[];
EndPackage[];
