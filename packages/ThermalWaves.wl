(* ::Package:: *)

(* ::Code::Bold:: *)
(**)


BeginPackage["ThermalWaves`"];

ClearAll @@ Names["ThermalWaves`*"];

$Path = $Path~Join~{NotebookDirectory[],NotebookDirectory[]<>"packages/"};
$Path = $Path~Join~{NotebookDirectory[]<>"../packages/"};
Get["ThermalPotentials.wl"];
(*Needs["NumericalCalculus`"];*)

ThermalWaves::usage = "ThermalWaves[subN,waveDirection, rng\[Omega]] = {TW\[CurlyPhi], TW\[CurlyTheta], TW\[Psi]}, where TW is a thermowave. See section \"Reflection from a halfspace\" in the document TVE-equations.pdf for more details.  ";
DisplacementTW::usage = "DisplacementTW[thermalWave, vecX]";
\[Theta]TW::usage = "\[Theta]TW[thermalWave, vecX]";
StressTW::usage = "StressTW[thermalWave, vecX] gives the Cauchy stress at vecX from the origin of the plane wave.";
EnergyFluxTW::usage = "EnergyFluxTW[thermalWave]";
(*ScatteredTW::usage = "ScatteredTW[TWin,TWRs,TWTs,opts:OptionsPattern[]]";*)
TWQ::usage = "TWQ[TW] returns true if TW is a ThermalWave in the right format";
PlotTW::usage = "PlotThermalWave[thermalWave] plots the wave given in the association";
PlusTW::usage = "PlusTW[TW1_,TW2_] returns a thermal wave with both wavevectors and wavenumbers summed together"


(* ::Code::Bold:: *)
(**)


Begin["`Private`"]

{{\[Lambda][\[CurlyPhi]],\[Lambda][\[CurlyTheta]],\[Lambda][\[Psi]]},{V[\[CurlyPhi]],V[\[CurlyTheta]]}} = TWEigensystem;
(* where \[CurlyPhi] \[Equal] \[ExponentialE]^(\[PlusMinus] \[ImaginaryI] Sqrt[\[Lambda][\[CurlyPhi]]]x.n), \[CurlyTheta] \[Equal] \[ExponentialE]^(\[PlusMinus] \[ImaginaryI] Sqrt[\[Lambda][\[CurlyTheta]]]x.n), \[VerticalSeparator]n\[VerticalSeparator] \[Equal] 1, 
and 
{\[Phi], \[Theta]} \[Equal]  \[CurlyPhi] V[\[CurlyPhi]] + \[CurlyTheta] V[\[CurlyTheta]] *)
TWQ[TW_]:= TrueQ[Length@TW["rng\[Omega]"] == Length@TW["WaveNumbers"] == Length@TW["WaveVectors"]];

vec3D = If[Length@# == 2, #~Join~{0.}, #, #]&;
Unitvec3D = If[N@Norm@# == 0., 1., 1./Sqrt[#.#], 1./Sqrt[#.#]] vec3D[#]&;
OutgoingWaveNumber = 
  If[ Im@#<0, 
       Print["Wave amplitude increases in time: are you using negative viscosity or numerical precision too low?"];
       Print["To remedy, we will return only the real part of wavenumber."];
       Re[#],
       #
  ] &@ If[Re[#]>0, #,-#] &;

ThermalWaves[subN_,waveDirection_,rng\[Omega]_]:=Module[{NSqrt\[Lambda]\[CurlyPhi]\[Omega],NSqrt\[Lambda]\[CurlyTheta]\[Omega],NSqrt\[Lambda]\[Psi]\[Omega],NV\[CurlyPhi],NV\[CurlyTheta], vecNs = Unitvec3D/@ArrayReshape[waveDirection,{3,3},"Periodic"]},
  NSqrt\[Lambda]\[CurlyPhi]\[Omega][\[Omega]_] = OutgoingWaveNumber[Sqrt[\[Lambda][\[CurlyPhi]]]/.subN]/.\[FormalW]-> \[Omega];
  NSqrt\[Lambda]\[CurlyTheta]\[Omega][\[Omega]_] = OutgoingWaveNumber[Sqrt[\[Lambda][\[CurlyTheta]]]/.subN]/.\[FormalW]-> \[Omega];
  NSqrt\[Lambda]\[Psi]\[Omega][\[Omega]_] = OutgoingWaveNumber[Sqrt[\[Lambda][\[Psi]]]/.subN]/.\[FormalW]-> \[Omega];
  NV\[CurlyPhi][\[Omega]_] = V[\[CurlyPhi]]/.subN/.\[FormalW]-> \[Omega];
  NV\[CurlyTheta][\[Omega]_] = V[\[CurlyTheta]]/.subN/.\[FormalW]-> \[Omega];

  Return@{<|"Name" -> "Pressure dominated Wave \[CurlyPhi]", "rng\[Omega]"-> rng\[Omega], "WaveNumbers"->(vecNs[[1]] NSqrt\[Lambda]\[CurlyPhi]\[Omega][#] &/@rng\[Omega]), "WaveVectors"->({#/Norm@#,{0.,0.,0.}} &/@(NV\[CurlyPhi]/@rng\[Omega])), "parameters"->AddToNsub@subN |>,
           <|"Name" -> "Thermal dominated Wave \[CurlyTheta]", "rng\[Omega]"-> rng\[Omega], "WaveNumbers"->(vecNs[[2]] NSqrt\[Lambda]\[CurlyTheta]\[Omega][#] &/@rng\[Omega]), "WaveVectors"->({#/Norm@#,{0.,0.,0.}} &/@(NV\[CurlyTheta]/@rng\[Omega])), "parameters"->AddToNsub@subN |>,
           <|"Name" -> "Shear Wave", "rng\[Omega]"-> rng\[Omega], "WaveNumbers"->(vecNs[[3]] NSqrt\[Lambda]\[Psi]\[Omega][#] &/@rng\[Omega]), "WaveVectors"->ConstantArray[{{0.,0.},{1.,0.,0.}},Length@rng\[Omega]], "parameters"->AddToNsub@subN|> 
         }
];



ClearAll[PlusTW]
PlusTW[TW1_,TW2_]:= Module[ {TW = TW1+TW2},
  If[!TrueQ[TW1["parameters"] == TW2["parameters"]], Print["Combining thermal waves with different mateiral parameters does not make sense. Averaged vlaues will be used."]];
  If[!TrueQ[TW1["rng\[Omega]"] == TW2["rng\[Omega]"]], Print["Combining thermal waves with different frequencies does not make sense. Averaged vlaues will be used."]];  
  TW["Name"] = ToString@TW["Name"];
  TW["rng\[Omega]"] =TW["rng\[Omega]"]/2;
  TW["parameters"]=(#[[1]]/2->#[[2]]/2) &/@(TW["parameters"]/.Rule-> List);
TW]

\[Theta]TW[thermalWave_, vectorX_]:=
Module[{vecX = vec3D@vectorX},
   E^(I vecX.#[[1]]) #[[2,1,2]]  &/@Transpose[{thermalWave["WaveNumbers"], thermalWave["WaveVectors"]}]
];
DisplacementTW[thermalWave_, vectorX_]:=
Module[{vecX = vec3D@vectorX},
   I E^(I vecX.#[[1]])(#[[1]] #[[2,1,1]] + #[[1]]\[Cross]#[[2,2]]) &/@Transpose[{thermalWave["WaveNumbers"], thermalWave["WaveVectors"]}]
];

StressTW[thermalWave_, vectorX_]:=
  Module[{vecX = vec3D@vectorX, subN = thermalWave["parameters"]},
    E^(I vecX.#[[2]]) (
       - IdentityMatrix[3]( ("\[Lambda]" - I #[[1]] "\[Eta]\[Lambda]") #[[3,1,1]] - "pT"#[[3,1,2]] /.subN) 
      - ("\[Mu]" - I #[[1]] "\[Eta]\[Mu]"/.subN) ( 2 KroneckerProduct[#[[2]],#[[2]]]  #[[3,1,1]] 
      + KroneckerProduct[#[[2]]\[Cross]#[[3,2]],#[[2]]] + Transpose@KroneckerProduct[#[[2]]\[Cross]#[[3,2]],#[[2]]]  )
    )&/@Transpose[{thermalWave["rng\[Omega]"], thermalWave["WaveNumbers"], thermalWave["WaveVectors"]}] 
  ];

EnergyFluxTW[thermalWave_,vecX_:{0.,0.,0.}]:= Module[{v = -I thermalWave["rng\[Omega]"] DisplacementTW[thermalWave,vecX]},
     -0.5Re@MapThread[Dot,{v\[Conjugate],#}] &@ StressTW[thermalWave, {0.,0.,0.}]
  ];

colorF[{x_,y_}]:=Blend[{Blend[{White,Green},0.5],Blend[{Blue,Red},x]},y] ;

Options[PlotTW]= {"ColorStyle"->"EnergyFlux", "ReferenceTW" -> {}};
PlotTW[thermalWave_,opts:OptionsPattern[]]:=
Module[{rngv,rngDecay,rngColors,rng\[Omega] = thermalWave["rng\[Omega]"], waveNumber, EI\[Theta], EI\[Phi], colorData, TW\[Phi] = thermalWave, legend,pos,g, boolRefTW},
    waveNumber = OutgoingWaveNumber[Sqrt[#.#]] &/@thermalWave["WaveNumbers"] ;

    boolRefTW = Quiet@NumericQ[Plus@@(("ReferenceTW"/.opts)["WaveNumbers"][[1]])];
    If[ boolRefTW &&!TrueQ[("ColorStyle"/.opts)=="ThermalRatio"],
      colorData = Norm/@EnergyFluxTW[thermalWave] / Norm/@EnergyFluxTW["ReferenceTW"/.opts];
      legendlabel = "Relative energy flux";
      listlegend = {"0%","50%","100%"};
      rngDecay= ConstantArray[1.,Length@rng\[Omega]];
    ,
      TW\[Phi]["WaveVectors"] = {{#[[1,1]],0.},#[[2]]}&/@ TW\[Phi]["WaveVectors"];
      EI\[Phi] = EnergyFluxTW[TW\[Phi]];
      EI\[Theta] = Norm/@(EnergyFluxTW[thermalWave] - EI\[Phi]);
      EI\[Phi] = Norm/@EI\[Phi];
      colorData = If[#[[1]] == 0., 1., #[[2]]/(#[[1]] + #[[2]])]& /@ Transpose[{EI\[Phi],EI\[Theta]}];
      legendlabel = "\[Theta] energy flux";
      listlegend = {"0%","50%","100%", "Decay rate"};
      rngDecay=Clip[Rescale[-Im@waveNumber,{-4000,0},{0.,1}],{0.,1}];    
    ];
  
  legend = PointLegend[{colorF[{0.,1.}],colorF[{0.5,1.}],colorF[{1.,1.}],colorF[{1.,0.}]},listlegend, LegendLabel->legendlabel];

  
  rngColors=colorF/@Transpose[{colorData,rngDecay}];
  rngv= rng\[Omega]/Re@waveNumber;
  pos= {{0,0}}~Join~Transpose[{rng\[Omega],rngv}];
  
  colorData={White}~Join~rngColors;
  gopts = {Axes->True,AspectRatio->1,AxesLabel-> {"\[FormalW]","\!\(\*SubscriptBox[\(V\), \(1\)]\)"}};
  If[!TrueQ[thermalWave["Name"] == Missing["KeyAbsent","Name"]], gopts = gopts~Join~{PlotLabel->thermalWave["Name"] } ];
  g=Graphics[{PointSize [Large], Point[pos,VertexColors->colorData]},gopts];
   
  Return[Legended[g,legend]]
];

(*eqns = LoadTWScattering;*)

End[];
EndPackage[];













