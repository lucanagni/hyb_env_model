(* ::Package:: *)

(* SF side of the hybrid binding energy (Honet+2025 Sec. III C, eqs. 47-50).
   0PA only: the geodesic specific binding energy of a circular equatorial Kerr orbit, eq. 48,
     EhatGeo(x, chi) = -1 + (1 - 2 x f^(1/3)) / Sqrt[f (1 + chi x^(3/2) - 3 x f^(1/3))],  f = 1 - chi x^(3/2),
   with x = (M Omega)^(2/3), chi > 0 prograde, chi < 0 retrograde; E_0PA = nu M EhatGeo.
   The 1PA term nu^2 M E^SF_(1) (eq. 49, from the 1SF redshift z_(1), App. A2) is not implemented yet;
   its redshift grid would go in data/ next to this file.
   Usage:
     Get["/path/to/BindingEnergy.wl"]   (loads ../Fluxes/FluxInterpolation.wl for rISCO)
     EhatGeo[x, chi]                    closed form, also accepts symbolic arguments (series checks)
     xISCO[chi]                         ISCO value of x *)

With[{dir = DirectoryName[$InputFileName]},
  If[!MemberQ[$Packages, "FluxInterpolation`"], Get[FileNameJoin[{ParentDirectory[dir], "Fluxes", "FluxInterpolation.wl"}]]]
];

BeginPackage["BindingEnergy`", {"FluxInterpolation`"}];

EhatGeo::usage = "EhatGeo[x, chi]: geodesic specific binding energy E/(nu M) of a circular equatorial Kerr orbit (Honet eq. 48), x = (M Omega)^(2/3), signed spin chi.";
xISCO::usage = "xISCO[chi] = (rISCO[chi]^(3/2) + chi)^(-2/3).";

Begin["`Private`"];

EhatGeo[x_, chi_] := With[{f = 1 - chi x^(3/2)},
  -1 + (1 - 2 x f^(1/3))/Sqrt[f (1 + chi x^(3/2) - 3 x f^(1/3))]
];

xISCO[chi_?NumericQ] := (rISCO[chi]^(3/2) + chi)^(-2/3);

End[];

EndPackage[];
