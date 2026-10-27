(* ::Package:: *)

(* Chebyshev interpolant of the 0PA flux grid produced by FluxGrid_nonspinning.nb
   (Honet+2025 App. A, eqs. A9-A14).
   Grid: rC_i = Cos[(2i+1)Pi/(2nR)], chi_j = Cos[(2j+1)Pi/(2nA)] (Chebyshev-Gauss nodes),
   r0 = 2 r*/(rC+1), r* = ISCO radius; data rows {rC, chi, r0^5 FI, r0^(15/2) FH}.
   Usage:
     Get["/path/to/FluxInterpolation.wl"]   (loads data/fluxGrid_nR18_nA27_lmax8.m by default)
     fInf[rC, chi], fHor[rC, chi]           compactified A13/A14 functions, rC, chi in [-1,1]
     FluxInf[r0, chi], FluxHor[r0, chi]     physical fluxes (units m1 = 1), r0 >= r*(chi)
     LoadFluxGrid[file]                     switch to another grid file *)

BeginPackage["FluxInterpolation`"];

LoadFluxGrid::usage = "LoadFluxGrid[file] loads an {nR, nA, 4} flux grid and builds the Chebyshev coefficients.";
fInf::usage = "fInf[rC, chi] = r0^5 FI interpolated on the compactified grid.";
fHor::usage = "fHor[rC, chi] = r0^(15/2) FH interpolated on the compactified grid.";
FluxInf::usage = "FluxInf[r0, chi]: energy flux to infinity for a circular equatorial orbit at radius r0, signed spin chi.";
FluxHor::usage = "FluxHor[r0, chi]: energy flux through the horizon.";
rISCO::usage = "rISCO[chi]: Boyer-Lindquist ISCO radius; chi > 0 prograde, chi < 0 retrograde.";
rCOf::usage = "rCOf[r0, chi] = 2 rISCO[chi]/r0 - 1.";
FluxGridData::usage = "FluxGridData: the raw {nR, nA, 4} grid.";
FluxChebCoeffs::usage = "FluxChebCoeffs: {cInf, cHor}, nR x nA Chebyshev coefficient matrices.";

Begin["`Private`"];

(* 2D DCT-II on Chebyshev-Gauss nodes: f(x,y) = Sum c_kl T_k(x) T_l(y) *)
chebCoeffs[f_?MatrixQ] := Module[{nR, nA, tR, tA},
  {nR, nA} = Dimensions[f];
  tR = Table[(2 - KroneckerDelta[k, 0])/nR Cos[k (2 i + 1) Pi/(2 nR)], {k, 0, nR - 1}, {i, 0, nR - 1}];
  tA = Table[(2 - KroneckerDelta[l, 0])/nA Cos[l (2 j + 1) Pi/(2 nA)], {l, 0, nA - 1}, {j, 0, nA - 1}];
  N[tR . f . Transpose[tA]]
];

chebT[n_, x_] := Cos[Range[0, n - 1] ArcCos[x]];

chebEval[c_, x_?NumericQ, y_?NumericQ] /; -1 <= x <= 1 && -1 <= y <= 1 :=
  chebT[Length[c], x] . c . chebT[Length[First[c]], y];

LoadFluxGrid[file_String] := (
  FluxGridData = Get[file];
  FluxChebCoeffs = chebCoeffs /@ {FluxGridData[[All, All, 3]], FluxGridData[[All, All, 4]]};
  Dimensions[FluxGridData]
);

fInf[rC_?NumericQ, chi_?NumericQ] := chebEval[FluxChebCoeffs[[1]], rC, chi];
fHor[rC_?NumericQ, chi_?NumericQ] := chebEval[FluxChebCoeffs[[2]], rC, chi];

rISCO[chi_?NumericQ] := Module[{a = Abs[chi], z1, z2},
  z1 = 1 + (1 - a^2)^(1/3) ((1 + a)^(1/3) + (1 - a)^(1/3));
  z2 = Sqrt[3 a^2 + z1^2];
  3 + z2 - Sign[chi] Sqrt[(3 - z1) (3 + z1 + 2 z2)]
];

rCOf[r0_, chi_] := 2 rISCO[chi]/r0 - 1;

FluxInf[r0_?NumericQ, chi_?NumericQ] := fInf[rCOf[r0, chi], chi]/r0^5;
FluxHor[r0_?NumericQ, chi_?NumericQ] := fHor[rCOf[r0, chi], chi]/r0^(15/2);

End[];

LoadFluxGrid[FileNameJoin[{DirectoryName[$InputFileName], "data", "fluxGrid_nR18_nA27_lmax8.m"}]];

EndPackage[];
