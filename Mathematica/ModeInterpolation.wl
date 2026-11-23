(* ::Package:: *)

(* Chebyshev interpolant of the spherical-basis mode amplitudes cbar_lm produced by ModesGrid.nb
   (Honet+2025 App. A, eqs. A27-A29), on the same grid as FluxInterpolation.wl:
   rC_i = Cos[(2i+1)Pi/(2nR)], chi_j = Cos[(2j+1)Pi/(2nA)], r0 = 2 r*/(rC+1), r* = ISCO radius.
   Stored: cbar_lm = C_lm * pref, pref = r0^(l/2)/Omega^2 (l+m even), r0^((l+1)/2)/Omega^2 (l+m odd),
   Omega = 1/(chi + r0^(3/2)); C_lm = spherical-basis Teukolsky infinity amplitude (m >= 1 only).
   Usage:
     Get["/path/to/ModeInterpolation.wl"]   (loads FluxInterpolation.wl and data/m*/amplGrid_*.m)
     cBar[l, m, rC, chi]                    interpolated compactified cbar_lm, rC, chi in [-1,1]
     ModeAmp[l, m, r0, chi]                 physical C_lm (units m1 = 1), r0 >= r*(chi)
     ModePrefactor[l, m, r0, chi]           pref above
     ModeList[]                             loaded {l, m} pairs
     LoadModeGrid[dir]                      reload from another data directory *)

With[{dir = DirectoryName[$InputFileName]},
  If[!MemberQ[$Packages, "FluxInterpolation`"], Get[FileNameJoin[{dir, "FluxInterpolation.wl"}]]]
];

BeginPackage["ModeInterpolation`", {"FluxInterpolation`"}];

LoadModeGrid::usage = "LoadModeGrid[dir] loads every m*/amplGrid_nR*_nA*_l*m*.m under dir and builds the Chebyshev coefficients; returns the loaded {l, m} pairs.";
ModeList::usage = "ModeList[] gives the sorted list of loaded {l, m} pairs.";
cBar::usage = "cBar[l, m, rC, chi]: interpolated compactified amplitude cbar_lm = C_lm * pref on the Chebyshev grid.";
ModeAmp::usage = "ModeAmp[l, m, r0, chi]: spherical-basis infinity amplitude C_lm for a circular equatorial orbit at radius r0, signed spin chi.";
ModePrefactor::usage = "ModePrefactor[l, m, r0, chi]: r0^(l/2)/Omega^2 (l+m even) or r0^((l+1)/2)/Omega^2 (l+m odd), Omega = 1/(chi + r0^(3/2)).";
ModeGridData::usage = "ModeGridData[{l, m}]: the raw {nR, nA, 3} grid {rC, chi, cbar_lm}.";
ModeChebCoeffs::usage = "ModeChebCoeffs[{l, m}]: complex nR x nA Chebyshev coefficient matrix of cbar_lm.";
cBar::nomode = "Mode {l, m} = `1` not loaded; available: ModeList[].";

Begin["`Private`"];

(* 2D DCT-II on Chebyshev-Gauss nodes: f(x,y) = Sum c_kl T_k(x) T_l(y) (copy of FluxInterpolation`Private`) *)
chebCoeffs[f_?MatrixQ] := Module[{nR, nA, tR, tA},
  {nR, nA} = Dimensions[f];
  tR = Table[(2 - KroneckerDelta[k, 0])/nR Cos[k (2 i + 1) Pi/(2 nR)], {k, 0, nR - 1}, {i, 0, nR - 1}];
  tA = Table[(2 - KroneckerDelta[l, 0])/nA Cos[l (2 j + 1) Pi/(2 nA)], {l, 0, nA - 1}, {j, 0, nA - 1}];
  N[tR . f . Transpose[tA]]
];

chebT[n_, x_] := Cos[Range[0, n - 1] ArcCos[x]];

chebEval[c_, x_?NumericQ, y_?NumericQ] /; -1 <= x <= 1 && -1 <= y <= 1 :=
  chebT[Length[c], x] . c . chebT[Length[First[c]], y];

modeOfFile[f_String] := First[StringCases[FileBaseName[f],
  "_l" ~~ l : DigitCharacter .. ~~ "m" ~~ m : DigitCharacter .. ~~ EndOfString :> {ToExpression[l], ToExpression[m]}]];

LoadModeGrid[dir_String] := Module[{files},
  files = FileNames["amplGrid_nR*_nA*_l*m*.m", Select[FileNames["m*", dir], DirectoryQ]];
  Clear[ModeGridData, ModeChebCoeffs];
  Do[With[{lm = modeOfFile[f], g = Get[f]},
      ModeGridData[lm] = g;
      ModeChebCoeffs[lm] = chebCoeffs[N[g[[All, All, 3]]]]],
    {f, files}];
  ModeList[]
];

ModeList[] := Sort[Cases[DownValues[ModeChebCoeffs][[All, 1]], Verbatim[HoldPattern][_[{l_Integer, m_Integer}]] :> {l, m}]];

cBar[l_Integer, m_Integer, rC_?NumericQ, chi_?NumericQ] :=
  If[MatrixQ[ModeChebCoeffs[{l, m}]],
    chebEval[ModeChebCoeffs[{l, m}], rC, chi],
    Message[cBar::nomode, {l, m}]; $Failed];

ModePrefactor[l_Integer, m_Integer, r0_?NumericQ, chi_?NumericQ] :=
  If[EvenQ[l + m], r0^(l/2), r0^((l + 1)/2)] (chi + r0^(3/2))^2;

ModeAmp[l_Integer, m_Integer, r0_?NumericQ, chi_?NumericQ] :=
  With[{c = cBar[l, m, rCOf[r0, chi], chi]},
    If[c === $Failed, $Failed, c/ModePrefactor[l, m, r0, chi]]];

End[];

LoadModeGrid[FileNameJoin[{DirectoryName[$InputFileName], "data"}]];

EndPackage[];
