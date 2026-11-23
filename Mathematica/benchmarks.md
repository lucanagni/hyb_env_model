# FluxesGrid benchmarks (Honet+ App. A, f_inf grid; nR=18, nA=27 -> 486 nodes)
Setup: Mathematica 14.0, BHPToolkit KerrGeodesics+Teukolsky, $ProcessorCount=4 (nproc=8) -> LaunchKernels[] = 4 kernels.
Node: r0=2*rstar/(rC+1), chi=cos((2j+1)pi/(2nA)), orbit KerrGeoOrbit[a,r0,0,x]; modes l=2..lmax, m=1..l; F=2*sum flux. Notebook line 46-47 = lmax,prec.
Sparse grid S = i{0,4,9,14,17} x j{0,7,13,20,26} (25 nodes). Scripts in scratchpad: bench.wls, par.wls, par8.wls, adapt.wls.

## Precision (serial, lmax=4, S)
- prec30: 1.7 s/node, 6/25 fail (all i=17; (0,0)); surviving nodes keep 5-21 digits. Cause: r0>~2500 (i=17 -> r0~520*rstar) and near-extremal chi.
- prec60: 4.0 s/node, 0 fail, min 28 digits. Full grid ~33 min serial. prec100 gives identical values.
- i=17 nodes: f_inf=6.37-6.40 (limit 32/5=6.4 OK).

## Parallel (ParallelTable, Method->"FinestGrained", DistributeDefinitions[nR,nA,lmax,prec,node,fluxAtNode], ParallelEvaluate Needs)
- lmax=4 prec60 S: 100 s serial -> 33.2 s (3.0x), identical results. Full ~11 min.
- lmax=8 prec60 S: 160 s (4.8x lmax=4; modes 9->35 = 3.9x). Full ~52 min.
  - node (0,0) (ISCO, chi=+0.998) -> Indeterminate: FindRoot::cvmit, SpinWeightedSpheroidalEigenvalue::findroot, SpinWeightedSpheroidalHarmonicS::maxterms (Leaver 200 terms). Returns Indeterminate, not $Failed -> test with NumberQ/Precision, not MemberQ[$Failed].
  - min precision ~16 digits (28 at lmax=4).
  - Quiet does not silence subkernel messages.
- lmax=8 prec75 S: 223 s (+39% vs prec60), 0 bad, min 28 digits. Full ~72 min. (0,0) alone serial: 36 s; f_inf(0,0)=0.25213148.
- lmax=8 adaptive (try prec60, if non-finite or Precision<=8 rerun prec75) S: 221 s, 0 bad, matches prec75 to 8e-22 rel. No gain on S: only 1 failing node, its retry is serial tail. Not tested on full grid.

## Reference values f_inf (lmax=8, S rows i=0,... cols j=0,7,13,20,26)
i=0: 0.25213, 4.82083, 7.27824, 8.55102, 8.86837; i=17 row ~6.3655-6.3953.

## Open / untested
- Adaptive fallback on full grid; schedule hard nodes (i=0 row, j=0 and j=26 cols) first to hide tail.
- Lower lmax for hardest nodes (high-l flux contribution likely small).
- LaunchKernels[8]; Leaver term/Method option for spheroidal harmonics.
- Paper uses l' up to 50; near-extremal nodes will worsen with lmax.

## Per-mode precision scan (lmax=8, S, every mode at prec 20/30/40/60 vs prec90 reference; times = summed kernel-s)
- Prec needed for 1e-12 on node totals: 30 -> 16 nodes; 40 -> (0,7),(14,26),(17,0/7/20); 60 -> (9,0),(17,26); 90 -> (0,0),(4,0). Uniform prec60 wrong at (0,0),(4,0).
- prec20 can hang: (0,0),(4,0) ~2400 s each. Start ladder at >=30.
- Trap: at prec 20-40 near chi=+0.998 some modes return exact 0 (Precision=Infinity), up to 6% of f_inf. Otherwise gate "finite, Precision<Infinity, >=12" never accepted a node total off by >1e-12.
- Cost: uniform p90 1111 s; p60 638 s; node-level ladder 40/60/90 652 s (max err 3e-19); mode-level ladder 40/60/90 427 s (max err 5e-15, all 25 nodes ok). Oracle 400 s.
- Notebook now uses mode-level ladder 40/60/90, gate 12. fluxAtNode[0,0] serial 21 s, f_inf=0.252131478538500872 (matches prec75). Full grid not yet run.

## Chebyshev interpolant (FluxInterpolation.wl, nR=18, nA=27, lmax=8)
- Reproduces nodes to 4e-12 rel. fInf(rC=-1) = 6.3993-6.3998 (Newtonian 6.4). FluxInf[10, 0] = 6.15035e-5, FluxHor[10, 0] = 1.259e-8.
- Coefficient tail (max |c| ~ 6-7): last 3 in rC ~ 6e-5, 4e-5, 2e-5 (f_inf), 1.5e-4..4e-5 (f_H); in chi ~ 1.8e-3, 1.2e-3, 6e-4 (f_inf), 1.4e-3..4.5e-4 (f_H). Chi direction converges slowly -> interpolation error ~1e-4 rel; Honet+ A11 log-spin compactification targets this. Not yet checked against an off-grid Teukolsky evaluation.
- Evaluation ~27 us per call.
- Nested-subgrid test (6x9 subgrid i=3i'+1, j=3j'+1, evaluated on all 18x27 nodes): max rel err 0.82 at chi=0.998 column (outside the coarse grid's outermost node 0.985), 5% elsewhere; coarse coefficient tail only ~1.5e-2. f_inf is steep as chi->1 (ISCO has (1-chi)^(1/3) branch point), so the full-grid error near chi=0.998 may exceed the 1e-4 coefficient estimate. Argues for the A12 log(1-chi) variable.

## ModesGrid amplitude grid (nR=18, nA=27, lmax=8, nExtra=2 -> l' up to 10), 2026-11-21
- Gotcha: m=3, l'=10 at nodes (i=0..3, j=0) (ISCO side, chi=0.998) fails the gate at prec 40/60/90 -> modeAdaptive Indeterminate -> mu . vec makes every c_l3 of the node Indeterminate. Fixed with precLadder = {40, 60, 90, 120}.
- Digits of A_inf(l'=10, m=3) at prec 90: i=0 3.2, i=1 3.5, i=2 5.4, i=3 7.6 (l'=9: 14-17.5). At prec 60, l'>=8 gives Indeterminate/ComplexInfinity/exact 0; at prec 40 already from l'=7.
- i=0, l'=10: prec 120 -> 32 digits, 150 -> 62, 200 -> 111, values agree. ~2 s/mode.
- Loss is in MST (TeukolskyRadial In/Up amplitudes and R(r0) keep 5-6 of 90 digits); lambda and nu fine (nu 73 digits). Same node m=2/m=4 l'=10 keep 62-76 digits; m=3 at i=4 keeps 70. Why m=3 specifically: not established.
- m=3 rerun with the new ladder: 514 s on 4 kernels, 0 bad nodes; other 482 nodes unchanged to ~12 digits; c_l3(0,0) carries 16-20 digits.

## SF vs PN mode conventions (HybridModes.nb Hybrid chapter), 2026-11-23
- Amplitudes: |h_lm|/nu = 2|C_lm|/(m Omega)^2 with Omega = x^(3/2) (Honet eq. 53); PNpedia h_lm is already the strain mode (no 2/(m Omega) factor). 0PA vs PN (nu->0) agree to 1e-6..2e-4 at x=0.01 for (2,2),(2,1),(3,3),(3,2),(4,4), a=0,+-0.5; ~1% at x=0.1 ((4,4) 15%).
- Phase sign: Arg[(-1)^(m+1) C_lm] = Arg[h^PN_lm] - m c(x) for all 35 modes l=2..8, x=0.005,0.02, a=0,+-0.6 (offset exactly 0 for odd m, pi for even m). (-1)^m: PN phase follows body 1 (n=(y1-y2)/r12, delta>0 -> primary), SF phi_p follows the secondary. Global -1: polarisation triad or psi4->h sign (Honet A24 has +2/(m Omega)^2); not distinguishable numerically. Honet+2025 does not state either. hybSF uses Arg[(-1)^(m+1) cSF].
- c(x) (0.0047 at x=0.005, 0.026 at 0.02, 0.139 at 0.1) is a phase-origin shift proportional to m; it cancels only in dpsi_lm = dphi_lm - dphi_22 (Honet eqs. 15, 56), so HybridInfAmplPhase hybridises dpsi. The (-1)^m part also cancels in dpsi; the global -1 does not (offset pi/2 - pi/m).
- Residual after removing sign and shift: <=1e-4 rad for l<=5 at x=0.005; ~2e-3 for l=6..8 (all l=7,8; l=6 odd m), growing ~x^(3/2). Probably PN truncation of high-l modes; not checked.
- Open: which convention Honet flipped (ask authors).
