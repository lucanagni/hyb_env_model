# FluxGrid_nonspinning benchmarks (Honet+ App. A, f_inf grid; nR=18, nA=27 -> 486 nodes)
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
- Notebook currently: lmax=8, prec=60, no fallback (node (0,0) still Indeterminate).
