<!-- verified: cc916f6 2026-10-27 -->
# hyb_env_model
Hybrid GSF+PN waveform model after Honet+2025 (`biblio/Honet et al. - 2025 - ...pdf`), later extended with environmental effects. Current stage: reproduce the interpolated 0PA flux grid of Honet+ App. A (eqs. A9-A14), then build a vacuum hybrid model. Read the Honet+ PDF for equations and definitions.

## Commands
- Run notebook code headless with `wolframscript -file script.wls` (needs BHPToolkit `KerrGeodesics`` and `Teukolsky``). To evaluate the notebook's own cells: `Get` the `.nb`, pull `Cases[nb, Cell[b_, "Input", ___] :> b, Infinity]`, `ToExpression` them in order, keeping harness variables in a non-`Global`` context (cell 1 is `Clear["Global`*"]`). Shrink `nR`, `nA`, `lmax`, `modes` before the parallel cell for a quick check.

## Layout
- `Mathematica/FluxGrid_nonspinning.nb` - the only code. Energy fluxes at infinity/horizon for circular equatorial Kerr orbits on the Honet+ Chebyshev grid: rows `i` = rC nodes (i=0 at ISCO, i=nR-1 far field, r0 ~ 500 r*), columns `j` = χ nodes (j=0 at χ=+0.998 prograde, j=nA-1 retrograde). Returns `{rC, chi, r0^5 FI, r0^(15/2) FH}` per node (A13/A14). Section "Tests and debugging (ignore)" can be skipped.
- `Mathematica/benchmarks.md` - timings, precision scan, reference values, open questions. Append new measurements there.

## Where to change X
- Mode range or grid size: `nR`, `nA`, `lmax` in the definitions cell (`modes` derives from `lmax`).
- Precision: `precLadder = {40, 60, 90}` and `precGoal = 12`. `modeAdaptive` accepts a mode at the first prec where both fluxes pass `fluxOkQ`; a mode failing at all prec returns `{Indeterminate, Indeterminate, None}`, and the parallel cell prints `Counts` of prec used.
- Parallel run: last input cell of the main section, `ParallelMap` over (i, j, l, m) tasks, then `nodeFromModes` per node. Full grid at lmax=8 ≈ 36 min on 4 kernels.

## Conventions
- Edit `.nb` files through a headless front end (`UsingFrontEnd`, `NotebookOpen`, `NotebookWrite` with boxes from `FrontEnd`UndocumentedTestFEParserPacket`), then `SetOptions[nb, Visible -> Inherited]` before `NotebookSave` - keeps cache offsets and typesetting valid; text edits break both.
- Read `.nb` by grepping `Cell[` / `RowBox` input cells or dumping inputs via `wolframscript`; files carry cached outputs.

## Gotchas
- A notebook saved after `NotebookOpen[..., Visible -> False]` stores `Visible->False` and opens as an invisible window; reset it as above.
- Load the packages on the master kernel too before `DistributeDefinitions`; otherwise `KerrGeoOrbit` etc. are distributed as `Global`` shadows and everything comes back unevaluated.
- Low precision near χ=0.998 returns an exact `0` flux (`Precision = Infinity`), a silent wrong answer; `fluxOkQ` rejects it. prec 20 can hang for ~40 min per node there.
- Node (0,0) needs prec 90 at lmax=8 (fails at 60 with `SpinWeightedSpheroidalHarmonicS::maxterms`).
- Some modes emit `Infinite expression`/`Indeterminate expression` warnings during the Teukolsky mode computation; `Quiet` does not silence subkernel messages.
