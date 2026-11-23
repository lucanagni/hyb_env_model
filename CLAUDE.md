<!-- verified: 039bb54 2026-11-23 -->
# hyb_env_model
Hybrid GSF+PN waveform model after Honet+2025 (`biblio/Honet et al. - 2025 - ...pdf`), later extended with environmental effects. Current stage: vacuum hybrid model; 0PA flux and mode-amplitude grids are done, PN hybridisation of fluxes done, of amplitudes in progress. Read the Honet+ PDF for equations and definitions.

## Commands
- Run notebook code headless with `wolframscript -file script.wls` (needs BHPToolkit `KerrGeodesics`` and `Teukolsky``). To evaluate the notebook's own cells: `Get` the `.nb`, pull `Cases[nb, Cell[b_, "Input", ___] :> b, Infinity]`, `ToExpression` them in order, keeping harness variables in a non-`Global`` context (cell 1 is `Clear["Global`*"]`). Shrink `nR`, `nA`, `lmax`, `modes` before the parallel cell for a quick check.

## Layout
- `Mathematica/FluxesGrid.nb` - flux computation. Energy fluxes at infinity/horizon for circular equatorial Kerr orbits on the Honet+ Chebyshev grid: rows `i` = rC nodes (i=0 at ISCO, i=nR-1 far field, r0 ~ 500 r*), columns `j` = χ nodes (j=0 at χ=+0.998 prograde, j=nA-1 retrograde). Returns `{rC, chi, r0^5 FI, r0^(15/2) FH}` per node (A13/A14). Section "Tests and debugging (ignore)" can be skipped.
- `Mathematica/data/fluxGrid_nR18_nA27_lmax8.m` - full grid (18x27x4, lmax=8) extracted from the notebook's cached `data` output; `Get` returns the array.
- `Mathematica/FluxInterpolation.wl` - package: loads the grid, builds 2D Chebyshev (DCT-II) coefficients; `fInf`/`fHor[rC, chi]`, `FluxInf`/`FluxHor[r0, chi]`, `rISCO[chi]` (closed form, no KerrGeodesics needed), `LoadFluxGrid[file]` for other grids.
- `Mathematica/ModesGrid.nb` - spherical-basis mode amplitudes c̄_ℓm (Honet A27-A29) on the same Chebyshev grid: spheroidal C^up up to ℓ' = lmax + `nExtra` (=2), change of basis via `muMatrix`, PN prefactor; section "Amplitude grid" runs `ParallelMap` over nodes per m and writes the files below. Retrograde orbits use a<0 with x=1 (Honet frame, continuous at χ=0; x=-1 would give the complex conjugate).
- `Mathematica/data/m<m>/amplGrid_nR18_nA27_l<l>m<m>.m` - one file per (ℓ, m), nR x nA array of `{rC, chi, c̄_ℓm}`.
- `Mathematica/ModeInterpolation.wl` - package (loads `FluxInterpolation.wl` first, reuses its `rISCO`/`rCOf`): Chebyshev interpolant of the `data/m<m>/amplGrid_*.m` grids, m ≥ 1 only; `cBar[l, m, rC, chi]` (stored c̄_ℓm), `ModeAmp[l, m, r0, chi]` = c̄_ℓm / `ModePrefactor` (physical spherical C_ℓm), `ModeList[]`, `LoadModeGrid[dir]`.
- `Mathematica/HybridEnFluxes.nb` - PN ∞-flux assembly (nonspinning 4.5PN + SO/SS/SSS, sources in `biblio/`), SF+PN hybrid ∞-flux (Honet eq. 43) and horizon flux (Saketh+ PN, Honet eq. 44). SF fluxes use `FluxInf`/`FluxHor` from `FluxInterpolation.wl` (loaded by absolute path) with r0(x) from eq. A15.
- `Mathematica/HybridModes.nb` - PN mode amplitudes h_ℓm (ℓ = 2..8, m ≥ 1, 3.5PN; 4PN nonspinning h22) from PNpedia (`~/repos/PNpedia`), split into `hNS`/`hSO`/`hSS`/`hSSS`, total `hPN[l, m][x, ν, χ]` with Honet's spin map; Checks section compares with the PNpedia files; hand-written Hybrid chapter (`hybSF`/`hybPN`/`hybPNSF`, `HybridInfAmplPhase`, `PlotWF`) for the SF+PN mode hybrid (Honet eqs. 53-56). The PN definitions come from `Mathematica/tools/make_PNModes_nb.wls`, which writes a fresh notebook (default `PNModes.nb`); regenerate into a separate file and copy the definition cells over, because writing to `HybridModes.nb` would wipe the Hybrid chapter.
- `Mathematica/tuorials/ChebyshevTutorial.nb` - tutorial: Chebyshev interpolation theory + exercises, rebuilds `FluxInterpolation.wl` from the grid file.
- `Mathematica/tuorials/SpheroidalToSphericalTutorial.nb` - tutorial: spheroidal→spherical change of basis for Teukolsky amplitudes (Honet A24–A31), exercises + collapsed solutions. Defines `muCol`/`muMatrix` (μ from the package's `"SphericalExpansion"` coefficients). Sections 6–7 run from the precomputed amplitudes via a toolbox cell (`amps`, `toSpherical`, `cTrunc`, `errTrunc`); section 8 is the grid-rerun recipe.
- `Mathematica/tuorials/data/sphdAmps_examples.m` - Teukolsky ∞-amplitudes C^up for 6 circular equatorial orbits (keys like `"a0.9_isco"`), m = 2,4,6,8, l' = max(2,m)..14, ≥12 digits; `Get` returns `<|key -> <|"a","r0","x","modes" -> <|m -> <|"omega","gamma","lp","Cup","prec"|>|>|>|>`.
- `notes/spheroidal_to_spherical.tex` (+ PDF) - companion note to the tutorial: physics of the basis change, the three ℓ truncations, measured truncation errors, what the grid rerun needs. Build with `latexmk -pdf`.
- `Mathematica/benchmarks.md` - timings, precision scan, reference values, open questions. Append new measurements there.
- `Mathematica/findings.md` - derivations/convention findings (why, not timings); append new ones there.

## Where to change X
- Flux grid (`FluxesGrid.nb`) mode range or grid size: `nR`, `nA`, `lmax` in the definitions cell (`modes` derives from `lmax`).
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
- Never evaluate `SpheroidalToSphericalTutorial.nb` whole headless: it filled 15 GB RAM and froze the machine (2026-11-19; culprit cell not isolated). Run selected cells under `systemd-run --user --scope -p MemoryMax=3G -p MemorySwapMax=0` plus `MemoryConstrained`/`TimeConstrained` per cell.
