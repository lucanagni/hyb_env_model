# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Goal

Hybrid GSF+PN waveform model based on Honet+2025 (`biblio/Honet et al. - 2025 - ...pdf`), later extended with environmental effects. Current stage: reproduce Honet+ results and build a vacuum hybrid model. Currently reproducing the interpolated flux grid of Honet+ Appendix A (eqs. A13/A14 are referenced in the notebook).

## Layout

- `biblio/`: reference papers (PDFs). Read the Honet+ PDF for equations/definitions rather than guessing.
- `Mathematica/FluxGrid_nonspinning.nb`: the only code so far. Non-spinning (Schwarzschild) flux grid reproduction. Loads `KerrGeodesics`` and `Teukolsky`` (Black Hole Perturbation Toolkit) and computes energy fluxes at infinity (`"\[ScriptCapitalI]"`) and horizon (`"\[ScriptCapitalH]"`) via `TeukolskyPointParticleMode`, then builds the grid following A13/A14. Written with Mathematica 14.0 on Linux.

## Working with the notebook

- `.nb` files are large Wolfram expression text with cached outputs; don't read whole. Use `grep -a` for the `Cell[` / `RowBox` input cells, or ask the user for the relevant cell.
- Run headless with `wolframscript -file`/`-code` (needs the BHPToolkit packages installed); the notebook has a section "Tests and debugging (ignore)" that can be skipped.
- Settings: `lmax = 8`, `prec = 60` (prec 30 fails for r0 ≳ 2500 and near-extremal spin). The grid is built with `ParallelTable` over 4 kernels (`DistributeDefinitions` needed). Full grid ≈ 52 min at lmax=8 (≈ 11 min at lmax=4). Node (i=0, j=0) (ISCO, χ≈0.998) returns `Indeterminate` at lmax=8 (`SpinWeightedSpheroidalHarmonicS::maxterms`, 200-term Leaver limit); unresolved.
- Known issue: some modes at certain (p, e) emit `Infinite expression`/`Indeterminate expression` messages (1/0 warnings from the Teukolsky mode computation).
