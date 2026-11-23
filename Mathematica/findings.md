# Findings

## SF vs PN mode-phase convention (2026-11-23)

`hybSF` in `HybridModes.nb` uses `δφ_lm = -1/m Arg[(-1)^(m+1) C_lm]`, not Honet eq. 54 (`-1/m Arg C_lm`).
Eq. 54 is correct only in Honet's sign convention; with BHPToolkit C^up it puts every even-m mode off by π in Arg.

Check: Arg[C_lm]/π from `ModeAmp` at r0 = 200, 400 (χ = 0, 0.5) vs leading PN Ĥ_lm (Blanchet/PNpedia):

| lm | SF raw | SF × (−1)^(m+1) | leading PN |
|----|--------|-----------------|------------|
| 22 | 1    | 0    | 0    |
| 21 | 1/2  | 1/2  | 1/2  |
| 33 | −1/2 | −1/2 | −1/2 |
| 32 | 1    | 0    | 0    |
| 31 | 1/2  | 1/2  | 1/2  |
| 44 | 0    | 1    | 1    |
| 43 | −1/2 | −1/2 | −1/2 |
| 55 | 1/2  | 1/2  | 1/2  |

(deviations from the listed values are O(1/r0), ≲ 0.006.)

The factor (−1)^(m+1) has two parts:
- **The overall −1** comes from the sign convention in the h–C relation: Honet eq. 52 has h_lm = +2 C_lm/(mΩ)², BHPToolkit's C^up (ψ₄ = ½ḧ, modes e^{−iωt}) gives h_lm = −2 C_lm/ω².
- **The (−1)^m** aligns the orbital-phase origin, which fixes the sign of the odd-m modes. These are the modes proportional to δ = (m₁−m₂)/M, so it amounts to choosing which body sits at φ = 0.

The (−1)^m does not cancel in δψ_lm = δφ_lm − δφ_22: δφ_22 is defined only mod π (h_22 cannot tell ψ from ψ+π),
so it is unchanged, while δφ_lm moves by π mod 2π/m — nothing for even m, π/m for odd m (flips the sign of e^{−imδψ_lm}).
Hence `Arg[-cSF]` alone is wrong for odd m; the factor is needed both in the hybrid and in direct PN/SF comparisons.
Without any factor, δψ_lm is off by π/2 − π/m (even m) or π/2 (odd m).
