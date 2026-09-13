# An exact certificate for the Heilbronn-type bound

*One page; the statement in the README (§3(d)) summarizes it. Found during an adversarial read on 2026-09-12 and checked by hand; the record is `RESULT_heilbronn_certificate_statement_errata_citation_cards_20260912.md` in the archive (DOI 10.5281/zenodo.22727504).*

## Setting

Let ρ be the conductor-N lift (N = 1951 or 2141) with determinant χ of order 5, K̃ the field it cuts out (Gal(K̃/Q) ≅ SL₂(F₅) × C₅, order 600), and K̃₀ ⊂ K̃ the degree-120 field cut out by the determinant-one twist ρ ⊗ χ², whose image is G = SL₂(F₅) with irreducible characters 1, 2, 2′, 3, 3′, 4, 4′, 5, 6 (the dimension names the character; "2" is ρ ⊗ χ² restricted to G).

For s₀ ≠ 1 the Heilbronn character θ = Σ_π n_π π, n_π = ord_{s₀} L(s, π, K̃₀/Q), has θ(1) = ord_{s₀} ζ_{K̃₀}. For a subgroup H and an irreducible character τ of H whose Artin L-function over K̃₀^H is entire, ⟨θ, Ind_H τ⟩ = ord_{s₀} L(s, τ, K̃₀/K̃₀^H) ≥ 0 (Foote–Murty, Math. Proc. Cambridge Philos. Soc. 105 (1989) 5–11).

## Five rows

Induced characters, with decompositions as printed by `scripts/heilbronn_coin.py` (output in `results/heilbronn_coin_output.txt`):

| row | decomposition | dim | entire because |
|---|---|---|---|
| Ind_{C₅} λ, λ of order 5 | 2 + 3′ + 4 + 4′ + 5 + 6 | 24 | Hecke (monomial) |
| Ind_{C₁₀} 1 | 1 + 3 + 3′ + 5 | 12 | Hecke |
| Ind_{C₁₀} ξ, ξ faithful of order 10 | 2 + 4′ + 6 | 12 | Hecke |
| Ind_{C₁₀} ξ′, ξ′ of order 5 | 3 + 4 + 5 | 12 | Hecke |
| Ind_{SL₂(F₃)} τ, τ a faithful 2-dim character of 2.A₄ | 2 + 2′ + 6 | 10 | Langlands (tetrahedral), over the quintic field K̃₀^{SL₂(F₃)} |

## The identity

With reg = 1 + 2·2 + 2·2′ + 3·3 + 3·3′ + 4·4 + 4·4′ + 5·5 + 6·6 the regular character,

  2·Ind_{C₅}λ + Ind_{C₁₀}1 + 2·Ind_{C₁₀}ξ + 2·Ind_{C₁₀}ξ′ + 2·Ind_{SL₂(F₃)}τ = reg + 4·[2].

Tally by character: [1]: 1 · [2]: 2+2+2 = 6 = 2+4 · [2′]: 2 · [3]: 1+2 = 3 · [3′]: 2+1 = 3 · [4]: 2+2 = 4 · [4′]: 2+2 = 4 · [5]: 2+1+2 = 5 · [6]: 2+2+2 = 6. ✓

## Consequences

1. Pairing with θ and using ⟨θ, reg⟩ = θ(1): **θ(1) + 4·n₂ ≥ 0**. A pole of L(s, 2) of order m at s₀ (n₂ = −m) forces **ord_{s₀} ζ_{K̃₀} ≥ 4m**. Without the tetrahedral row the linear program's minimum drops to 2, so Langlands' theorem over the quintic field is load-bearing.

2. Over the full group G × C₅: each row tensored with a character ψʲ of the C₅ factor is again a known-entire induced row (monomial, or tetrahedral over the same quintic field), so the identity holds sector by sector: M_j + 4·n_{2⊠ψʲ} ≥ 0, where M_j = ⟨θ, reg ⊠ ψʲ⟩ is the mass of sector j. Each M_j ≥ 0, since reg ⊠ ψʲ = Ind_{1×C₅}(1 ⊠ ψʲ) is a monomial row (a Hecke L-function of K̃₀). With ρ = 2 ⊠ ψ^{j₀}, j₀ ≠ 0: **ord_{s₀} ζ_{K̃} = Σ_j M_j ≥ M_{j₀} ≥ 4m.**

3. For **real** s₀, the conjugate ρ̄ = 2 ⊠ ψ^{−j₀} (the characters 2, 2′ are real) has a pole of the same order in the sector −j₀ ≠ j₀, so **ord_{s₀} ζ_{K̃} ≥ M_{j₀} + M_{−j₀} ≥ 8m.**

Hence ord_{s₀} ζ_{K̃} ≤ 3 (≤ 7 for real s₀) at s₀ ≠ 1 implies L(s, ρ) holomorphic there, for poles of any order. Stark's theorem (Invent. Math. 23 (1974) 135–152) gives this from ord ≤ 1 for any Galois extension; Booker (Experiment. Math. 15 (2006) 385–407, §2) observed that the monomial rows alone do not force n₂ ≥ 0 and named the configurations −[2] + [6], −[2′] + [6], −[4′] + [6].

## What remains computational

The *list* of minimal configurations at m = 1 — exactly three over SL₂(F₅): {2: −1, 6: +1}, {2: −1, 2′: +1, 4′: +1}, {2: −1, 2′: −1, 4′: −1, 6: +2} (`heilbronn_coin.py`, exhaustive at total order ≤ 4; replicated by `heilbronn_lp_replication.py`, a separate implementation in the same model family) — and the statement that the same three return over the full group (single implementation, 228 rows, not included here).

*Row counting.* The statement says "33 distinct rows (32 monomial, one tetrahedral)" following the 2026-09-04 record; `heilbronn_lp_replication.py` prints "31 monomial, 2 tetrahedral" — the same 33 rows, split differently because one twisted tetrahedral induction coincides with a monomial row. The bound does not depend on the split.
