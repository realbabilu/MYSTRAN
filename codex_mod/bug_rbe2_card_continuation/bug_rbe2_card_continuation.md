# RBE2 Continuation Truncation

This package captures the proven `RBE2` card-continuation bug from deck
generation and explains why it is not an active parser bug in the local
MYSTRAN source.

## Problem

The generated benchmark decks originally wrote too many slave grids on a single free-field `RBE2` line, for example:

```text
RBE2,2801,28,126,10,11,12,13,14,15,16,17,18
```

That format is not equivalent to "keep reading dependent grids until end of
line" in MYSTRAN.

The MYSTRAN bulk-data reader only accepts:

- first line: dependent grids in fields `5-9`
- continuation lines: dependent grids in fields `2-9`

Source proof:

- [BD_RBE2.f90](C:/mystran3/MYSTRANSolver-18.0.0/Source/LK1/L1A-BD/BD_RBE2.f90)

Relevant source snippet for this issue:

```fortran
! BD_RBE2 field contract:
!   first line:        fields 5-9 -> dependent grids
!   continuation line: fields 2-9 -> dependent grids
```

See:

- [bd_rbe2_continuation_snippet.f90](C:/mystran3/codex_mod/bug_rbe2_card_continuation/bd_rbe2_continuation_snippet.f90)
- [BD_RBE2.f90](C:/mystran3/codex_mod/bug_rbe2_card_continuation/BD_RBE2.f90)

## Wrong vs Correct Deck Form

Wrong single-line form:

```text
RBE2,2801,28,126,10,11,12,13,14,15,16,17,18
```

Correct continued form:

```text
RBE2,2801,28,126,10,11,12,13,14
,15,16,17,18
```

For the full floor diaphragms in `Problem 1-024`, continuation cards are mandatory because each diaphragm has nine slave grids.

## Local source status

The important distinction for this package is:

- the local `BD_RBE2` parser already handles continuation correctly
- the benchmark bug was caused by emitting too many dependent grids on one
  physical parent line
- so this item is a **deck-generation/input-format bug**, not a new source-code
  parser patch like the `SPC1` hardening case

## Proof It Truncated

The numerical proof is stronger than the syntax proof:

- before fixing continuation, `RBE2` and explicit `MPC` gave different static and modal results
- after fixing continuation, `RBE2 126` and explicit diaphragm `MPC` became numerically identical

The DOF accounting also changed in the expected direction:

- incorrect deck path previously showed `M-set = 30`
- corrected deck path shows `M-set = 54`

That is exactly what should happen for `18` slave grids with `3` constrained components each.

## Numerical Evidence

After fixing continuation:

- `RBE2 126`: `0.274964, 0.268458, 0.087143, 0.085224 s`
- explicit `MPC`: `0.274964, 0.268458, 0.087143, 0.085224 s`

Static comparison after fixing continuation:

- joint `28`, `Ux = 1.234744000e-04`, `Rz = -1.176715000e-06` for both `RBE2` and `MPC`
- joint `29`, `Ux = 2.265431000e-04`, `Rz = -1.136647000e-06` for both `RBE2` and `MPC`

## Conclusion

The proven `RBE2` bug for this benchmark was not a rigid-body math error in the
active source path.

It was a deck-generation bug:

- the emitted `RBE2` cards overflowed the allowed dependent-grid fields
- MYSTRAN then read only part of the intended slave list
- that truncation contaminated earlier conclusions about diaphragm behavior

So the correct local takeaway is:

- no new parser patch is required in `BD_RBE2`
- deck generators must emit continuation cards whenever the dependent-grid list
  exceeds fields `5-9` on the parent line

## Included Files

- [bd_rbe2_continuation_snippet.f90](C:/mystran3/codex_mod/bug_rbe2_card_continuation/bd_rbe2_continuation_snippet.f90)
- [bad_vs_good_rbe2_cards.txt](C:/mystran3/codex_mod/bug_rbe2_card_continuation/bad_vs_good_rbe2_cards.txt)
- [BD_RBE2.f90](C:/mystran3/codex_mod/bug_rbe2_card_continuation/BD_RBE2.f90)
